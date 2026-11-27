#!/usr/bin/env python3
"""Extract NetMHCpan 4.1b PDF outputs and evaluate the exact BA holdout targets.

Run from a cloned project after placing:
- eight server-result PDFs in data/raw/reference_predictions/netmhcpan41_ba/
- each batch's target_pairs.csv in data/interim/netmhcpan41_ba_submissions/<job_name>/

Dependencies:
    python -m pip install pymupdf scikit-learn

The script keeps the full server cross-product predictions separate from the
exact benchmark-target prediction table. It uses -%Rank_BA as the continuous
score because lower BA percentile ranks indicate stronger predicted binding.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

try:
    import fitz  # PyMuPDF
except ImportError as exc:
    raise SystemExit("Missing PyMuPDF. Install with: python -m pip install pymupdf") from exc

try:
    from sklearn.metrics import (
        average_precision_score,
        confusion_matrix,
        f1_score,
        matthews_corrcoef,
        roc_auc_score,
    )
except ImportError as exc:
    raise SystemExit(
        "Missing scikit-learn. Install with: python -m pip install scikit-learn"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
PDF_DIR = ROOT / "data/raw/reference_predictions/netmhcpan41_ba"
TARGET_ROOT = ROOT / "data/interim/netmhcpan41_ba_submissions"
RAW_OUT = ROOT / "data/processed/reference_predictions"
TARGET_OUT = ROOT / "data/processed/binding"
RESULTS = ROOT / "results/tables"

ROW_RE = re.compile(r"^\s*(\d+)\s+(HLA-[A-Z]\*\d{2}:\d{2})\s+([A-Z]+)\s+")
HEADER_RE = re.compile(
    r"^\s*(HLA-[ABC]\d{2}:\d{2})\s*:\s*Distance to training data\s+"
    r"([0-9.]+)\s*\(using nearest neighbor\s+(HLA-[ABC]\d{2}:\d{2})\)"
)
JOBID_RE = re.compile(r"jobid=([0-9A-Fa-f]+)")
MIN_ACCEPTED_PEPTIDE_LENGTH = 8
STRONG_RANK_CUTOFF = 0.5
WEAK_RANK_CUTOFF = 2.0


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_target_pairs(job_name: str) -> list[dict]:
    path = TARGET_ROOT / job_name / "target_pairs.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing target-pair file: {path}")
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"hla_allele", "peptide", "binary_label"}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(f"{path} must contain {sorted(required)}")
        rows = list(reader)
    for row in rows:
        row["hla_allele"] = row["hla_allele"].strip()
        row["peptide"] = row["peptide"].strip().upper()
        row["binary_label"] = int(row["binary_label"])
        if row["binary_label"] not in (0, 1):
            raise RuntimeError(f"Invalid binary label in {path}: {row}")
    return rows


def parse_pdf(path: Path) -> tuple[list[dict], dict, str, int]:
    records: list[dict] = []
    allele_headers: dict[str, dict] = {}
    versions: set[str] = set()
    modes: set[str] = set()
    job_id = ""
    doc = fitz.open(path)
    for page_number, page in enumerate(doc, start=1):
        text = page.get_text("text")
        if not job_id:
            match = JOBID_RE.search(text)
            if match:
                job_id = match.group(1)
        for line in text.splitlines():
            if "NetMHCpan version" in line:
                versions.add(line.split("NetMHCpan version", 1)[1].strip())
            if "Make both EL and BA predictions" in line:
                modes.add("both EL and BA")
            header_match = HEADER_RE.match(line)
            if header_match:
                allele_name, distance, nearest = header_match.groups()
                allele_headers[allele_name] = {
                    "distance_to_training_data": float(distance),
                    "nearest_neighbor_allele": nearest,
                }
            if ROW_RE.match(line):
                tokens = line.split()
                if len(tokens) < 16:
                    raise RuntimeError(
                        f"Could not parse prediction row in {path.name}, page {page_number}: {line}"
                    )
                try:
                    records.append({
                        "position": int(tokens[0]),
                        "hla_allele": tokens[1],
                        "peptide": tokens[2],
                        "core": tokens[3],
                        "score_el": float(tokens[11]),
                        "rank_el_percent": float(tokens[12]),
                        "score_ba": float(tokens[13]),
                        "rank_ba_percent": float(tokens[14]),
                        "affinity_nm": float(tokens[15]),
                        "server_bindlevel": " ".join(tokens[16:]) if len(tokens) > 16 else "",
                        "source_job": path.stem,
                        "source_pdf": path.name,
                        "source_pdf_page": page_number,
                    })
                except (ValueError, IndexError) as exc:
                    raise RuntimeError(
                        f"Bad prediction row in {path.name}, page {page_number}: {line}"
                    ) from exc
    if not records:
        raise RuntimeError(f"No prediction rows found in {path}")
    if not versions:
        raise RuntimeError(f"No NetMHCpan version header found in {path}")
    if not allele_headers:
        raise RuntimeError(f"No allele nearest-neighbour headers found in {path}")
    for row in records:
        server_allele = row["hla_allele"].replace("*", "")
        if server_allele not in allele_headers:
            raise RuntimeError(
                f"No nearest-neighbour header for {server_allele} in {path.name}"
            )
        row.update(allele_headers[server_allele])
    meta = {
        "versions": sorted(versions),
        "modes": sorted(modes),
        "allele_job_sections": len(allele_headers),
        "distinct_alleles": len({r["hla_allele"] for r in records}),
    }
    return records, meta, job_id, len(doc)


def category_from_rank(rank: float) -> str:
    if rank <= STRONG_RANK_CUTOFF:
        return "strong"
    if rank <= WEAK_RANK_CUTOFF:
        return "weak"
    return "non_binder"


def binary_from_server_level(level: str) -> str:
    level = level.strip()
    if "<= SB" in level:
        return "SB"
    if "<= WB" in level:
        return "WB"
    return "NB"


def metric_record(rows: list[dict], group_name: str, group_value: str) -> dict:
    y = [int(r["binary_label"]) for r in rows]
    rank = [float(r["rank_ba_percent"]) for r in rows]
    n = len(rows)
    positive_n = sum(y)
    negative_n = n - positive_n
    result = {
        group_name: group_value,
        "n_target_pairs": n,
        "positive_n": positive_n,
        "negative_n": negative_n,
        "predicted_n": sum(r["rank_ba_percent"] is not None for r in rows),
        "prediction_coverage_percent": 100.0,
        "auc_roc": "",
        "average_precision": "",
        "auc_status": "estimable" if positive_n and negative_n else "not_estimable_single_class",
        "median_rank_ba_percent": sorted(rank)[n // 2] if n % 2 else (
            sorted(rank)[n // 2 - 1] + sorted(rank)[n // 2]
        ) / 2,
        "median_affinity_nm": sorted(float(r["affinity_nm"]) for r in rows)[n // 2],
    }
    if positive_n and negative_n:
        result["auc_roc"] = roc_auc_score(y, [-x for x in rank])
        result["average_precision"] = average_precision_score(y, [-x for x in rank])

    for cutoff, prefix in ((0.5, "strong_0_5"), (2.0, "binder_2")):
        pred = [int(x <= cutoff) for x in rank]
        tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
        result[f"{prefix}_predicted_n"] = int(sum(pred))
        result[f"{prefix}_tp"] = int(tp)
        result[f"{prefix}_fp"] = int(fp)
        result[f"{prefix}_tn"] = int(tn)
        result[f"{prefix}_fn"] = int(fn)
        result[f"{prefix}_sensitivity"] = tp / (tp + fn) if tp + fn else ""
        result[f"{prefix}_specificity"] = tn / (tn + fp) if tn + fp else ""
        result[f"{prefix}_balanced_accuracy"] = (
            (tp / (tp + fn) + tn / (tn + fp)) / 2
            if tp + fn and tn + fp else ""
        )
        result[f"{prefix}_precision"] = tp / (tp + fp) if tp + fp else ""
        result[f"{prefix}_f1"] = f1_score(y, pred, zero_division=0) if positive_n else ""
        result[f"{prefix}_mcc"] = (
            matthews_corrcoef(y, pred) if positive_n and negative_n else ""
        )
    return result


def main() -> int:
    pdf_paths = sorted(PDF_DIR.glob("job_*.pdf"))
    if len(pdf_paths) != 8:
        raise RuntimeError(f"Expected 8 batch PDFs in {PDF_DIR}; found {len(pdf_paths)}")

    all_predictions: list[dict] = []
    target_predictions: list[dict] = []
    job_qc: list[dict] = []
    versions_seen: set[str] = set()
    job_modes: list[list[str]] = []
    failed_pairs: list[dict] = []
    target_keys_seen: set[tuple[str, str]] = set()
    source_headers_by_allele: dict[str, dict] = {}

    for pdf_path in pdf_paths:
        job_name = pdf_path.stem
        target_rows = read_target_pairs(job_name)
        pdf_rows, meta, job_id, pages = parse_pdf(pdf_path)
        versions_seen.update(meta["versions"])
        job_modes.append(meta["modes"])

        prediction_map: dict[tuple[str, str], dict] = {}
        for pred in pdf_rows:
            key = (pred["hla_allele"], pred["peptide"])
            if key in prediction_map:
                raise RuntimeError(f"Duplicate prediction pair in {job_name}: {key}")
            prediction_map[key] = pred
        all_predictions.extend(pdf_rows)

        matched = 0
        for target in target_rows:
            key = (target["hla_allele"], target["peptide"])
            if key in target_keys_seen:
                raise RuntimeError(f"Duplicate target pair across batches: {key}")
            target_keys_seen.add(key)
            prediction = prediction_map.get(key)
            if prediction is None:
                failed_pairs.append({
                    "source_job": job_name,
                    "hla_allele": target["hla_allele"],
                    "peptide": target["peptide"],
                    "binary_label": target["binary_label"],
                    "failure_reason": "target_pair_missing_from_server_pdf",
                })
                continue
            matched += 1
            record = {
                "source_job": job_name,
                "source_pdf": pdf_path.name,
                "source_pdf_page": prediction["source_pdf_page"],
                "hla_allele": target["hla_allele"],
                "locus": target["hla_allele"].split("-")[1][0],
                "peptide": target["peptide"],
                "peptide_length": len(target["peptide"]),
                "binary_label": target["binary_label"],
                "score_ba": prediction["score_ba"],
                "rank_ba_percent": prediction["rank_ba_percent"],
                "affinity_nm": prediction["affinity_nm"],
                "score_el": prediction["score_el"],
                "rank_el_percent": prediction["rank_el_percent"],
                "server_bindlevel": prediction["server_bindlevel"],
                "ba_category_0_5_2": category_from_rank(prediction["rank_ba_percent"]),
                "nearest_neighbor_allele": prediction["nearest_neighbor_allele"],
                "distance_to_training_data": prediction["distance_to_training_data"],
            }
            target_predictions.append(record)

        job_qc.append({
            "batch": job_name,
            "job_id": job_id,
            "pdf_file": pdf_path.name,
            "pdf_pages": pages,
            "netmhcpan_version": ";".join(meta["versions"]),
            "prediction_mode": ";".join(meta["modes"]),
            "allele_job_sections": meta["allele_job_sections"],
            "distinct_alleles": meta["distinct_alleles"],
            "server_prediction_rows": len(pdf_rows),
            "target_pairs": len(target_rows),
            "target_pairs_matched": matched,
            "target_pairs_missing": len(target_rows) - matched,
            "target_pair_coverage_percent": 100 * matched / len(target_rows) if target_rows else "",
        })

    if len(versions_seen) != 1 or "4.1b" not in versions_seen:
        raise RuntimeError(f"Unexpected/mixed NetMHCpan versions: {sorted(versions_seen)}")
    if len(failed_pairs):
        raise RuntimeError(f"{len(failed_pairs)} target predictions are missing; see failure report.")

    # Build the full cross-product dataframe and target table outputs.
    raw_columns = [
        "source_job", "source_pdf", "source_pdf_page", "hla_allele", "peptide",
        "core", "position", "score_el", "rank_el_percent", "score_ba",
        "rank_ba_percent", "affinity_nm", "server_bindlevel",
        "nearest_neighbor_allele", "distance_to_training_data",
    ]
    target_columns = [
        "source_job", "source_pdf", "source_pdf_page", "hla_allele", "locus",
        "peptide", "peptide_length", "binary_label", "score_ba",
        "rank_ba_percent", "affinity_nm", "score_el", "rank_el_percent",
        "server_bindlevel", "ba_category_0_5_2", "nearest_neighbor_allele",
        "distance_to_training_data",
    ]
    RAW_OUT.mkdir(parents=True, exist_ok=True)
    TARGET_OUT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    raw_path = RAW_OUT / "netmhcpan41_ba_el_all_server_pairs.csv.gz"
    with gzip.open(raw_path, "wt", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=raw_columns, extrasaction="ignore")
        writer.writeheader()
        for row in all_predictions:
            writer.writerow(row)

    target_path = TARGET_OUT / "netmhcpan41_ba_target_pair_predictions.csv"
    write_csv(target_path, target_predictions, target_columns)

    allele_to_rows: dict[str, list[dict]] = defaultdict(list)
    locus_to_rows: dict[str, list[dict]] = defaultdict(list)
    for row in target_predictions:
        allele_to_rows[row["hla_allele"]].append(row)
        locus_to_rows[row["locus"]].append(row)

    allele_metrics = []
    for allele in sorted(allele_to_rows):
        group_rows = allele_to_rows[allele]
        metric = metric_record(group_rows, "hla_allele", allele)
        metric["locus"] = group_rows[0]["locus"]
        metric["nearest_neighbor_allele"] = group_rows[0]["nearest_neighbor_allele"]
        metric["distance_to_training_data"] = group_rows[0]["distance_to_training_data"]
        metric["nearest_neighbor_is_exact_allele"] = (
            allele.replace("*", "") == group_rows[0]["nearest_neighbor_allele"]
        )
        metric["sample_size_note"] = (
            "n<20; descriptive only" if len(group_rows) < 20
            else "n<50; interpret cautiously" if len(group_rows) < 50
            else "n>=50"
        )
        allele_metrics.append(metric)
    write_csv(RESULTS / "netmhcpan41_ba_per_allele_metrics.csv", allele_metrics)

    locus_metrics = []
    for locus in ("A", "B", "C"):
        group_rows = locus_to_rows.get(locus, [])
        if not group_rows:
            continue
        metric = metric_record(group_rows, "locus", f"HLA-{locus}")
        metric["alleles_in_locus"] = len({r["hla_allele"] for r in group_rows})
        metric["unique_peptides_in_locus"] = len({r["peptide"] for r in group_rows})
        locus_metrics.append(metric)
    write_csv(RESULTS / "netmhcpan41_ba_per_locus_metrics.csv", locus_metrics)

    overall = metric_record(target_predictions, "scope", "All eligible target pairs")
    overall["alleles"] = len(allele_to_rows)
    overall["unique_peptides"] = len({r["peptide"] for r in target_predictions})
    write_csv(RESULTS / "netmhcpan41_ba_overall_metrics.csv", [overall])

    mapping = []
    for allele in sorted(allele_to_rows):
        row = allele_to_rows[allele][0]
        mapping.append({
            "hla_allele": allele,
            "nearest_neighbor_allele": row["nearest_neighbor_allele"],
            "distance_to_training_data": row["distance_to_training_data"],
            "nearest_neighbor_is_exact_allele": (
                allele.replace("*", "") == row["nearest_neighbor_allele"]
            ),
            "target_pair_rows": len(allele_to_rows[allele]),
            "predicted_target_rows": len(allele_to_rows[allele]),
            "coverage_percent": 100.0,
        })
    write_csv(RESULTS / "netmhcpan41_ba_allele_nearest_neighbor_mapping.csv", mapping)

    write_csv(RESULTS / "netmhcpan41_ba_per_job_qc.csv", job_qc)
    write_csv(
        RESULTS / "netmhcpan41_ba_failed_predictions.csv",
        failed_pairs,
        ["source_job", "hla_allele", "peptide", "binary_label", "failure_reason"],
    )
    # Explicit empty report: all target alleles had returned prediction sections.
    write_csv(
        RESULTS / "netmhcpan41_ba_unsupported_alleles.csv",
        [],
        ["hla_allele", "locus", "reason"],
    )
    write_csv(
        RESULTS / "netmhcpan41_ba_length_eligibility_summary.csv",
        [{
            "original_reference_disjoint_holdout_rows": 3299,
            "scored_rows": len(target_predictions),
            "length_ineligible_rows": 3299 - len(target_predictions),
            "scored_fraction_of_original_percent": (
                100 * len(target_predictions) / 3299
            ),
            "note": (
                "Rows below the DTU minimum peptide length were excluded before "
                "submission; the original holdout was not modified."
            ),
        }],
        [
            "original_reference_disjoint_holdout_rows",
            "scored_rows",
            "length_ineligible_rows",
            "scored_fraction_of_original_percent",
            "note",
        ],
    )

    server_key_counts = Counter(
        (row["hla_allele"], row["peptide"]) for row in all_predictions
    )
    if any(count > 1 for count in server_key_counts.values()):
        raise RuntimeError("Duplicate HLA-peptide predictions found across PDF outputs.")

    bindlevel_el_agree = 0
    bindlevel_ba_agree = 0
    for row in target_predictions:
        level = row["server_bindlevel"].strip()
        level = "SB" if "<= SB" in level else "WB" if "<= WB" in level else "NB"
        el_rank = row["rank_el_percent"]
        el_class = "SB" if el_rank <= 0.5 else "WB" if el_rank <= 2.0 else "NB"
        ba_rank = row["rank_ba_percent"]
        ba_class = "SB" if ba_rank <= 0.5 else "WB" if ba_rank <= 2.0 else "NB"
        bindlevel_el_agree += level == el_class
        bindlevel_ba_agree += level == ba_class

    versions = sorted(versions_seen)
    unique_alleles = len(allele_to_rows)
    nearest_diff = sum(not row["nearest_neighbor_is_exact_allele"] for row in mapping)
    nonzero_distance = sum(float(row["distance_to_training_data"]) > 0 for row in mapping)
    auc_estimable = sum(row["auc_status"] == "estimable" for row in allele_metrics)

    qc = [
        {"check": "PDF files processed", "observed": len(pdf_paths), "status": "PASS" if len(pdf_paths)==8 else "FAIL"},
        {"check": "NetMHCpan version", "observed": ";".join(versions), "status": "PASS"},
        {"check": "Server prediction rows parsed", "observed": len(all_predictions), "status": "PASS"},
        {"check": "Unique HLA-peptide server prediction pairs", "observed": len(server_key_counts), "status": "PASS" if len(server_key_counts)==len(all_predictions) else "FAIL"},
        {"check": "Target peptide-allele pairs", "observed": len(target_predictions), "status": "PASS"},
        {"check": "Target prediction coverage (%)", "observed": 100.0, "status": "PASS"},
        {"check": "Missing target predictions", "observed": len(failed_pairs), "status": "PASS" if not failed_pairs else "FAIL"},
        {"check": "Distinct target alleles", "observed": unique_alleles, "status": "PASS"},
        {"check": "Alleles with nearest neighbour different from requested allele", "observed": nearest_diff, "status": "INFO"},
        {"check": "Alleles with non-zero training distance", "observed": nonzero_distance, "status": "INFO"},
        {"check": "Per-allele AUC estimable (both classes present)", "observed": auc_estimable, "status": "INFO"},
        {"check": "Per-allele AUC not estimable (single class)", "observed": unique_alleles-auc_estimable, "status": "INFO"},
        {"check": "Server BindLevel agrees with EL rank categories (%)", "observed": round(100*bindlevel_el_agree/len(target_predictions),3), "status": "INFO"},
        {"check": "Server BindLevel agrees with BA rank categories (%)", "observed": round(100*bindlevel_ba_agree/len(target_predictions),3), "status": "INFO"},
    ]
    write_csv(RESULTS / "netmhcpan41_ba_extraction_qc.csv", qc)

    manifest = {
        "project": "SparQ-HLA-Prediction",
        "predictor": "NetMHCpan 4.1b",
        "service": "DTU Health Tech NetMHCpan 4.1 web server",
        "input_type": "PEPTIDE",
        "prediction_modes": ["EL", "BA"],
        "primary_benchmark_mode": "BA",
        "primary_continuous_score": "-%Rank_BA (lower reported rank means stronger predicted binding)",
        "secondary_scores": ["Score_BA", "Aff(nM)"],
        "server_reported_strong_rank_threshold_percent": 0.5,
        "server_reported_weak_rank_threshold_percent": 2.0,
        "ba_threshold_classification": {
            "strong": "%Rank_BA <= 0.5",
            "weak": "0.5 < %Rank_BA <= 2.0",
            "non_binder": "%Rank_BA > 2.0"
        },
        "rank_filter_setting": (
            "The batch preparation instructions specified -99 (retain all). "
            "The PDF header itself does not echo this UI setting."
        ),
        "original_reference_disjoint_holdout_rows": 3299,
        "server_length_ineligible_rows": 14,
        "eligible_target_rows_scored": len(target_predictions),
        "unique_target_peptides": len({r["peptide"] for r in target_predictions}),
        "distinct_hla_alleles": unique_alleles,
        "server_prediction_cross_product_rows": len(all_predictions),
        "exact_target_pairs_matched": len(target_predictions),
        "target_pair_prediction_coverage_percent": 100.0,
        "failed_or_missing_target_predictions": len(failed_pairs),
        "alleles_nearest_neighbor_differs_from_requested": nearest_diff,
        "alleles_with_nonzero_distance_to_training_data": nonzero_distance,
        "per_allele_auc_estimable": auc_estimable,
        "per_allele_auc_not_estimable_single_class": unique_alleles - auc_estimable,
        "server_bindlevel_note": (
            "Across the 3,285 target rows, the server BindLevel label agrees with "
            "the EL rank thresholds, not the BA rank thresholds. For BA classification, "
            "derive strong/weak/non-binder labels from %Rank_BA explicitly."
        ),
        "model_files_note": (
            "Predictions were obtained from DTU's hosted 4.1b service. The hosted "
            "model parameter files are not included or directly inspectable in the PDF "
            "outputs. Repository NetMHCpan_train files are reference training data, "
            "not the prediction executable/model parameter package."
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "results": {
            "target_prediction_table": str(target_path.relative_to(ROOT)),
            "all_server_predictions": str(raw_path.relative_to(ROOT)),
            "per_allele_metrics": str((RESULTS / "netmhcpan41_ba_per_allele_metrics.csv").relative_to(ROOT)),
            "per_locus_metrics": str((RESULTS / "netmhcpan41_ba_per_locus_metrics.csv").relative_to(ROOT)),
            "overall_metrics": str((RESULTS / "netmhcpan41_ba_overall_metrics.csv").relative_to(ROOT)),
            "allele_mapping": str((RESULTS / "netmhcpan41_ba_allele_nearest_neighbor_mapping.csv").relative_to(ROOT)),
            "failed_predictions": str((RESULTS / "netmhcpan41_ba_failed_predictions.csv").relative_to(ROOT)),
            "qc": str((RESULTS / "netmhcpan41_ba_extraction_qc.csv").relative_to(ROOT)),
        }
    }
    (RESULTS / "netmhcpan41_ba_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print("\nNETMHCPAN 4.1b BA BENCHMARK CONSOLIDATION")
    print(f"PDF jobs processed: {len(pdf_paths)}")
    print(f"Server prediction rows parsed: {len(all_predictions):,}")
    print(f"Exact target rows matched: {len(target_predictions):,}")
    print(f"Prediction coverage: 100.0%")
    print(f"Failed/missing target predictions: {len(failed_pairs)}")
    print(f"Alleles with exact nearest-neighbour identity: {unique_alleles-nearest_diff}/{unique_alleles}")
    print(f"Alleles with alternate nearest neighbour: {nearest_diff}")
    print(f"Alleles with non-zero distance: {nonzero_distance}")
    print(f"Per-allele AUROC estimable: {auc_estimable}/{unique_alleles}")
    print(f"Results written under: {ROOT / 'results'}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
