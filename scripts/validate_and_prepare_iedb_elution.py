from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1] if Path(__file__).parent.name == "scripts" else Path.cwd()
RAW_ROOT = ROOT / "data" / "raw" / "processing" / "elution"
PROCESSED_DIR = ROOT / "data" / "processed" / "elution"
RESULTS_DIR = ROOT / "results" / "tables"
MANIFEST_PATH = RAW_ROOT / "manifest.json"

OUTPUT_PATH = PROCESSED_DIR / "iedb_hla_i_elution_observations.csv.gz"
TEMP_OUTPUT_PATH = PROCESSED_DIR / "iedb_hla_i_elution_observations.csv.gz.tmp"
SUMMARY_PATH = RESULTS_DIR / "iedb_elution_validation_summary.csv"
ERRORS_PATH = RESULTS_DIR / "iedb_elution_validation_errors.csv"

SOURCE_FIELDS = [
    "elution_id",
    "linear_sequence",
    "mhc_class",
    "mhc_allele_name",
    "mhc_allele_resolution",
    "assay_names",
    "qualitative_measure",
    "reference_id",
    "reference_dates",
    "parent_source_antigen_name",
    "parent_source_antigen_iri",
]
OUTPUT_FIELDS = [
    "elution_id",
    "linear_sequence",
    "peptide_length",
    "hla_locus",
    "mhc_class",
    "mhc_allele_name",
    "mhc_allele_resolution",
    "assay_names",
    "qualitative_measure",
    "reference_id",
    "reference_dates",
    "parent_source_antigen_name",
    "parent_source_antigen_iri",
]
MAX_ERROR_ROWS = 5000


def csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value)


def text_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False)


def id_sort_key(value: Any) -> tuple[int, Any]:
    value_text = str(value).strip()
    try:
        return (0, int(value_text))
    except (TypeError, ValueError):
        return (1, value_text)


def normalise_token(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


def add_error(
    errors: list[dict[str, str]],
    counters: Counter,
    locus: str,
    filename: str,
    line_number: str,
    category: str,
    detail: str,
) -> None:
    counters[category] += 1
    if len(errors) < MAX_ERROR_ROWS:
        errors.append(
            {
                "locus": locus,
                "file": filename,
                "line": line_number,
                "error": category,
                "detail": detail[:500],
            }
        )


def main() -> None:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Manifest not found: {MANIFEST_PATH}")

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    expected_by_locus = manifest.get("download_state", {}).get("next_offset", {})
    chunk_size = int(manifest.get("chunk_size", 10000))

    if chunk_size <= 0:
        raise ValueError(f"Invalid chunk_size in manifest: {chunk_size}")

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_OUTPUT_PATH.unlink(missing_ok=True)

    errors: list[dict[str, str]] = []
    summary_rows: list[dict[str, Any]] = []
    all_checks_passed = True

    totals = {
        "chunks": 0,
        "records": 0,
        "parse_errors": 0,
        "validation_errors": 0,
    }

    with gzip.open(TEMP_OUTPUT_PATH, "wt", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()

        for locus in ("A", "B", "C"):
            folder = RAW_ROOT / f"HLA_{locus}_2_chain"
            expected_records = int(expected_by_locus.get(locus, -1))
            expected_chunks = math.ceil(expected_records / chunk_size) if expected_records >= 0 else -1
            counters: Counter = Counter()
            previous_id_key = None
            previous_id_text = None
            observed_records = 0

            if expected_records < 0:
                add_error(
                    errors, counters, locus, str(MANIFEST_PATH.name), "", 
                    "missing_manifest_count", f"No expected record count for HLA-{locus}"
                )
                all_checks_passed = False

            if not folder.is_dir():
                add_error(
                    errors, counters, locus, str(folder), "", 
                    "missing_locus_directory", "Expected locus directory does not exist"
                )
                all_checks_passed = False
                chunk_files: list[Path] = []
            else:
                chunk_files = sorted(folder.glob("chunk_*.jsonl.gz"))

            expected_names = {
                f"chunk_{index:05d}.jsonl.gz" for index in range(max(expected_chunks, 0))
            }
            observed_names = {path.name for path in chunk_files}
            missing_names = sorted(expected_names - observed_names)
            extra_names = sorted(observed_names - expected_names) if expected_chunks >= 0 else []

            for missing_name in missing_names:
                add_error(
                    errors, counters, locus, missing_name, "", 
                    "missing_chunk", "Chunk required by manifest is missing"
                )
            for extra_name in extra_names:
                add_error(
                    errors, counters, locus, extra_name, "", 
                    "unexpected_chunk", "Chunk is not expected from manifest offsets"
                )

            if missing_names or extra_names:
                all_checks_passed = False

            for chunk_path in chunk_files:
                totals["chunks"] += 1
                chunk_records_before = observed_records

                try:
                    with gzip.open(chunk_path, "rt", encoding="utf-8") as chunk_file:
                        for line_number, line in enumerate(chunk_file, start=1):
                            if not line.strip():
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "blank_jsonl_line", "Blank line encountered"
                                )
                                continue

                            try:
                                record = json.loads(line)
                            except json.JSONDecodeError as exc:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "invalid_json", str(exc)
                                )
                                totals["parse_errors"] += 1
                                continue

                            if not isinstance(record, dict):
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "invalid_record_type", f"Expected JSON object, got {type(record).__name__}"
                                )
                                continue

                            observed_records += 1
                            totals["records"] += 1

                            missing_fields = [field for field in SOURCE_FIELDS if field not in record]
                            if missing_fields:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "missing_fields", ", ".join(missing_fields)
                                )

                            record_id = record.get("elution_id")
                            if record_id is None or str(record_id).strip() == "":
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "missing_elution_id", "elution_id is empty"
                                )
                            else:
                                current_id_key = id_sort_key(record_id)
                                if previous_id_key is not None and current_id_key <= previous_id_key:
                                    category = (
                                        "duplicate_elution_id"
                                        if current_id_key == previous_id_key
                                        else "non_increasing_elution_id"
                                    )
                                    add_error(
                                        errors, counters, locus, chunk_path.name, str(line_number),
                                        category,
                                        f"Current ID {record_id!r} follows {previous_id_text!r}; expected strictly increasing IDs"
                                    )
                                previous_id_key = current_id_key
                                previous_id_text = str(record_id)

                            sequence = str(record.get("linear_sequence") or "").strip()
                            if not sequence:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "missing_peptide_sequence", "linear_sequence is empty"
                                )

                            allele = str(record.get("mhc_allele_name") or "").strip()
                            if not allele.upper().startswith(f"HLA-{locus}"):
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "allele_locus_mismatch", f"Expected HLA-{locus} allele, got {allele!r}"
                                )

                            mhc_class = normalise_token(record.get("mhc_class"))
                            if mhc_class not in {"i", "mhci", "mhcclassi"}:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "unexpected_mhc_class", f"Expected class I, got {record.get('mhc_class')!r}"
                                )

                            resolution = normalise_token(record.get("mhc_allele_resolution"))
                            if resolution not in {"2chain", "twochain"}:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "unexpected_allele_resolution",
                                    f"Expected 2 chain, got {record.get('mhc_allele_resolution')!r}"
                                )

                            assay_text = text_value(record.get("assay_names")).lower()
                            if "mass spectrometry" not in assay_text:
                                add_error(
                                    errors, counters, locus, chunk_path.name, str(line_number),
                                    "assay_filter_mismatch",
                                    f"Assay does not mention mass spectrometry: {assay_text[:200]!r}"
                                )

                            output_row = {
                                "elution_id": csv_value(record.get("elution_id")),
                                "linear_sequence": sequence,
                                "peptide_length": len(sequence),
                                "hla_locus": locus,
                            }
                            for field in SOURCE_FIELDS[2:]:
                                output_row[field] = csv_value(record.get(field))
                            writer.writerow(output_row)

                except (OSError, EOFError, UnicodeDecodeError) as exc:
                    add_error(
                        errors, counters, locus, chunk_path.name, "", 
                        "chunk_read_error", f"{type(exc).__name__}: {exc}"
                    )
                    totals["parse_errors"] += 1

                if observed_records == chunk_records_before:
                    add_error(
                        errors, counters, locus, chunk_path.name, "", 
                        "empty_chunk", "No valid JSON objects were read from this chunk"
                    )

            if observed_records != expected_records:
                add_error(
                    errors, counters, locus, folder.name, "", 
                    "record_count_mismatch",
                    f"Expected {expected_records:,} records; observed {observed_records:,}"
                )

            locus_validation_errors = sum(
                count for category, count in counters.items()
                if category not in {"missing_chunk", "unexpected_chunk", "missing_locus_directory", "missing_manifest_count"}
            )
            if counters:
                all_checks_passed = False

            summary_row = {
                "locus": f"HLA-{locus}",
                "expected_records": expected_records,
                "observed_records": observed_records,
                "expected_chunks": expected_chunks,
                "observed_chunks": len(chunk_files),
                "missing_chunks": len(missing_names),
                "unexpected_chunks": len(extra_names),
                "validation_errors": sum(counters.values()),
            }
            for category, count in sorted(counters.items()):
                summary_row[category] = count
            summary_rows.append(summary_row)

            totals["validation_errors"] += sum(counters.values())
            print(
                f"HLA-{locus}: {observed_records:,}/{expected_records:,} records | "
                f"{len(chunk_files)}/{expected_chunks} chunks | "
                f"{sum(counters.values())} validation issues"
            )

    # Save reports regardless of whether the processed dataset passes validation.
    summary_fields = list(dict.fromkeys(key for row in summary_rows for key in row))
    with SUMMARY_PATH.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=summary_fields, restval="")
        writer.writeheader()
        writer.writerows(summary_rows)

    error_fields = ["locus", "file", "line", "error", "detail"]
    with ERRORS_PATH.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=error_fields)
        writer.writeheader()
        writer.writerows(errors)

    expected_total = sum(int(expected_by_locus.get(locus, 0)) for locus in ("A", "B", "C"))
    observed_total = sum(int(row["observed_records"]) for row in summary_rows)

    print("\n========== IEDB ELUTION VALIDATION ==========")
    print(f"Chunks read:       {totals['chunks']:,}")
    print(f"Records parsed:    {observed_total:,}")
    print(f"Expected records:  {expected_total:,}")
    print(f"Validation issues: {totals['validation_errors']:,}")
    print(f"Summary:           {SUMMARY_PATH}")
    print(f"Error details:     {ERRORS_PATH}")

    if all_checks_passed and observed_total == expected_total and totals["validation_errors"] == 0:
        TEMP_OUTPUT_PATH.replace(OUTPUT_PATH)
        print("\nPASS: All chunks validated; processed elution dataset assembled.")
        print(f"Dataset: {OUTPUT_PATH}")
    else:
        TEMP_OUTPUT_PATH.unlink(missing_ok=True)
        print("\nFAIL: Validation issues found. No final processed dataset was published.")
        print("Review the validation summary and error report before proceeding.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
