#!/usr/bin/env python3
"""Create sequence-grouped IEDB holdouts and audit them against NetMHCpan 4.1.

Run from the SparQ-HLA-Prediction project root:
    python scripts/create_holdouts_and_audit_netmhcpan.py

Uses only the Python standard library. Existing source datasets are not modified.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

ROOT = Path.cwd()
BINDING_INPUT = ROOT / "data/processed/binding/iedb_2chain_training_candidates.csv"
ELUTION_INPUT = ROOT / "data/processed/elution/iedb_hla_i_elution_observations.csv.gz"
REFERENCE_DIR = ROOT / "data/reference/netmhcpan41/NetMHCpan_train"
ALLELELIST = REFERENCE_DIR / "allelelist"

BINDING_DIR = ROOT / "data/processed/binding"
ELUTION_DIR = ROOT / "data/processed/elution"
TABLES_DIR = ROOT / "results/tables"

SEED = 42
HOLDOUT_PERCENT = 20
NETMHCPAN_CONFIGS = ["c000", "c001", "c002", "c003", "c004"]

# Accept HLA-A02:01, HLA-A*02:01, and HLA-A0201 forms.
HLA_RE = re.compile(r"^(?:HLA[-_ ]*)?([ABC])\*?(\d{2}):?(\d{2})(?::\d{2,3})?$", re.I)


def normalise_peptide(value: str | None) -> str:
    return (value or "").strip().upper().replace(" ", "")


def normalise_allele(value: str | None) -> str | None:
    if not value:
        return None
    value = value.strip().strip(",;\"'")
    match = HLA_RE.match(value)
    if not match:
        return None
    locus, field1, field2 = match.groups()
    return f"HLA-{locus.upper()}*{field1}:{field2}"


def is_holdout(peptide: str) -> bool:
    digest = hashlib.sha256(f"{SEED}:{peptide}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % 100 < HOLDOUT_PERCENT


def open_csv_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def open_csv_write(path: Path, compressed: bool = False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        return gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=6)
    return path.open("w", encoding="utf-8-sig", newline="")


def require_inputs() -> None:
    missing = [p for p in (BINDING_INPUT, ELUTION_INPUT, ALLELELIST) if not p.exists()]
    missing += [REFERENCE_DIR / f"{config}_{kind}" for config in NETMHCPAN_CONFIGS
                for kind in ("ba", "el") if not (REFERENCE_DIR / f"{config}_{kind}").exists()]
    if missing:
        print("Missing required input files:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        raise SystemExit(2)


def split_binding() -> Tuple[Dict[str, int], Dict[str, Set[str]]]:
    train_path = BINDING_DIR / "iedb_2chain_train.csv"
    holdout_path = BINDING_DIR / "iedb_2chain_holdout.csv"
    counts = {"train_rows": 0, "holdout_rows": 0}
    unique = {"train": set(), "holdout": set()}
    with open_csv_text(BINDING_INPUT) as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or not {"peptide", "hla_allele"}.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected binding columns in {BINDING_INPUT}: {reader.fieldnames}")
        with open_csv_write(train_path) as train_f, open_csv_write(holdout_path) as hold_f:
            train_writer = csv.DictWriter(train_f, fieldnames=reader.fieldnames)
            hold_writer = csv.DictWriter(hold_f, fieldnames=reader.fieldnames)
            train_writer.writeheader()
            hold_writer.writeheader()
            for row in reader:
                peptide = normalise_peptide(row.get("peptide"))
                if not peptide:
                    continue
                if is_holdout(peptide):
                    hold_writer.writerow(row)
                    counts["holdout_rows"] += 1
                    unique["holdout"].add(peptide)
                else:
                    train_writer.writerow(row)
                    counts["train_rows"] += 1
                    unique["train"].add(peptide)
    if unique["train"] & unique["holdout"]:
        raise RuntimeError("Binding split leakage: a peptide sequence occurs in both partitions")
    counts["train_unique_peptides"] = len(unique["train"])
    counts["holdout_unique_peptides"] = len(unique["holdout"])
    counts["total_unique_peptides"] = len(unique["train"] | unique["holdout"])
    print(f"Binding split: {counts['train_rows']:,} train / {counts['holdout_rows']:,} holdout rows")
    return counts, unique


def split_elution() -> Tuple[Dict[str, int], Dict[str, Set[str]], Dict[str, Set[str]]]:
    train_path = ELUTION_DIR / "iedb_hla_i_elution_train.csv.gz"
    holdout_path = ELUTION_DIR / "iedb_hla_i_elution_holdout.csv.gz"
    counts = {"train_rows": 0, "holdout_rows": 0}
    unique = {"train": set(), "holdout": set()}
    # Holdout sequence -> alleles observed in IEDB. This supports allele-specific auditing.
    holdout_by_sequence: Dict[str, Set[str]] = defaultdict(set)
    with open_csv_text(ELUTION_INPUT) as source:
        reader = csv.DictReader(source)
        needed = {"linear_sequence", "mhc_allele_name"}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected elution columns in {ELUTION_INPUT}: {reader.fieldnames}")
        with open_csv_write(train_path, compressed=True) as train_f, open_csv_write(holdout_path, compressed=True) as hold_f:
            train_writer = csv.DictWriter(train_f, fieldnames=reader.fieldnames)
            hold_writer = csv.DictWriter(hold_f, fieldnames=reader.fieldnames)
            train_writer.writeheader()
            hold_writer.writeheader()
            for row_number, row in enumerate(reader, 1):
                peptide = normalise_peptide(row.get("linear_sequence"))
                if not peptide:
                    continue
                allele = normalise_allele(row.get("mhc_allele_name"))
                if is_holdout(peptide):
                    hold_writer.writerow(row)
                    counts["holdout_rows"] += 1
                    unique["holdout"].add(peptide)
                    if allele:
                        holdout_by_sequence[peptide].add(allele)
                else:
                    train_writer.writerow(row)
                    counts["train_rows"] += 1
                    unique["train"].add(peptide)
                if row_number % 250000 == 0:
                    print(f"  Elution rows processed: {row_number:,}")
    if unique["train"] & unique["holdout"]:
        raise RuntimeError("Elution split leakage: a peptide sequence occurs in both partitions")
    counts["train_unique_peptides"] = len(unique["train"])
    counts["holdout_unique_peptides"] = len(unique["holdout"])
    counts["total_unique_peptides"] = len(unique["train"] | unique["holdout"])
    counts["holdout_unique_peptide_allele_pairs"] = sum(len(v) for v in holdout_by_sequence.values())
    print(f"Elution split: {counts['train_rows']:,} train / {counts['holdout_rows']:,} holdout rows")
    return counts, unique, holdout_by_sequence


def load_cell_line_alleles() -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with ALLELELIST.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            label = parts[0]
            allele_blob = " ".join(parts[1:])
            alleles = {a for token in allele_blob.split(",") if (a := normalise_allele(token))}
            if alleles:
                mapping[label] = alleles
    print(f"NetMHCpan cell-line mappings loaded: {len(mapping):,}")
    return mapping


def audit_reference(modality: str, holdout_by_sequence: Dict[str, Set[str]],
                    cell_line_map: Dict[str, Set[str]]) -> Tuple[Dict[str, object], List[Tuple[str, ...]], List[Tuple[str, ...]], List[Tuple[str, ...]]]:
    """Audit held-out pairs against union of NetMHCpan c000-c004 training files."""
    exact: Set[Tuple[str, str, str, str]] = set()
    possible: Set[Tuple[str, str, str, str]] = set()
    sequence_hits: Set[str] = set()
    files_scanned = 0
    rows_scanned = 0
    direct_hla_rows = 0
    cell_line_rows = 0
    unrecognised_rows = 0

    for config in NETMHCPAN_CONFIGS:
        path = REFERENCE_DIR / f"{config}_{modality}"
        files_scanned += 1
        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                fields = text.split()
                if len(fields) < 3:
                    continue
                peptide = normalise_peptide(fields[0])
                if not peptide or peptide not in holdout_by_sequence:
                    continue
                sequence_hits.add(peptide)
                raw_allele = fields[-1]
                allele = normalise_allele(raw_allele)
                if allele:
                    direct_hla_rows += 1
                    if allele in holdout_by_sequence[peptide]:
                        exact.add((path.name, peptide, allele, raw_allele))
                elif raw_allele in cell_line_map:
                    cell_line_rows += 1
                    matching = holdout_by_sequence[peptide] & cell_line_map[raw_allele]
                    for hla_allele in matching:
                        possible.add((path.name, peptide, hla_allele, raw_allele))
                else:
                    unrecognised_rows += 1
        print(f"  Audited {path.name}")

    exact_pairs = {(peptide, allele) for _, peptide, allele, _ in exact}
    possible_pairs = {(peptide, allele) for _, peptide, allele, _ in possible}
    exact_sequences = {peptide for peptide, _ in exact_pairs}
    possible_sequences = {peptide for peptide, _ in possible_pairs}
    no_allele_specific_match = sequence_hits - exact_sequences - possible_sequences
    summary = {
        "modality": "binding_BA" if modality == "ba" else "elution_EL",
        "reference_files_scanned": files_scanned,
        "reference_rows_with_heldout_sequence": direct_hla_rows + cell_line_rows + unrecognised_rows,
        "direct_HLA_reference_rows_with_heldout_sequence": direct_hla_rows,
        "cell_line_reference_rows_with_heldout_sequence": cell_line_rows,
        "unmapped_reference_rows_with_heldout_sequence": unrecognised_rows,
        "heldout_unique_peptides": len(holdout_by_sequence),
        "heldout_unique_peptide_allele_pairs": sum(len(v) for v in holdout_by_sequence.values()),
        "shared_peptide_sequences_any_reference_row": len(sequence_hits),
        "exact_peptide_allele_overlaps": len(exact_pairs),
        "possible_multiallelic_peptide_allele_overlaps": len(possible_pairs - exact_pairs),
        "shared_sequences_without_allele_specific_match": len(no_allele_specific_match),
        "exact_overlap_pairs_are_confirmed_allele_specific": True,
        "possible_multiallelic_matches_are_not_confirmed": True,
    }
    exact_rows = sorted(exact)
    possible_rows = sorted(possible)
    sequence_rows = [(modality, peptide, ";".join(sorted(holdout_by_sequence[peptide])),
                      "exact" if peptide in exact_sequences else "possible_multiallelic" if peptide in possible_sequences else "sequence_only")
                     for peptide in sorted(sequence_hits)]
    return summary, exact_rows, possible_rows, sequence_rows


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[Tuple[object, ...]]) -> None:
    with open_csv_write(path) as handle:
        writer = csv.writer(handle)
        writer.writerow(fieldnames)
        writer.writerows(rows)


def main() -> None:
    require_inputs()
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    print("========== SEQUENCE-GROUPED HOLDOUTS ==========")
    print(f"Split: {100-HOLDOUT_PERCENT}% train / {HOLDOUT_PERCENT}% holdout | seed={SEED}")
    print("All records for the same peptide sequence are assigned to the same partition.")

    binding_counts, _ = split_binding()
    elution_counts, _, elution_holdout_by_sequence = split_elution()

    # For binding, load held-out pairs after writing the split.
    binding_holdout_path = BINDING_DIR / "iedb_2chain_holdout.csv"
    binding_by_sequence: Dict[str, Set[str]] = defaultdict(set)
    with open_csv_text(binding_holdout_path) as handle:
        for row in csv.DictReader(handle):
            peptide = normalise_peptide(row.get("peptide"))
            allele = normalise_allele(row.get("hla_allele"))
            if peptide and allele:
                binding_by_sequence[peptide].add(allele)

    cell_line_map = load_cell_line_alleles()
    print("\n========== NETMHCPAN HELD-OUT OVERLAP AUDIT ==========")
    ba_summary, ba_exact, ba_possible, ba_sequences = audit_reference("ba", binding_by_sequence, cell_line_map)
    el_summary, el_exact, el_possible, el_sequences = audit_reference("el", elution_holdout_by_sequence, cell_line_map)

    # Output split report.
    split_rows = [("all", "split_method", "SHA-256 peptide-sequence grouping"),
                  ("all", "holdout_percent", HOLDOUT_PERCENT),
                  ("all", "seed", SEED),
                  ("all", "grouping_scope", "Same peptide sequence receives the same partition across binding and elution") ]
    for modality, counts in (("binding", binding_counts), ("elution", elution_counts)):
        for metric, value in counts.items():
            split_rows.append((modality, metric, value))
    write_csv(TABLES_DIR / "heldout_split_summary.csv", ["modality", "metric", "value"], split_rows)

    overlap_summaries = [ba_summary, el_summary]
    summary_fields = list(overlap_summaries[0].keys())
    write_csv(TABLES_DIR / "netmhcpan_heldout_overlap_summary.csv", summary_fields,
              [tuple(row.get(field, "") for field in summary_fields) for row in overlap_summaries])

    exact_rows = [("binding_BA", *row) for row in ba_exact] + [("elution_EL", *row) for row in el_exact]
    possible_rows = [("binding_BA", *row) for row in ba_possible] + [("elution_EL", *row) for row in el_possible]
    sequence_rows = ba_sequences + el_sequences
    write_csv(TABLES_DIR / "netmhcpan_heldout_exact_overlaps.csv",
              ["modality", "reference_file", "peptide", "iedb_heldout_allele", "reference_allele"], exact_rows)
    write_csv(TABLES_DIR / "netmhcpan_heldout_possible_multiallelic_overlaps.csv",
              ["modality", "reference_file", "peptide", "iedb_heldout_allele", "reference_cell_line"], possible_rows)
    write_csv(TABLES_DIR / "netmhcpan_heldout_sequence_overlaps.csv",
              ["reference_modality", "peptide", "iedb_heldout_alleles", "overlap_class"], sequence_rows)

    print("\n========== FINAL SUMMARY ==========")
    print(f"Binding holdout: {binding_counts['holdout_rows']:,} rows; {binding_counts['holdout_unique_peptides']:,} unique peptides")
    print(f"Elution holdout: {elution_counts['holdout_rows']:,} rows; {elution_counts['holdout_unique_peptides']:,} unique peptides")
    for summary in overlap_summaries:
        print(f"{summary['modality']}: exact={summary['exact_peptide_allele_overlaps']:,}; "
              f"possible multi-allelic={summary['possible_multiallelic_peptide_allele_overlaps']:,}; "
              f"sequence-only={summary['shared_sequences_without_allele_specific_match']:,}")
    print("\nOutputs:")
    for path in (
        BINDING_DIR / "iedb_2chain_train.csv",
        BINDING_DIR / "iedb_2chain_holdout.csv",
        ELUTION_DIR / "iedb_hla_i_elution_train.csv.gz",
        ELUTION_DIR / "iedb_hla_i_elution_holdout.csv.gz",
        TABLES_DIR / "heldout_split_summary.csv",
        TABLES_DIR / "netmhcpan_heldout_overlap_summary.csv",
        TABLES_DIR / "netmhcpan_heldout_exact_overlaps.csv",
        TABLES_DIR / "netmhcpan_heldout_possible_multiallelic_overlaps.csv",
        TABLES_DIR / "netmhcpan_heldout_sequence_overlaps.csv",
    ):
        print(f"  {path.relative_to(ROOT)}")
    print("\nPASS: Holdouts created and audited. Review exact and possible overlap categories separately.")


if __name__ == "__main__":
    main()
