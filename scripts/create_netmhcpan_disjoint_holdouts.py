#!/usr/bin/env python3
"""Create conservative NetMHCpan-sequence-disjoint subsets from IEDB holdouts.

Run from the SparQ-HLA-Prediction project root:
    python scripts/create_netmhcpan_disjoint_holdouts.py

This script DOES NOT alter the original IEDB datasets or the original holdouts.
It creates additional benchmark subsets by excluding every peptide sequence
found in the corresponding NetMHCpan 4.1 reference modality (BA for binding,
EL for elution), including exact, possible multi-allelic, and sequence-only
hits. This is a conservative sequence-level exclusion policy.
"""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path.cwd()
BINDING_INPUT = ROOT / "data/processed/binding/iedb_2chain_holdout.csv"
ELUTION_INPUT = ROOT / "data/processed/elution/iedb_hla_i_elution_holdout.csv.gz"
OVERLAP_INPUT = ROOT / "results/tables/netmhcpan_heldout_sequence_overlaps.csv"

BINDING_OUTPUT = ROOT / "data/processed/binding/iedb_2chain_holdout_netmhcpan_disjoint.csv"
ELUTION_OUTPUT = ROOT / "data/processed/elution/iedb_hla_i_elution_holdout_netmhcpan_disjoint.csv.gz"
SUMMARY_OUTPUT = ROOT / "results/tables/netmhcpan_reference_disjoint_holdout_summary.csv"


def normalise_peptide(value: str | None) -> str:
    return (value or "").strip().upper().replace(" ", "")


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def open_output(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=6)
    return path.open("w", encoding="utf-8-sig", newline="")


def load_exclusions() -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    if not OVERLAP_INPUT.exists():
        raise FileNotFoundError(f"Missing overlap report: {OVERLAP_INPUT}")

    by_modality: dict[str, set[str]] = {
        "binding_BA": set(),
        "elution_EL": set(),
    }
    classes: dict[str, dict[str, str]] = {
        "binding_BA": {},
        "elution_EL": {},
    }
    # The upstream audit writes these values as "ba" and "el" in the
    # reference_modality column (not "binding_BA" / "elution_EL").
    modality_aliases = {
        "ba": "binding_BA",
        "binding_ba": "binding_BA",
        "binding_BA": "binding_BA",
        "el": "elution_EL",
        "elution_el": "elution_EL",
        "elution_EL": "elution_EL",
    }
    with OVERLAP_INPUT.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"reference_modality", "peptide", "overlap_class"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Unexpected columns in {OVERLAP_INPUT}: {reader.fieldnames}"
            )
        for row in reader:
            raw_modality = (row.get("reference_modality") or "").strip()
            modality = modality_aliases.get(raw_modality, modality_aliases.get(raw_modality.lower()))
            peptide = normalise_peptide(row.get("peptide"))
            category = (row.get("overlap_class") or "").strip()
            if modality in by_modality and peptide:
                by_modality[modality].add(peptide)
                classes[modality][peptide] = category

    empty_modalities = [name for name, sequences in by_modality.items() if not sequences]
    if empty_modalities:
        raise ValueError(
            "No overlap sequences were loaded for: " + ", ".join(empty_modalities)
            + f". Check the modality values in {OVERLAP_INPUT}."
        )
    return by_modality, classes


def filter_dataset(
    input_path: Path,
    output_path: Path,
    modality: str,
    excluded_sequences: set[str],
    category_by_sequence: dict[str, str],
) -> dict[str, int]:
    if not input_path.exists():
        raise FileNotFoundError(f"Missing holdout dataset: {input_path}")

    counts: Counter[str] = Counter()
    # Initialise keys explicitly so a modality with zero excluded rows still
    # produces a complete, printable summary instead of a KeyError.
    counts["excluded_rows"] = 0
    counts["retained_rows"] = 0
    counts["rows_without_peptide"] = 0
    original_sequences: set[str] = set()
    retained_sequences: set[str] = set()
    excluded_sequence_rows: Counter[str] = Counter()
    excluded_category_rows: Counter[str] = Counter()

    with open_text(input_path) as source, open_output(output_path) as destination:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError(f"No CSV header found in {input_path}")
        writer = csv.DictWriter(destination, fieldnames=reader.fieldnames)
        writer.writeheader()

        sequence_column = "peptide" if modality == "binding_BA" else "linear_sequence"
        if sequence_column not in reader.fieldnames:
            raise ValueError(
                f"Missing {sequence_column!r} in {input_path}: {reader.fieldnames}"
            )

        for row in reader:
            counts["original_rows"] += 1
            peptide = normalise_peptide(row.get(sequence_column))
            if not peptide:
                counts["rows_without_peptide"] += 1
                continue
            original_sequences.add(peptide)
            if peptide in excluded_sequences:
                counts["excluded_rows"] += 1
                category = category_by_sequence.get(peptide, "unclassified")
                excluded_category_rows[category] += 1
                excluded_sequence_rows[peptide] += 1
                continue

            writer.writerow(row)
            counts["retained_rows"] += 1
            retained_sequences.add(peptide)

    residual = retained_sequences & excluded_sequences
    if residual:
        raise RuntimeError(
            f"{modality}: {len(residual)} reference-overlapping sequences remain"
        )

    counts["original_unique_peptides"] = len(original_sequences)
    counts["excluded_unique_peptides"] = len(excluded_sequence_rows)
    counts["retained_unique_peptides"] = len(retained_sequences)
    counts["excluded_exact_overlap_rows"] = excluded_category_rows["exact"]
    counts["excluded_possible_multiallelic_rows"] = excluded_category_rows["possible_multiallelic"]
    counts["excluded_sequence_only_rows"] = excluded_category_rows["sequence_only"]
    counts["residual_shared_sequences"] = len(residual)
    return dict(counts)


def main() -> None:
    modality_exclusions, category_maps = load_exclusions()

    print("========== CREATE REFERENCE-DISJOINT HOLDOUTS ==========")
    print(f"Binding BA overlap sequences loaded: {len(modality_exclusions['binding_BA']):,}")
    print(f"Elution EL overlap sequences loaded:  {len(modality_exclusions['elution_EL']):,}")
    print("Policy: exclude every held-out peptide found in the corresponding")
    print("NetMHCpan training files, regardless of allele-match certainty.")

    binding = filter_dataset(
        BINDING_INPUT, BINDING_OUTPUT, "binding_BA",
        modality_exclusions["binding_BA"], category_maps["binding_BA"]
    )
    elution = filter_dataset(
        ELUTION_INPUT, ELUTION_OUTPUT, "elution_EL",
        modality_exclusions["elution_EL"], category_maps["elution_EL"]
    )

    rows = []
    for modality, metrics in (("binding_BA", binding), ("elution_EL", elution)):
        for metric, value in metrics.items():
            rows.append((modality, metric, value))
        rows.append((modality, "excluded_reference_sequences_in_audit", len(modality_exclusions[modality])))
        rows.append((modality, "exclusion_policy", "sequence-level; exact, possible multi-allelic, and sequence-only matches removed"))

    SUMMARY_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY_OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["modality", "metric", "value"])
        writer.writerows(rows)

    print("\n========== REFERENCE-DISJOINT SUMMARY ==========")
    for modality, metrics in (("binding_BA", binding), ("elution_EL", elution)):
        print(f"{modality}:")
        print(f"  Original rows:              {metrics['original_rows']:,}")
        print(f"  Excluded rows:              {metrics['excluded_rows']:,}")
        print(f"  Retained rows:              {metrics['retained_rows']:,}")
        print(f"  Original unique peptides:   {metrics['original_unique_peptides']:,}")
        print(f"  Excluded unique peptides:   {metrics['excluded_unique_peptides']:,}")
        print(f"  Retained unique peptides:   {metrics['retained_unique_peptides']:,}")
        print(f"  Residual shared sequences:  {metrics['residual_shared_sequences']:,}")

    print("\nOutputs:")
    print(f"  {BINDING_OUTPUT.relative_to(ROOT)}")
    print(f"  {ELUTION_OUTPUT.relative_to(ROOT)}")
    print(f"  {SUMMARY_OUTPUT.relative_to(ROOT)}")

    if binding["retained_rows"] == 0 or elution["retained_rows"] == 0:
        print("\nWARNING: One disjoint holdout is empty; review the exclusion policy.")
    else:
        print("\nPASS: Both reference-disjoint holdouts created; residual sequence overlaps = 0.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
