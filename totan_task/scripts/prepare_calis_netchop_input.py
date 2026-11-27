#!/usr/bin/env python3
"""
Build the NetChop input FASTA from the CALIS 2015 cleavage labels.

Week 2, antigen-processing reference benchmark.

Why this script exists rather than feeding the CSV straight in:

1. Peptide names in Table S1 contain spaces (`p seq1`). NetChop echoes only the
   first whitespace-delimited token of a FASTA header in its `Ident` column, so
   `p seq1` and `p seq2` would both come back as `p` and become unjoinable.
   Every peptide therefore gets a safe surrogate id (CALIS01..CALIS16) and an
   explicit id map is written alongside.
2. `nat10mod` contains a literal `X` as printed in the source table. It is kept
   (the project does not silently alter source data) but flagged here and in the
   id map so downstream steps can account for it.
3. The two proteasome conditions share one sequence per peptide, so NetChop is
   run once per peptide, not once per condition. NetChop has no
   immunoproteasome/constitutive model - the I/C split belongs to the ground
   truth, not to the predictor.

Outputs:
    data/interim/calis_peptides.fasta
    data/interim/calis_peptide_id_map.csv

Run from the project root:
    python ./scripts/prepare_calis_netchop_input.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

TASK_ROOT = Path(__file__).resolve().parents[1]   # totan_task/
REPO_ROOT = Path(__file__).resolve().parents[2]   # repository root (Week 1 inputs)
ROOT = TASK_ROOT                                  # this task's outputs
LABELS = REPO_ROOT / "data" / "processed" / "cleavage" / "calis_2015_cleavage_site_labels.csv"
OUT_DIR = ROOT / "data" / "interim"
FASTA = OUT_DIR / "calis_peptides.fasta"
ID_MAP = OUT_DIR / "calis_peptide_id_map.csv"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def load_peptides() -> dict[str, str]:
    """peptide_name -> sequence, asserting one sequence per peptide."""
    if not LABELS.exists():
        sys.exit(f"missing input: {LABELS}\nRun scripts/prepare_calis_cleavage.py first.")

    seqs: dict[str, set[str]] = {}
    with LABELS.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            seqs.setdefault(row["peptide_name"], set()).add(row["sequence"])

    out: dict[str, str] = {}
    for name, variants in seqs.items():
        if len(variants) != 1:
            sys.exit(f"peptide {name!r} has {len(variants)} different sequences: {variants}")
        out[name] = variants.pop()
    return out


def main() -> int:
    peptides = load_peptides()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    with FASTA.open("w", encoding="utf-8", newline="\n") as fh:
        for i, (name, seq) in enumerate(sorted(peptides.items()), start=1):
            safe_id = f"CALIS{i:02d}"
            nonstandard = sorted(set(seq) - STANDARD_AA)
            # The header holds ONLY the surrogate id. NetChop truncates the
            # identifier to a fixed width, so any description bleeds into its
            # `Ident` column (">CALIS01 nat10" comes back as "CALIS01 na").
            # The human-readable name lives in the id map instead.
            fh.write(f">{safe_id}\n{seq}\n")
            rows.append({
                "netchop_id": safe_id,
                "peptide_name": name,
                "sequence": seq,
                "sequence_length": len(seq),
                "n_internal_bonds": len(seq) - 1,
                "nonstandard_residues": "".join(nonstandard),
                "has_nonstandard_residue": int(bool(nonstandard)),
            })

    with ID_MAP.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    total_bonds = sum(r["n_internal_bonds"] for r in rows)
    flagged = [r for r in rows if r["has_nonstandard_residue"]]

    print(f"wrote {FASTA.relative_to(ROOT)}  ({len(rows)} peptides)")
    print(f"wrote {ID_MAP.relative_to(ROOT)}")
    print(f"internal bonds per proteasome condition: {total_bonds}")
    if flagged:
        print("non-standard residues retained as printed in Table S1:")
        for r in flagged:
            print(f"  {r['peptide_name']}: {r['nonstandard_residues']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
