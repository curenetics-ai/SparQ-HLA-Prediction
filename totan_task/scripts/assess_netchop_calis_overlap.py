#!/usr/bin/env python3
"""
Assess overlap between the CALIS 2015 evaluation set and NetChop's training data.

Week 2, antigen-processing reference benchmark.

The binding arm could demonstrate "0 shared sequences" by exact comparison,
because NetMHCpan ships its training partitions. NetChop does not: its training
sets are described in the literature but are not distributed. The claim this
benchmark can support is therefore different in kind, and the table this script
writes records that explicitly rather than implying an equivalent audit.

What CAN be established: NetChop 20S-3.0 was trained on degradation data for
enolase and beta-casein. CALIS 2015 is HIV-1 GAG and TAT. Those are disjoint
source proteins, so no CALIS cleavage site can have been a 20S training example.

What CANNOT be established: site-level non-overlap, because no site list exists
publicly. And for Cterm-3.0 the position is weaker still - its ligand corpus is
described inconsistently (SYFPEITHI/MHCPEP, which contain HIV ligands, versus
"188 human proteins", which would not), so HIV GAG contamination can be neither
confirmed nor excluded. That asymmetry is itself a reason to benchmark on 20S.

Outputs:
    results/tables/netchop_calis_overlap_assessment.csv
    results/tables/netchop_calis_source_protein_inventory.csv

Run from the project root:
    python ./scripts/assess_netchop_calis_overlap.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

TASK_ROOT = Path(__file__).resolve().parents[1]   # totan_task/
REPO_ROOT = Path(__file__).resolve().parents[2]   # repository root (Week 1 inputs)
ROOT = TASK_ROOT                                  # this task's outputs
ID_MAP = ROOT / "data" / "interim" / "calis_peptide_id_map.csv"
TABLES = ROOT / "results" / "tables"

# Source proteins of the CALIS 2015 Table S1 peptides, per the article: HIV-1
# GAG (p17/p24 and derived natural/modified variants) and HIV-1 TAT.
CALIS_SOURCE_ORGANISM = "Human immunodeficiency virus 1"
CALIS_SOURCE_PROTEINS = ["HIV-1 GAG (p17, p24)", "HIV-1 TAT"]

# NetChop training data as described by the distributors. Deliberately recorded
# with the conflicting descriptions intact.
NETCHOP_TRAINING = [
    {
        "network": "20S-3.0",
        "training_description": "In vitro 20S proteasome degradation data",
        "named_source_proteins": "enolase; beta-casein",
        "source_organisms": "non-HIV (yeast/bovine model substrates)",
        "distributed_with_tool": "no",
        "documentation_source": "DTU NetChop page; IEDB tool documentation",
        "documentation_conflict": (
            "DTU describes the 20S training data as in vitro; IEDB describes it "
            "as in vivo degradation data for enolase and beta-casein"
        ),
    },
    {
        "network": "Cterm-3.0",
        "training_description": "C-terminal cleavage sites of ~1260 MHC class I ligands",
        "named_source_proteins": "not enumerated publicly",
        "source_organisms": "disputed - see conflict",
        "distributed_with_tool": "no",
        "documentation_source": "DTU NetChop page; Kesmir et al. 2002; IEDB tool documentation",
        "documentation_conflict": (
            "Kesmir 2002 states ligands were drawn from SYFPEITHI and MHCPEP, "
            "which include viral (incl. HIV) ligands; IEDB describes the set as "
            "ligands from 188 human proteins, which would exclude HIV"
        ),
    },
]


def main() -> int:
    if not ID_MAP.exists():
        sys.exit(f"missing {ID_MAP}\nRun scripts/prepare_calis_netchop_input.py first.")
    TABLES.mkdir(parents=True, exist_ok=True)

    with ID_MAP.open(encoding="utf-8-sig", newline="") as fh:
        peptides = list(csv.DictReader(fh))

    # Inventory of what is actually being evaluated.
    inv = TABLES / "netchop_calis_source_protein_inventory.csv"
    with inv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "peptide_name", "netchop_id", "sequence_length", "n_internal_bonds",
            "source_organism", "source_proteins", "has_nonstandard_residue"])
        w.writeheader()
        for p in peptides:
            w.writerow({
                "peptide_name": p["peptide_name"],
                "netchop_id": p["netchop_id"],
                "sequence_length": p["sequence_length"],
                "n_internal_bonds": p["n_internal_bonds"],
                "source_organism": CALIS_SOURCE_ORGANISM,
                "source_proteins": "; ".join(CALIS_SOURCE_PROTEINS),
                "has_nonstandard_residue": p["has_nonstandard_residue"],
            })

    rows = []
    for t in NETCHOP_TRAINING:
        if t["network"] == "20S-3.0":
            verdict = "DISJOINT_AT_SOURCE_PROTEIN_LEVEL"
            basis = (
                "Named training substrates (enolase, beta-casein) share no source "
                "protein with HIV-1 GAG/TAT, so no CALIS site can be a training "
                "example. Site-level verification not possible."
            )
            risk = "low"
        else:
            verdict = "CANNOT_BE_ESTABLISHED"
            basis = (
                "Training ligand list is not published and its source corpora are "
                "described inconsistently. If drawn from SYFPEITHI/MHCPEP, HIV GAG "
                "ligands may be present, which would place CALIS C-terminal sites "
                "inside the training data. Neither confirmable nor excludable."
            )
            risk = "unquantifiable"
        rows.append({
            "evaluation_set": "CALIS_2015_Supplementary_Table_S1",
            "evaluation_source_organism": CALIS_SOURCE_ORGANISM,
            "evaluation_source_proteins": "; ".join(CALIS_SOURCE_PROTEINS),
            "n_peptides_evaluated": len(peptides),
            "n_sites_per_proteasome_condition": sum(
                int(p["n_internal_bonds"]) for p in peptides),
            **t,
            "exact_site_level_check_possible": "no",
            "overlap_verdict": verdict,
            "basis": basis,
            "leakage_risk": risk,
            "comparable_to_netmhcpan_audit": (
                "no - the NetMHCpan audit compared against distributed training "
                "partitions; no equivalent artefact exists for NetChop"
            ),
        })

    out = TABLES / "netchop_calis_overlap_assessment.csv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {inv.relative_to(ROOT)}  ({len(peptides)} peptides)")
    print(f"wrote {out.relative_to(ROOT)}")
    for r in rows:
        print(f"  {r['network']:10} {r['overlap_verdict']:34} risk={r['leakage_risk']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
