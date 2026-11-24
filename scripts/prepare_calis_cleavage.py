#!/usr/bin/env python3
"""Prepare a source-faithful CALIS 2015 Supporting Table S1 dataset.

Run from the SparQ-HLA-Prediction project root:
    python .\\scripts\\prepare_calis_cleavage.py

The vertical bars in the supplied table denote inferred cleavage boundaries.
The script does not add cleavage sites to force agreement with article-level
summary statistics; those discrepancies are recorded in a reconciliation file.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

PROJECT_ROOT = Path.cwd()
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "cleavage"
RESULTS_DIR = PROJECT_ROOT / "results" / "tables"

# Fragments transcribed from the user-provided CALIS 2015 Supporting Table S1.
# I = immunoproteasome; C = constitutive proteasome.
# The fragment boundaries are the printed vertical bars in the source table.
SOURCE_PATTERNS = {
    "tat2": {
        "I": ["CFHCQVC", "FI", "TK", "GLGISY", "G", "RKK", "RR"],
        "C": ["CFHCQVC", "FITK", "GLGISY", "GRKK", "RR"],
    },
    "p seq3": {
        "I": ["RLIY", "A", "T", "R", "QL", "Q", "R", "F", "A", "V", "N", "PGL", "LI", "T"],
        "C": ["RLIYA", "T", "R", "QL", "Q", "R", "F", "A", "V", "N", "PGL", "LI", "T"],
    },
    "tat1": {
        # Spaces within these source cells are typesetting spaces, not cuts.
        "I": ["MEPVD", "PRLEPWKHPG", "SQPKTA", "CTNC", "Y", "C", "K"],
        "C": ["MEPVD", "PRLEPWKH", "PG", "SQPKTA", "C", "TN", "C", "Y", "C", "K"],
    },
    "nat12mod": {
        "I": ["RWL", "L", "L", "GL", "NPLV", "G", "GGR", "L", "Y", "SPTSI", "L", "G"],
        "C": ["RWLLLGLNPLV", "GGGRLYSPTSILG"],
    },
    "nat12": {
        "I": ["KRWIILGL", "NKIVRMYSPVSIL", "D"],
        "C": ["K", "R", "WIIL", "G", "L", "NK", "I", "V", "R", "M", "Y", "SPV", "S", "I", "L", "D"],
    },
    "p seq2": {
        "I": ["YVL", "F", "L", "T", "K", "GL", "SI", "SY", "L", "GKK"],
        "C": ["Y", "V", "L", "F", "L", "T", "K", "GL", "S", "I", "SYL", "GKK"],
    },
    "p24": {
        "I": ["ALSEGATPQD", "LNTML", "NTV", "GGHQA", "AMQML"],
        "C": ["ALSE", "GATP", "Q", "DLNTM", "L", "NTVGGHQA", "AMQML"],
    },
    "p seq1": {
        "I": ["FVIH", "R", "L", "EPWL", "HPG", "SQHI", "TA", "S", "TN"],
        "C": ["FVIH", "R", "L", "E", "PWL", "H", "PGSQ", "H", "I", "T", "A", "STN"],
    },
    "nat10mod": {
        "I": ["R", "FII", "PXF", "T", "A", "L", "SGGRR", "A", "L", "L", "Y", "GATPY", "AI", "G"],
        "C": ["R", "FII", "PXFT", "A", "L", "SGGRR", "A", "L", "L", "Y", "GA", "TPYA", "I", "G"],
    },
    "p seq4": {
        "I": ["YAIP", "Q", "A", "L", "N", "T", "L", "L", "N", "TV", "GGHQAA"],
        "C": ["YAIP", "Q", "A", "L", "N", "T", "L", "L", "N", "TV", "GGHQAA"],
    },
    "p17": {
        "I": ["YKLK", "HI", "VW", "A", "S", "RELER", "F", "AVNPGL", "L", "E", "V", "TS", "E", "GC"],
        "C": ["YKLKHIVWASREL", "ERFAVNPGLLEVTSEGC"],
    },
    "nat13mod": {
        "I": ["R", "AL", "GPA", "A", "TL", "QTPWTA", "SL", "GV", "G"],
        "C": ["R", "AL", "G", "PA", "A", "T", "L", "Q", "T", "P", "WTA", "SLG", "VG"],
    },
    "nat13": {
        "I": ["KALGPA", "A", "TL", "EEMM", "T", "A", "CQGVGGPGH"],
        "C": ["KALGPAATL", "EEMM", "T", "A", "CQGVGGPGH"],
    },
    "nat11mod": {
        "I": ["RAIPIPA", "GTL", "L", "SGGGR", "AIYK", "R", "W", "AI", "L", "G"],
        "C": ["R", "A", "I", "PIP", "A", "GT", "LL", "SGGGR", "A", "IY", "K", "R", "W", "A", "I", "L", "G"],
    },
    "nat11": {
        "I": ["NNPP", "IPVG", "E", "I", "Y", "K", "R", "W", "II", "L", "G", "L", "N", "KI", "V"],
        "C": ["NNPPIPVG", "EIY", "K", "R", "W", "IIL", "G", "L", "NKI", "V"],
    },
    "nat10": {
        "I": ["PEVIPMF", "S", "AL", "SE", "GATPQ", "D", "L", "NTML", "NTVGGH"],
        "C": ["PEVIPMFS", "ALSE", "GATPQ", "D", "LNTML", "NTVGGH"],
    },
}

TYPE_LABELS = {
    "I": "immunoproteasome",
    "C": "constitutive_proteasome",
}
VALID_RESIDUES = set("ACDEFGHIKLMNPQRSTVWYX")

# Article-level counts are recorded as published, not treated as labels to impose.
PUBLISHED_SUMMARIES = {
    "results_section": {"I": 148, "C": 150, "shared": 103},
    "methods_section": {"I": 150, "C": 148, "shared": 103},
}


def get_cut_positions(fragments: list[str]) -> list[int]:
    """1-based residue position immediately before each | cleavage boundary."""
    positions: list[int] = []
    running = 0
    for fragment in fragments[:-1]:
        running += len(fragment)
        positions.append(running)
    return positions


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(dict.fromkeys(key for row in rows for key in row.keys()))
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    patterns: list[dict] = []
    cleaved_sites: list[dict] = []
    all_site_labels: list[dict] = []
    errors: list[str] = []
    peptide_sequence_by_name: dict[str, str] = {}
    cut_positions_by_condition: dict[tuple[str, str], set[int]] = {}

    for peptide_name, by_type in SOURCE_PATTERNS.items():
        if set(by_type) != {"I", "C"}:
            errors.append(f"{peptide_name}: missing I or C pattern")
            continue

        seq_i = "".join(by_type["I"])
        seq_c = "".join(by_type["C"])
        if seq_i != seq_c:
            errors.append(f"{peptide_name}: I/C sequence mismatch: {seq_i} != {seq_c}")
            continue

        sequence = seq_i
        peptide_sequence_by_name[peptide_name] = sequence
        invalid = sorted(set(sequence) - VALID_RESIDUES)
        if invalid:
            errors.append(f"{peptide_name}: invalid residue symbols {invalid}")
        if not 17 <= len(sequence) <= 30:
            errors.append(f"{peptide_name}: sequence length {len(sequence)} is outside 17-30")

        for code in ("I", "C"):
            fragments = by_type[code]
            cuts = get_cut_positions(fragments)
            cut_positions_by_condition[(peptide_name, code)] = set(cuts)
            patterns.append({
                "source": "CALIS_2015_Supplementary_Table_S1",
                "peptide_name": peptide_name,
                "proteasome_code": code,
                "proteasome_type": TYPE_LABELS[code],
                "sequence": sequence,
                "sequence_length": len(sequence),
                "cleavage_pattern": "|".join(fragments),
                "fragments_json": json.dumps(fragments),
                "cleavage_after_residue_positions_1based": ";".join(map(str, cuts)),
                "cleavage_site_count": len(cuts),
                "digestion_time_hours": 4,
                "assay_method": "mass spectrometry",
                "source_note": "Each printed vertical bar in Supporting Table S1 marks an inferred cleavage site.",
            })

            for site_index, cut_after in enumerate(cuts, start=1):
                cleaved_sites.append({
                    "source": "CALIS_2015_Supplementary_Table_S1",
                    "peptide_name": peptide_name,
                    "proteasome_code": code,
                    "proteasome_type": TYPE_LABELS[code],
                    "sequence": sequence,
                    "cleavage_site_index": site_index,
                    "cleavage_after_residue_position_1based": cut_after,
                    "left_residue": sequence[cut_after - 1],
                    "right_residue": sequence[cut_after],
                    "cleavage_context": f"{sequence[cut_after - 1]}|{sequence[cut_after]}",
                    "is_cleaved": 1,
                })

            # Every adjacent residue pair is a possible site; the printed bars
            # determine the positive/negative label in the supplied table.
            for cut_after in range(1, len(sequence)):
                is_cleaved = int(cut_after in set(cuts))
                all_site_labels.append({
                    "source": "CALIS_2015_Supplementary_Table_S1",
                    "peptide_name": peptide_name,
                    "proteasome_code": code,
                    "proteasome_type": TYPE_LABELS[code],
                    "sequence": sequence,
                    "sequence_length": len(sequence),
                    "candidate_cleavage_after_residue_position_1based": cut_after,
                    "left_residue": sequence[cut_after - 1],
                    "right_residue": sequence[cut_after],
                    "cleavage_context": f"{sequence[cut_after - 1]}|{sequence[cut_after]}",
                    "is_cleaved": is_cleaved,
                    "site_label": "cleaved" if is_cleaved else "not_cleaved",
                })

    # Aggregate only the boundaries explicitly present in Table S1.
    table_counts = {
        code: sum(len(cut_positions_by_condition.get((name, code), set())) for name in SOURCE_PATTERNS)
        for code in ("I", "C")
    }
    shared_count = sum(
        len(cut_positions_by_condition[(name, "I")] & cut_positions_by_condition[(name, "C")])
        for name in SOURCE_PATTERNS
    )
    possible_sites_per_condition = sum(len(seq) - 1 for seq in peptide_sequence_by_name.values())

    validation_rows = [
        {"check": "source_peptides", "observed": len(SOURCE_PATTERNS), "expected": 16, "status": "PASS" if len(SOURCE_PATTERNS) == 16 else "FAIL"},
        {"check": "peptide_condition_pattern_rows", "observed": len(patterns), "expected": 32, "status": "PASS" if len(patterns) == 32 else "FAIL"},
        {"check": "all_sequences_length_17_to_30", "observed": sum(17 <= len(s) <= 30 for s in peptide_sequence_by_name.values()), "expected": 16, "status": "PASS" if len(peptide_sequence_by_name) == 16 and all(17 <= len(s) <= 30 for s in peptide_sequence_by_name.values()) else "FAIL"},
        {"check": "I_C_sequence_pairs_match", "observed": len(peptide_sequence_by_name), "expected": 16, "status": "PASS" if len(peptide_sequence_by_name) == 16 and not errors else "FAIL"},
        {"check": "possible_sites_per_proteasome_condition", "observed": possible_sites_per_condition, "expected": 368, "status": "PASS" if possible_sites_per_condition == 368 else "FAIL"},
        {"check": "all_candidate_site_labels", "observed": len(all_site_labels), "expected": 736, "status": "PASS" if len(all_site_labels) == 736 else "FAIL"},
        {"check": "source_transcription_checks", "observed": len(errors), "expected": 0, "status": "PASS" if not errors else "FAIL"},
        {"check": "Table_S1_immunoproteasome_cut_count", "observed": table_counts["I"], "expected": "See source reconciliation; table-derived count", "status": "REVIEW"},
        {"check": "Table_S1_constitutive_cut_count", "observed": table_counts["C"], "expected": "See source reconciliation; table-derived count", "status": "REVIEW"},
        {"check": "Table_S1_shared_cut_count", "observed": shared_count, "expected": "See source reconciliation; table-derived count", "status": "REVIEW"},
    ]

    if errors:
        (RESULTS_DIR / "calis_2015_cleavage_validation_errors.txt").write_text("\n".join(errors) + "\n", encoding="utf-8")
        write_csv(RESULTS_DIR / "calis_2015_cleavage_validation_summary.csv", validation_rows)
        raise SystemExit("Source transcription validation failed. See results/tables/calis_2015_cleavage_validation_errors.txt")

    reconciliation_rows = []
    for section_name, section_counts in PUBLISHED_SUMMARIES.items():
        for code, label in TYPE_LABELS.items():
            reported = section_counts[code]
            observed = table_counts[code]
            reconciliation_rows.append({
                "comparison": "Table_S1_vs_article_summary",
                "article_section": section_name,
                "proteasome_code": code,
                "proteasome_type": label,
                "table_S1_bar_count": observed,
                "article_reported_count": reported,
                "difference_table_minus_article": observed - reported,
                "status": "MATCH" if observed == reported else "DISCREPANCY_REQUIRES_SOURCE_REVIEW",
            })
        reconciliation_rows.append({
            "comparison": "Table_S1_vs_article_summary",
            "article_section": section_name,
            "proteasome_code": "BOTH",
            "proteasome_type": "shared_cleavage_positions",
            "table_S1_bar_count": shared_count,
            "article_reported_count": section_counts["shared"],
            "difference_table_minus_article": shared_count - section_counts["shared"],
            "status": "MATCH" if shared_count == section_counts["shared"] else "DISCREPANCY_REQUIRES_SOURCE_REVIEW",
        })

    write_csv(OUTPUT_DIR / "calis_2015_cleavage_patterns.csv", patterns)
    write_csv(OUTPUT_DIR / "calis_2015_cleavage_sites.csv", cleaved_sites)
    write_csv(OUTPUT_DIR / "calis_2015_cleavage_site_labels.csv", all_site_labels)
    write_csv(RESULTS_DIR / "calis_2015_cleavage_validation_summary.csv", validation_rows)
    write_csv(RESULTS_DIR / "calis_2015_cleavage_source_reconciliation.csv", reconciliation_rows)

    provenance = {
        "source": "Calis et al. (2015), user-provided Supporting Table S1 PDF, page 1",
        "doi": "10.1007/s00251-014-0815-0",
        "source_table_title": "In vitro proteasome digestion cleavage patterns",
        "source_peptide_count": 16,
        "proteasome_conditions": {"I": "immunoproteasome", "C": "constitutive proteasome"},
        "pattern_records": len(patterns),
        "explicit_cleavage_annotations_from_bars": len(cleaved_sites),
        "candidate_site_label_records": len(all_site_labels),
        "table_S1_cleavage_counts": {"I": table_counts["I"], "C": table_counts["C"], "shared": shared_count},
        "possible_sites_per_condition": possible_sites_per_condition,
        "sequence_length_range": "17-30 amino acids, as stated in the table caption",
        "digestion_time_hours": 4,
        "cleavage_position_convention": "1-based residue position immediately before the cleavage boundary",
        "important_reconciliation_note": "The bar counts transcribed from Table S1 do not match the main article's aggregate counts; the article Results and Methods sections also swap the I/C aggregate totals. See calis_2015_cleavage_source_reconciliation.csv. No sites were invented to force agreement.",
        "transcription_notes": [
            "Each printed vertical bar is preserved as a cleavage boundary.",
            "Spaces inside the tat1 immunoproteasome cell were treated as typesetting spaces where no vertical bar is printed.",
            "The X symbol in nat10mod is preserved as shown in the source.",
        ],
    }
    (OUTPUT_DIR / "calis_2015_cleavage_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    print("========== CALIS CLEAVAGE PREPARATION ==========")
    print(f"Source peptides:                     {len(SOURCE_PATTERNS):,}")
    print(f"Peptide-condition patterns:          {len(patterns):,}")
    print(f"Possible sites per condition:        {possible_sites_per_condition:,}")
    print(f"All candidate-site labels:           {len(all_site_labels):,}")
    print(f"Table S1 cleavage annotations:       I={table_counts['I']:,}; C={table_counts['C']:,}")
    print(f"Table S1 shared cleavage positions:  {shared_count:,}")
    print(f"Source transcription errors:         {len(errors):,}")
    print(f"Patterns: {OUTPUT_DIR / 'calis_2015_cleavage_patterns.csv'}")
    print(f"Sites:    {OUTPUT_DIR / 'calis_2015_cleavage_sites.csv'}")
    print(f"Labels:   {OUTPUT_DIR / 'calis_2015_cleavage_site_labels.csv'}")
    print(f"Summary:  {RESULTS_DIR / 'calis_2015_cleavage_validation_summary.csv'}")
    print(f"Audit:    {RESULTS_DIR / 'calis_2015_cleavage_source_reconciliation.csv'}")
    print("PASS: Structural checks passed. Review the source reconciliation before model use.")


if __name__ == "__main__":
    main()
