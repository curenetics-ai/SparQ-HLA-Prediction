# CALIS 2015 cleavage dataset preparation

## Source

User-provided `calis_2015_supp_table_s1.pdf`, page 1, Supporting Table S1: *In vitro proteasome digestion cleavage patterns* (Calis et al., 2015; DOI: 10.1007/s00251-014-0815-0).

The table caption describes 16 HIV-1 GAG/TAT peptides, 17-30 amino acids long, digested in vitro and analyzed by mass spectrometry after four hours. `I` means immunoproteasome and `C` means constitutive proteasome. Each printed vertical bar is treated as a cleavage boundary between adjacent residues.

## Files

- `scripts/prepare_calis_cleavage.py` - reproducible, standard-library-only script. Run from the project root with `python .\scripts\prepare_calis_cleavage.py`.
- `data/processed/cleavage/calis_2015_cleavage_patterns.csv` - 32 peptide/condition records, with sequence, fragment pattern and 1-based cleavage boundary positions.
- `data/processed/cleavage/calis_2015_cleavage_sites.csv` - one row per cleavage boundary explicitly marked in Table S1.
- `data/processed/cleavage/calis_2015_cleavage_site_labels.csv` - all 368 possible adjacent-residue positions for each of the two proteasome conditions (736 rows), labelled `cleaved` or `not_cleaved` according to the bars printed in Table S1.
- `data/processed/cleavage/calis_2015_cleavage_provenance.json` - source and processing metadata.
- `results/tables/calis_2015_cleavage_validation_summary.csv` - structural validation summary.
- `results/tables/calis_2015_cleavage_source_reconciliation.csv` - comparison of Table S1 bar counts with article aggregate counts.

## Position convention

`cleavage_after_residue_position_1based = 7` means a cut between residue 7 and residue 8. It is a 1-based residue boundary, not a zero-based Python string index.

## Important source discrepancy

Literal transcription of the bars in the supplied Supporting Table S1 yields:

- Immunoproteasome (`I`): 140 annotated cuts.
- Constitutive proteasome (`C`): 139 annotated cuts.
- Shared cut positions: 94.

The main article's Results section reports 148 (`I`), 150 (`C`) and 103 shared, while its Methods section reports 150 (`I`), 148 (`C`) and 103 shared. Thus, the supplementary table and main-text aggregates do not agree, and the Results and Methods totals are themselves swapped. The outputs preserve the supplied table faithfully; no additional cleavage sites have been invented to force a match. Review `calis_2015_cleavage_source_reconciliation.csv` before using these annotations as definitive training labels.

## Transcription notes

- The `X` in peptide `nat10mod` is retained as printed.
- Whitespace within the immunoproteasome `tat1` cell is treated as typesetting whitespace where no vertical bar is printed.
- The two proteasome-condition sequences are checked to reconstruct the same peptide sequence for each peptide name.
