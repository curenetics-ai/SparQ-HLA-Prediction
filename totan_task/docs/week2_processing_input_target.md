# Week 2 — Processing benchmark: input, target and holdout definition

**Workstream:** antigen processing (Totan)
**Date:** 2026-11-25

This document fixes what is being predicted, from what, and what counts as held
out. It exists because three of these choices are silent failure modes: get them
wrong and the pipeline still runs and still produces plausible numbers.

---

## 1. Prediction task

> Given a peptide sequence, predict for each **internal peptide bond** whether
> the proteasome cleaves there.

Binary classification over positions, not over peptides.

---

## 2. Input

**What NetChop receives:** each CALIS peptide as a standalone sequence, 16
sequences, 17–30 aa, in one FASTA
(`data/interim/calis_peptides.fasta`).

**Why standalone and not the parent protein.** The CALIS experiment digested
isolated synthetic peptides in vitro. Scoring the peptide alone is the faithful
analogue. Embedding each peptide back into HIV-1 GAG or TAT would give NetChop
flanking context the experiment did not have, and would change the prediction at
every position near the termini.

**Why surrogate identifiers.** NetChop truncates the FASTA identifier to a fixed
width in its `Ident` column. Table S1 peptide names include spaces (`p seq1`,
`p seq2`), and a header of `>CALIS01 nat10` returns as `CALIS01 na`. Peptides are
therefore written with id-only headers `CALIS01`–`CALIS16`, with the mapping in
`data/interim/calis_peptide_id_map.csv`. Without this, several peptides would
collide in the join or be silently dropped.

**Non-standard residue.** `nat10mod` contains a literal `X` at position 6, as
printed in Table S1. It is retained rather than corrected or dropped, consistent
with the project's treatment of source data. NetChop accepts it and scores the
sequence. All benchmark results are additionally reported with the
`X`-containing peptide excluded, as a sensitivity check.

**One run per peptide, not per condition.** The I and C rows of Table S1 share
one sequence per peptide. NetChop has no proteasome-type parameter, so running it
twice would produce identical output. It is run once per network over 16
sequences.

---

## 3. Target

**Source:** `data/processed/cleavage/calis_2015_cleavage_site_labels.csv`
(736 rows = 368 internal bonds × 2 proteasome conditions).

| Field | Meaning |
|---|---|
| `candidate_cleavage_after_residue_position_1based` | 1-based bond: `7` = between residues 7 and 8 |
| `is_cleaved` | 1 if Table S1 prints a bar at that bond, else 0 |
| `proteasome_code` | `I` immunoproteasome, `C` constitutive |

| Condition | Sites | Cleaved | Not cleaved | Prevalence |
|---|---|---|---|---|
| Immunoproteasome | 368 | 140 | 228 | 38.0% |
| Constitutive | 368 | 139 | 229 | 37.8% |

---

## 4. Position alignment — the critical join

**This is the step most likely to produce silently wrong results.**

NetChop scores **every residue 1..L**, each meaning "cleavage after this
residue". CALIS annotates **internal bonds 1..L−1 only**.

The difference is position `L`: the peptide's own C-terminus. That is an artefact
of how the peptide was synthesised, not an observed proteasomal cut. It has no
CALIS label and must not be scored.

**Rule applied:** NetChop position *i* maps to CALIS bond *i* for
*i* ∈ [1, L−1]; **position L is dropped.**

Across 16 peptides this drops 16 rows per network (32 in total), leaving exactly
368 positions per network — matching the CALIS label count exactly.

### How the alignment is verified rather than assumed

`scripts/evaluate_cleavage_benchmark.py` asserts, for every joined row, that
NetChop's residue at position *i* equals CALIS's `left_residue` at bond *i*. A
one-position offset would put a different amino acid on at least some rows and
abort the run. The evaluation also asserts that the join produces exactly
`n_labels × n_networks` rows, so neither duplication nor silent loss can pass.

Both assertions pass on the current data.

---

## 5. Proteasome type is a property of the labels, not the predictor

NetChop 3.1 provides exactly two networks, `Cterm-3.0` and `20S-3.0`. **Neither
is immunoproteasome-specific, and no immunoproteasome model exists.**

Therefore "run cleavage scoring by proteasome type (I vs C)" is implemented as:

> one set of NetChop 20S predictions, evaluated **twice** — once against the
> immunoproteasome ground truth, once against the constitutive ground truth.

The predictions are identical in both strata. What differs is the label column.
Any difference in measured performance reflects how well a single
proteasome-type-agnostic model matches each digestion pattern.

**What would be wrong.** Running `-v 1` and `-v 0` and labelling them
"immunoproteasome" and "constitutive" results. IEDB does describe 20S's training
data as constitutive, so a 20S→C association has some basis; there is **no**
basis whatsoever for Cterm→I, since C-term is trained on ligand C-termini with no
proteasome-type provenance. Published work exists that makes this mapping; it
should not be followed here. Doing so would manufacture an "immunoproteasome
model" that does not exist.

---

## 6. Holdout definition

**The entire CALIS set is the holdout. It is not split.**

No model is trained, fitted, calibrated or selected on this data. This is an
external reference benchmark of a frozen third-party predictor, so a train/test
split would serve no purpose — there is nothing to train.

Stated explicitly because the binding workstream's artefacts carry train/holdout
splits, and a reader moving between the two could reasonably assume the same
structure here. All 736 labelled sites are evaluation data.

**Threshold.** The NetChop default of 0.5 is used for all binarised metrics and
is recorded in `netchop_run_config.json`. It was not tuned on this data; tuning a
threshold on the evaluation set would convert the holdout into a fitting set.

**Nothing is excluded for performance reasons.** The only subsetting reported is
the `X`-peptide sensitivity check, which is declared in advance and reported
alongside the full set rather than replacing it.

---

## 7. Metrics

Primary: **AUC** (threshold-free, so unaffected by the 0.5 choice).
Secondary at threshold 0.5: MCC, sensitivity, specificity, PPV, NPV, F1,
accuracy, and the full confusion matrix.

Every row carries `n_sites`, `n_cleaved` and `n_peptides`, so no figure appears
without its support count.

**Confidence intervals use a cluster bootstrap over peptides, not positions.**
The 368 positions come from only 16 peptides and positions within a peptide are
not independent — they share a sequence, a digestion experiment, and overlapping
NetChop input windows. Resampling positions would treat 368 correlated
observations as 368 independent ones and produce intervals that are too narrow.
Resampling whole peptides is wider and honest. 2000 replicates, seed recorded in
the script.

---

## Summary

| | Definition |
|---|---|
| Input | 16 standalone CALIS peptides, id-only FASTA headers |
| Unit of prediction | internal peptide bond |
| Target | `is_cleaved` from Table S1 bars |
| Positions | NetChop *i* → CALIS bond *i*, **position L dropped** |
| Verification | residue-identity assertion + exact row-count assertion |
| I vs C | stratification of the **labels**; one prediction set |
| Holdout | the entire 736-row set; no split, no tuning |
| Uncertainty | cluster bootstrap over 16 peptides |
