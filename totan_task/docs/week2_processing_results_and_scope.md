# Week 2 — Processing benchmark: results, overlap and scope

**Workstream:** antigen processing (Totan)
**Date:** 2026-11-25
**Reference:** NetChop 3.1, network 20S-3.0 (`-v 1`), threshold 0.5
**Evaluation set:** CALIS 2015 Supplementary Table S1 — 16 HIV-1 GAG/TAT peptides

Companion documents: `week2_processing_reference_decision.md` (why this tool and
network), `week2_processing_input_target.md` (what is predicted and how it is
joined).

---

## 1. Headline result

NetChop 20S separates cleaved from non-cleaved bonds better than chance on both
proteasome conditions, but modestly.

| Proteasome | Sites | Cleaved | AUC | 95% CI | MCC | 95% CI |
|---|---|---|---|---|---|---|
| Immunoproteasome | 368 | 140 | **0.726** | 0.642 – 0.795 | 0.359 | 0.235 – 0.466 |
| Constitutive | 368 | 139 | **0.681** | 0.609 – 0.749 | 0.229 | 0.115 – 0.339 |

At the default threshold of 0.5:

| Proteasome | Sens | Spec | PPV | NPV | Acc | TP | FP | TN | FN |
|---|---|---|---|---|---|---|---|---|---|
| Immunoproteasome | 0.693 | 0.675 | 0.567 | 0.782 | 0.682 | 97 | 74 | 154 | 43 |
| Constitutive | 0.612 | 0.624 | 0.497 | 0.726 | 0.620 | 85 | 86 | 143 | 54 |

Full table, including the C-term comparison and the `X`-peptide sensitivity
check: `results/tables/netchop_cleavage_benchmark_by_proteasome_type.csv`.
Per-peptide breakdown: `..._by_peptide.csv`. Joined site-level data:
`netchop_cleavage_joined_sites.csv`.

### What this does and does not say

**Does:** NetChop 20S carries real signal on in vitro HIV peptide digestion —
both lower CI bounds sit above 0.5. On the constitutive condition its positive
predictions are right roughly half the time (PPV 0.497), so at the default
threshold a predicted cleavage site is close to a coin flip.

**Does not:** establish that NetChop performs better on immunoproteasome than
constitutive digestion. The intervals overlap substantially (0.642–0.795 against
0.609–0.749). The apparent 0.045 AUC gap is **not** statistically distinguishable
at this sample size and must not be reported as a finding.

**Does not:** establish a general cleavage-prediction performance level. Sixteen
peptides from one experiment on one organism's proteins is not a basis for a
general claim. See §3.

---

## 2. Overlap assessment

The binding workstream demonstrated "0 shared sequences" against NetMHCpan's
training partitions. **An equivalent audit is not possible here, and the claim
made is deliberately different in kind.**

NetChop's training data is described in the literature but is not distributed
with the tool or published anywhere.

| Network | Named training substrates | Verdict | Leakage risk |
|---|---|---|---|
| **20S-3.0** (used) | enolase, beta-casein | `DISJOINT_AT_SOURCE_PROTEIN_LEVEL` | low |
| Cterm-3.0 (not used) | ~1260 MHC-I ligands, not enumerated | `CANNOT_BE_ESTABLISHED` | unquantifiable |

**The claim supported for 20S:** no shared source protein. The training
substrates are enolase and beta-casein; the evaluation set is HIV-1 GAG and TAT.
No CALIS cleavage site can have been a 20S training example.

**The claim NOT supported:** site-level non-overlap. No site list exists publicly,
so this cannot be verified the way the binding arm verified it. The correct
wording is *"no shared source proteins; site-level verification not possible as
NetChop's training sets are not distributed"* — not *"0 shared sites"*.

**Why C-term would have been worse.** Its ligand corpus is described
inconsistently: Kesmir et al. 2002 say SYFPEITHI and MHCPEP, which contain HIV
ligands, while IEDB says "188 human proteins", which would not. HIV GAG epitopes
are among the most heavily characterised MHC-I ligands in existence. Contamination
can be neither confirmed nor excluded, which is precisely why 20S was chosen.

Machine-readable: `results/tables/netchop_calis_overlap_assessment.csv` and
`netchop_calis_source_protein_inventory.csv`.

---

## 3. Scope and limitations

Stated here rather than buried, because several would otherwise be invisible to a
reader of the results table alone.

### 3.1 Sixteen peptides

368 positions per condition sound like a reasonable sample. They are not
independent: they come from **16 peptides** in a single published experiment.
Confidence intervals throughout use a cluster bootstrap resampling whole peptides
for this reason. A naive per-position bootstrap would have produced intervals
roughly a third as wide and materially overstated precision.

This is a small-scale reference point, not a general performance estimate.

### 3.2 The CALIS source discrepancy is preserved, not resolved

Literal transcription of Table S1 gives I = 140, C = 139, 94 shared positions.
The article's Results section reports 148 / 150 / 103; its Methods section
reports 150 / 148 / 103 — the Results and Methods totals are themselves swapped.

No attempt has been made to reconcile these. The labels used here are the literal
table transcription. **If the article's totals are correct, roughly 6–8% of
cleavage sites are missing from the labels**, which would depress measured
sensitivity. This is unresolved and requires source-level review.

See `results/tables/calis_2015_cleavage_source_reconciliation.csv` and
`CALIS_CLEAVAGE_README.md`.

### 3.3 NetChop has no immunoproteasome model

The I/C stratification is of the **ground truth**, not the predictor. One
prediction set is evaluated against two label sets. NetChop cannot be said to
"predict immunoproteasome cleavage" — no such model exists in the distribution,
which contains exactly `Cterm.3.0` and `20S.3.0`.

### 3.4 This reproduces published work

The CALIS 2015 article (Immunogenetics, doi:10.1007/s00251-014-0815-0) already
benchmarks NetChop on this data and reports that in vitro digestion is best
captured by methods trained on in vitro data. This benchmark is an independent
reproduction establishing a local baseline with documented settings — not a novel
result. It should not be presented as one.

### 3.5 Non-standard residue retained

`nat10mod` contains a literal `X` as printed in Table S1. It is kept. All results
are reported both with and without that peptide; excluding it moves AUC by
roughly 0.01–0.02 and changes no conclusion.

### 3.6 Default threshold, not an optimised one

0.5 is NetChop's default and was not tuned on this data. Tuning it here would
convert the evaluation set into a fitting set. A tuned threshold would raise the
binarised metrics and would not be comparable to anything.

### 3.7 Upstream documentation is unreliable

Five contradictions were found in this tool's own documentation and packaging:

1. IEDB's `predict.py` implies `-v 0` selects 20S; the binary's man page and
   runtime banner say `-v 0` is C-term. **Following IEDB silently selects the
   wrong network.**
2. The default `-v 0` is C-term — the wrong network for in vitro data.
3. DTU describes 20S training data as *in vitro*; IEDB as *in vivo*.
4. DTU/Kesmir describe the C-term corpus as SYFPEITHI/MHCPEP; IEDB as 188 human
   proteins.
5. The shipped regression fixture is broken: `test/test.fsa` is Influenza A M2
   while `test/test.out` holds predictions for a different protein, dated 2006.

Mitigations: the flag is pinned explicitly, the banner-reported network is
recorded per run in `netchop_run_config.json`, and the container selftest asserts
the two networks produce different output so an ignored flag fails loudly.

---

## 4. Reproducing this

```bash
# 1. obtain IEDB_NetChop-3.0.tar.gz separately (licensed; not in this repo)
#    and place it in docker/netchop/
docker build -t sparq/netchop:3.1 docker/netchop
docker run --rm sparq/netchop:3.1 selftest        # expect SELFTEST PASS

# 2. pipeline
python ./scripts/prepare_calis_netchop_input.py   # FASTA + id map
python ./scripts/run_netchop_benchmark.py         # runs 20S and C-term
python ./scripts/evaluate_cleavage_benchmark.py   # stratified metrics
python ./scripts/assess_netchop_calis_overlap.py  # overlap tables
```

Run configuration, including container image id and the banner-reported network
for each run, is written to
`data/processed/cleavage/netchop_run_config.json`.

---

## 5. Open items

| Item | Owner | Why it matters |
|---|---|---|
| Confirm I/C means stratifying ground truth, not running two models | checkpoint | NetChop has no immunoproteasome model; the alternative reading is not executable |
| CALIS Table S1 vs article totals | source review | up to ~8% of positive labels may be missing |
| NetChop licensing status for a commercial entity | outside workstream | Non-Profit OSL is scoped by entity, not research intent |
| NetMHCpan training data committed to this **public** repo | binding workstream | same licence family; may be redistribution |
| Whether an integrated-pathway baseline is also wanted | project | NetCTLpan would answer a different question, against epitope labels |
