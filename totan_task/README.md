# Totan — Week 2: antigen-processing reference benchmark

Self-contained work for the processing workstream. Everything produced for this
task lives in this folder; **nothing outside it has been modified.**

The one input read from the main repository is Week 1's CALIS label file,
`data/processed/cleavage/calis_2015_cleavage_site_labels.csv`, which is read
only. Scripts resolve it via `REPO_ROOT`; all outputs are written under
`TASK_ROOT` (this folder).

---

## Headline result

NetChop 3.1, network **20S-3.0** (`-v 1`), threshold 0.5, against CALIS 2015
in vitro cleavage data — 16 HIV-1 GAG/TAT peptides.

| Proteasome | Sites | Cleaved | AUC | 95% CI | MCC |
|---|---|---|---|---|---|
| Immunoproteasome | 368 | 140 | 0.726 | 0.642 – 0.795 | 0.359 |
| Constitutive | 368 | 139 | 0.681 | 0.609 – 0.749 | 0.229 |

Real but modest signal: both lower bounds clear 0.5, while PPV on the
constitutive condition is 0.497. **The I-vs-C difference is not statistically
distinguishable** at this sample size and is not a finding.

---

## Layout

```
totan_task/
  reports/     the written report - PDF, Word and Markdown
  docs/        decision, input/target definition, results and scope
  scripts/     the pipeline (4 steps) plus the report builder
  results/     benchmark and overlap tables
  data/        this task's intermediate and prediction outputs
  docker/      pinned container for NetChop (binary NOT committed)
```

## Read in this order

1. `reports/Week2_Antigen_Processing_Benchmark_Report.pdf` — the whole thing, 11 pages
2. `docs/week2_processing_reference_decision.md` — why NetChop, and why the 20S network
3. `docs/week2_processing_input_target.md` — what is predicted, and the position-L join
4. `docs/week2_processing_results_and_scope.md` — results, overlap, limitations

## Reproducing

```bash
# obtain IEDB_NetChop-3.0.tar.gz separately (licensed) -> docker/netchop/
docker build -t sparq/netchop:3.1 docker/netchop
docker run --rm sparq/netchop:3.1 selftest      # expect SELFTEST PASS

python ./scripts/prepare_calis_netchop_input.py
python ./scripts/run_netchop_benchmark.py
python ./scripts/evaluate_cleavage_benchmark.py
python ./scripts/assess_netchop_calis_overlap.py
```

Run from this folder. Verified to reproduce byte-identically from a clean run.

## Three things to know before reading the numbers

**The `-v` flag is documented backwards by the vendor.** The binary's man page
and runtime banner say `-v 0` = Cterm-3.0 and `-v 1` = 20S-3.0, but IEDB's
bundled `predict.py` implies `0` means 20S. Following IEDB silently selects
C-term — the wrong network for this data, and the one with HIV leakage
exposure. The flag is pinned explicitly and the selftest asserts the two
networks differ.

**NetChop has no immunoproteasome model.** Only `Cterm.3.0` and `20S.3.0` exist.
The I/C split is a stratification of the ground truth; one prediction set is
evaluated against two label sets.

**The overlap claim is weaker than the binding arm's, by necessity.** NetChop's
training data is not distributed, so site-level verification is impossible. What
is shown is source-protein disjointness: 20S was trained on enolase and
beta-casein, the evaluation set is HIV-1 GAG/TAT.

## Open items

| # | Question | Owner |
|---|---|---|
| 1 | Confirm I/C means stratifying ground truth, not running two models | checkpoint |
| 2 | CALIS Table S1 vs article totals (~6–8% of positives may be missing) | source review |
| 3 | NetChop licensing status for a commercial entity | contracts |
| 4 | NetMHCpan training data committed to this public repo — redistribution? | binding workstream |
| 5 | Is an integrated-pathway (NetCTLpan) baseline also wanted? | project |
