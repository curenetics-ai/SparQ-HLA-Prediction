# Week 2 Antigen-Processing Benchmark Report

SparQ HLA Prediction · processing workstream · 27 November 2026

## Summary

The antigen-processing reference benchmark is complete. NetChop 3.1 (network 20S-3.0) was benchmarked against the CALIS 2015 in vitro cleavage dataset, stratified by proteasome type, in a reproducible containerised run.

| Proteasome | Sites | Cleaved | AUC | 95% CI | MCC |
| --- | --- | --- | --- | --- | --- |
| Immunoproteasome | 368 | 140 | 0.726 | 0.642 – 0.795 | 0.359 |
| Constitutive | 368 | 139 | 0.681 | 0.609 – 0.749 | 0.229 |

NetChop carries real but modest signal: both lower confidence bounds sit above 0.5, while at the default threshold a predicted constitutive cleavage site is correct roughly half the time (PPV 0.497). The apparent immunoproteasome advantage is **not** statistically distinguishable at this sample size and is not reported as a finding.

Two decisions were required where the plan named one. Beyond choosing NetChop over NetCTLpan, NetChop ships two networks and the default is the wrong one for this data. Selecting it correctly turned out to depend on a flag that the vendor's own wrapper documents backwards.

Five separate documentation and packaging defects were found in the NetChop distribution. One of them silently selects the wrong network for anyone who follows the vendor's instructions.

All four deliverables are written and the pipeline reproduces byte-identically from a clean run. Four items remain open that require decisions outside this workstream.

## The assignment

Week 2 ran two parallel workstreams, by design neither waiting on the other. The binding arm benchmarked NetMHCpan-4.1; this report covers the processing arm.

Three tasks were set:

1. **Confirm the comparator independently.** NetChop was preferred for cleavage-only separation; document the broader NetCTLpan trade-off.
2. **Define a cleavage-site holdout and check reference-training overlap.** With an explicit warning: do not use NetChop to score eluted-ligand data, because it predicts cleavage.
3. **Run cleavage scoring by proteasome type (I vs C).** Preserve the CALIS discrepancy and report the 16-peptide scope.

The project README specified four deliverables: a processing-reference decision and rationale, an input/target definition, an overlap and compatibility note, and stratified benchmark results.

The Week 2 target applied to both arms: benchmark each reference on its own held-out data, record coverage, sample counts and method settings, and do not report only a pooled score.

## Decision 1 — NetChop rather than NetCTLpan

NetCTLpan was listed as the alternative. Checking it independently showed it is not an independent one.

**NetCTLpan contains NetChop.** Its cleavage component is the NetChop C-term network, combined with TAP transport efficiency and MHC class I binding as a weighted average. Benchmarking NetCTLpan would mean benchmarking NetChop plus two further components, one of which is HLA binding — the other workstream's subject.

That makes the separation argument structural rather than a matter of preference:

- A NetCTLpan score cannot attribute performance to processing, because binding signal is mixed into it by construction.
- The two workstreams would cease to be independent. Any later comparison of a Curenetics processing method against this baseline would be partly a comparison against a binding predictor.

**The trade-off accepted.** NetCTLpan models more of the presentation pathway, so it is the better tool for predicting epitopes. It is the wrong tool for measuring cleavage prediction, which is this arm's responsibility. If an integrated-pathway baseline is wanted later, NetCTLpan is the right choice for that question and should be benchmarked separately against epitope-level labels, not cleavage-site labels.

## Decision 2 — the 20S network rather than C-term

This choice was not in the plan. NetChop 3.1 ships two networks, its `etc/` directory containing exactly `Cterm.3.0` and `20S.3.0`, and the default is the wrong one for this data.

### Grounds

**Data-type compatibility.** 20S was trained on proteasome degradation data for enolase and beta-casein; CALIS Table S1 records which internal bonds were cut during in vitro digestion. That is 20S's task. C-term is trained on the C-terminal ends of MHC class I ligands, which answers a different question. The CALIS article itself reports that in vitro digestion is best captured by methods trained on in vitro data, while in vivo proteolysis is best captured by C-term.

**Proteasome-type labelling.** IEDB describes the 20S training data as deriving from the constitutive proteasome, corresponding directly to the CALIS C condition. No equivalent correspondence exists for C-term, which has no stated proteasome-type provenance.

**Leakage exposure.** This is the ground most relevant to the Week 1 standard. 20S's named training substrates, enolase and beta-casein, share no source protein with HIV-1 GAG or TAT, so no CALIS site can have been a training example. C-term's ligand corpus is not published and is described inconsistently — Kesmir et al. 2002 say it was drawn from SYFPEITHI and MHCPEP, both of which contain HIV ligands, while IEDB describes it as ligands from 188 human proteins, which would exclude HIV. HIV GAG epitopes are among the best-characterised MHC-I ligands in existence. Contamination can therefore be neither confirmed nor excluded for C-term.

**Empirical confirmation.** Both networks were run so the decision could be shown rather than asserted. 20S scored higher on both conditions (AUC 0.681 and 0.726, against 0.636 and 0.720 for C-term), as the data-type argument predicts. Intervals overlap, so this supports the decision rather than proving it.

### The flag is documented backwards upstream

This is the most dangerous detail in the workstream. The authoritative mapping, from the binary's own man page and confirmed in its runtime banner, is:

```
-v value:  value = 0 use Cterm-3.0.  value = 1 use 20S-3.0.  Default 0.
```

IEDB's bundled `predict.py` contradicts this. Its validation error reads *"pred_method must be method id: '0', '1', corresponding method name: '20S', 'C term'"*, which parses as 0 meaning 20S. The wrapper passes that value straight through as `-v`, so `pred_method='0'` runs C-term.

Anyone following the vendor's own message therefore selects `0` believing they chose 20S, and silently gets C-term — which is both the wrong network for in vitro data and the one carrying unquantifiable HIV leakage risk. None of this is visible in the output unless the banner line is read.

The runner pins the flag explicitly, records the banner-reported network for every run, and the container selftest asserts that the two networks produce different scores, so an ignored flag fails loudly instead of producing plausible wrong numbers.

## Method

**Prediction task.** Given a peptide sequence, predict for each internal peptide bond whether the proteasome cleaves there. Binary classification over positions, not over peptides.

**Input.** The 16 CALIS peptides as standalone sequences, 17–30 aa, in one FASTA. Standalone rather than embedded in parent GAG/TAT protein, because the experiment digested isolated synthetic peptides; adding flanking context would give NetChop information the experiment did not have. Peptides use surrogate identifiers CALIS01–CALIS16, because NetChop truncates FASTA identifiers to a fixed width and several Table S1 names contain spaces — `>CALIS01 nat10` returns as `CALIS01 na`, which would collide in the join.

**Target.** The 736-row label set: 368 internal bonds × 2 proteasome conditions. Immunoproteasome 140 cleaved / 228 not; constitutive 139 / 229. Prevalence is close to 38% in both.

**Position alignment — the critical join.** NetChop scores every residue 1..L, each meaning "cleavage after this residue". CALIS annotates internal bonds 1..L-1 only. The difference is position L, the peptide's own C-terminus, which is an artefact of synthesis rather than an observed proteasomal cut. NetChop position i maps to CALIS bond i for i in 1..L-1, and **position L is dropped** — 16 rows per network, leaving exactly 368, matching the label count.

This is the step most likely to produce silently wrong results, so it is verified rather than assumed. The evaluator asserts for every joined row that NetChop's residue at position i equals CALIS's left residue at bond i; a one-position offset would mismatch amino acids and abort the run. It also asserts the join yields exactly the expected row count, so neither duplication nor silent loss can pass. Both assertions pass.

**Proteasome type is a property of the labels, not the predictor.** NetChop has no immunoproteasome model. One set of 20S predictions is evaluated twice, against the I ground truth and the C ground truth. The predictions are identical in both strata; only the label column differs.

**Holdout.** The entire CALIS set is the holdout and is not split. No model is trained, fitted, calibrated or threshold-tuned on it. This is stated explicitly because the binding workstream's artefacts carry train/holdout splits and a reader moving between the two could assume the same structure. The NetChop default threshold of 0.5 is used throughout; tuning it on this data would convert the holdout into a fitting set.

**Uncertainty.** Confidence intervals use a cluster bootstrap resampling whole peptides, 2000 replicates. The 368 positions come from only 16 peptides and positions within a peptide are not independent — they share a sequence, an experiment and overlapping NetChop input windows. Resampling positions would treat correlated observations as independent and produce intervals roughly a third as wide.

**Reproducibility.** The shipped netChop is an ELF 32-bit i386 statically linked Linux executable and cannot run natively on Windows, so the run is containerised against a pinned Debian image. The full pipeline reproduces byte-identically from a clean rerun.

## Results

Benchmark network 20S-3.0, threshold 0.5, all 16 peptides.

| Proteasome | Sites | Cleaved | AUC | AUC 95% CI | MCC | MCC 95% CI |
| --- | --- | --- | --- | --- | --- | --- |
| Immunoproteasome | 368 | 140 | 0.726 | 0.642 – 0.795 | 0.359 | 0.235 – 0.466 |
| Constitutive | 368 | 139 | 0.681 | 0.609 – 0.749 | 0.229 | 0.115 – 0.339 |

At the default threshold:

| Proteasome | Sens | Spec | PPV | NPV | Accuracy | TP | FP | TN | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Immunoproteasome | 0.693 | 0.675 | 0.567 | 0.782 | 0.682 | 97 | 74 | 154 | 43 |
| Constitutive | 0.612 | 0.624 | 0.497 | 0.726 | 0.620 | 85 | 86 | 143 | 54 |

The C-term network, run for comparison only, scored AUC 0.720 (immunoproteasome) and 0.636 (constitutive) — lower on both.

### What the numbers support

NetChop 20S carries real signal on in vitro HIV peptide digestion. Both lower confidence bounds sit above 0.5, so the separation is not attributable to chance at this sample size.

The signal is modest in practical terms. On the constitutive condition, PPV is 0.497: at the default threshold a predicted cleavage site is right about as often as a coin flip. Specificity is also low in both strata, around 0.62–0.68, meaning a substantial fraction of uncut bonds are predicted as cut.

### What the numbers do not support

**They do not establish that NetChop performs better on immunoproteasome than constitutive digestion.** The confidence intervals overlap substantially (0.642–0.795 against 0.609–0.749). The apparent 0.045 AUC gap is not statistically distinguishable at 16 peptides and must not be reported as a finding.

**They do not establish a general cleavage-prediction performance level.** Sixteen peptides from a single experiment on one organism's proteins is a local reference point, not a general estimate.

### Sensitivity check

One peptide, nat10mod, contains a literal X as printed in Table S1. It was retained rather than corrected or dropped. Excluding it moves AUC by 0.01–0.02 (20S constitutive 0.681 to 0.663; immunoproteasome 0.726 to 0.716) and changes no conclusion. Both variants are reported in the results table rather than one replacing the other.

## Overlap assessment

The binding workstream demonstrated zero shared sequences against NetMHCpan's training partitions, because NetMHCpan ships them. **NetChop does not.** Its training sets are described in the literature but are not distributed with the tool or published anywhere.

An equivalent audit is therefore impossible, and the claim made here is deliberately different in kind rather than a weaker version of the same claim.

| Network | Named training substrates | Verdict | Leakage risk |
| --- | --- | --- | --- |
| 20S-3.0 (used) | enolase, beta-casein | Disjoint at source-protein level | Low |
| Cterm-3.0 (not used) | ~1260 MHC-I ligands, not enumerated | Cannot be established | Unquantifiable |

**What is supported for 20S.** No shared source protein. The training substrates are enolase and beta-casein; the evaluation set is HIV-1 GAG and TAT. No CALIS cleavage site can have been a 20S training example. This argument is structural rather than empirical — it does not depend on sequences happening not to collide, but on the corpora being incapable of intersecting.

**What is not supported.** Site-level non-overlap. No site list exists publicly. The correct wording is "no shared source proteins; site-level verification not possible as NetChop's training sets are not distributed" — not "zero shared sites".

That limitation is itself a finding worth recording: the reference tool's training data being unavailable is a material difference between the two workstreams' audits, and anyone comparing the two arms' claims should know the stronger-sounding one rests on a distributed artefact that simply does not exist on this side.

## Defects found in the NetChop distribution

Five contradictions were found in the tool's own documentation and packaging. They are recorded because each one is capable of producing results that look correct and are not, and because anyone else in the group using NetChop will hit them.

| # | Defect | Consequence |
| --- | --- | --- |
| 1 | IEDB's predict.py implies -v 0 selects 20S; the binary's man page and runtime banner say -v 0 is C-term | Following the vendor silently selects the wrong network |
| 2 | The default is -v 0, i.e. C-term | The default is wrong for in vitro digestion data |
| 3 | DTU describes the 20S training data as in vitro; IEDB describes it as in vivo | The data-compatibility argument rests on a disputed description |
| 4 | DTU and Kesmir 2002 describe the C-term corpus as SYFPEITHI/MHCPEP; IEDB describes it as 188 human proteins | Determines whether HIV contamination is possible; unresolvable |
| 5 | The shipped regression fixture is broken: test.fsa is Influenza A M2, while test.out holds predictions for a different protein, dated 2006 | The vendor's own correctness check cannot be used |

Defect 1 is the serious one. It combines with defect 2 so that both the vendor's instructions and the vendor's default lead to C-term — the network that is wrong on data-type grounds and carries the HIV leakage exposure. Nothing in the output signals this except the banner line.

**Mitigations applied.** The -v flag is pinned explicitly in the runner. The banner-reported network is captured for every run in the run-configuration file, so the network actually used is recorded rather than assumed. Because the vendor fixture is unusable, the container selftest verifies seven properties instead: that predictions are produced, that -v 0 and -v 1 give different scores, that each banner names the network requested, that output is deterministic across runs, that scores lie within [0,1], and that the cleaved/not-cleaved flag agrees with the threshold. A silently ignored flag now fails loudly.

## Scope and limitations

Stated here rather than buried, because several would otherwise be invisible to a reader of the results table alone.

**Sixteen peptides.** 368 positions per condition sounds like a reasonable sample, but they are not independent: they come from 16 peptides in a single published experiment. This is why the confidence intervals use a cluster bootstrap. A naive per-position bootstrap would have produced intervals about a third as wide and would have materially overstated precision. This is a small-scale reference point, not a general performance estimate.

**The CALIS source discrepancy is preserved, not resolved.** Literal transcription of Table S1 gives I = 140, C = 139, 94 shared positions. The article's Results section reports 148 / 150 / 103; its Methods section reports 150 / 148 / 103 — the Results and Methods totals are themselves swapped. No attempt has been made to reconcile them, in line with the project's policy of not silently adjusting source data. If the article's totals are correct, roughly 6–8% of cleavage sites are missing from the labels used here, which would depress measured sensitivity. This remains unresolved and needs source-level review.

**NetChop has no immunoproteasome model.** The distribution contains exactly Cterm.3.0 and 20S.3.0. The I/C stratification is of the ground truth, not the predictor. NetChop cannot be described as predicting immunoproteasome cleavage. Published work exists that maps C-term to immunoproteasome and 20S to constitutive; the second half has some basis in IEDB's description, the first half has none, and the mapping was not followed here.

**This reproduces published work.** The CALIS 2015 article already benchmarks NetChop on this data and reports that in vitro digestion is best captured by methods trained on in vitro data. This exercise is an independent reproduction establishing a local baseline with documented settings and a verified pipeline. It is not a novel result and should not be presented as one.

**Default threshold, not an optimised one.** 0.5 is NetChop's default and was not tuned. Tuning it on this data would convert the evaluation set into a fitting set, and a tuned figure would not be comparable to anything.

**Eluted-ligand data was not used.** The Week 2 brief warned against scoring EL ligands with NetChop, and that warning is correct: NetChop scores positions in a sequence for whether a cut occurs there, whereas eluted ligands are complete peptides. Despite 1.49M elution records being available in the repository against only 736 cleavage labels, they answer a different question and were not used.

## Questions raised and open items

None of these blocks the Friday deliverables. All require a decision or review outside this workstream.

| # | Question | Owner | Why it matters |
| --- | --- | --- | --- |
| 1 | Does "run cleavage scoring by proteasome type (I vs C)" mean stratifying the ground truth, rather than running two models? | Checkpoint / Suraj | NetChop has no immunoproteasome model, so the alternative reading is not executable. The work assumed the former. |
| 2 | Should the CALIS Table S1 versus article discrepancy be resolved at source? | Source review | If the article totals are right, ~6–8% of positive labels are missing, depressing measured sensitivity. |
| 3 | Is Curenetics covered by the Non-Profit OSL 3.0 terms for NetChop? | Suraj / contracts | The licence is scoped by entity, not research intent. The vendor directs non-academic users to contact DTU. |
| 4 | Is committing NetMHCpan-4.1 training data to a public repository a redistribution problem? | Binding workstream | Same licence family. The reference directory is already public. |
| 5 | Is an integrated-pathway baseline also wanted? | Project | NetCTLpan would answer a different question and needs epitope-level labels, not cleavage-site labels. |

### Questions worth putting to the binding workstream

These arose from differences between the two arms rather than from the processing work itself:

- **How was NetMHCpan-4.1 licensed, and under what status?** Same DTU licence wall; the answer settles question 3 as well.
- **What metric and table schema is the binding arm reporting?** The two arms should stack rather than look like separate projects. This report uses AUC with MCC, sensitivity, specificity and PPV at a fixed threshold, with support counts on every row.
- **Is the binding arm reporting confidence intervals?** It matters more here: 736 rows from 16 peptides will look falsely precise as point estimates.
- **Note that the two arms' overlap claims are not equivalent.** The binding arm can say "zero shared sequences"; this arm can only say "no shared source proteins". That asymmetry should be understood before the two are presented side by side.

## Deliverables and reproduction

### Against the Friday output

| Required | Delivered |
| --- | --- |
| Reproducible processing run | Pinned Docker image, 7-check selftest, run config captured per run, byte-identical on clean rerun |
| I/C-stratified results | netchop_cleavage_benchmark_by_proteasome_type.csv, plus per-peptide and site-level tables |
| Overlap / scope note | netchop_calis_overlap_assessment.csv and the scope section of the results document |

### Files

Documents:

- docs/week2_processing_reference_decision.md — tool and network decision, rationale, licensing
- docs/week2_processing_input_target.md — input, target, position alignment, holdout
- docs/week2_processing_results_and_scope.md — results, overlap, limitations

Scripts:

- scripts/prepare_calis_netchop_input.py — FASTA and identifier map
- scripts/run_netchop_benchmark.py — runs both networks, parses output, records configuration
- scripts/evaluate_cleavage_benchmark.py — join assertions, stratified metrics, cluster bootstrap
- scripts/assess_netchop_calis_overlap.py — overlap tables

Results tables:

- netchop_cleavage_benchmark_by_proteasome_type.csv
- netchop_cleavage_benchmark_by_peptide.csv
- netchop_cleavage_joined_sites.csv
- netchop_calis_overlap_assessment.csv
- netchop_calis_source_protein_inventory.csv

Container:

- docker/netchop/Dockerfile and entrypoint.sh

### Reproducing the run

```
docker build -t sparq/netchop:3.1 docker/netchop
docker run --rm sparq/netchop:3.1 selftest

python ./scripts/prepare_calis_netchop_input.py
python ./scripts/run_netchop_benchmark.py
python ./scripts/evaluate_cleavage_benchmark.py
python ./scripts/assess_netchop_calis_overlap.py
```

The NetChop distribution must be obtained separately and placed in docker/netchop/. It is licensed software and is **not** committed to the repository, which is public; .gitignore excludes the archive and any extracted binaries. Run configuration, including the container image id and the banner-reported network for each run, is written to data/processed/cleavage/netchop_run_config.json.
