# Week 2 — Antigen-processing reference: decision and rationale

**Workstream:** antigen processing (Totan)
**Decision date:** 2026-11-25
**Status:** decided; both options were run so the choice is demonstrated, not asserted

---

## Decision

**Reference tool: NetChop 3.1.**
**Network: 20S-3.0, selected with `netChop -v 1`.**

Two separate choices are involved. The project plan named only the first. The
second is not optional and determines whether the benchmark is valid.

---

## 1. NetChop rather than NetCTLpan

NetCTLpan was listed as the alternative. It is not an independent one.

NetCTLpan's cleavage component **is** the NetChop C-term network. Its output is a
weighted average of three scores: MHC binding, TAP transport efficiency, and
NetChop C-terminal cleavage. Benchmarking NetCTLpan would therefore mean
benchmarking NetChop plus two other components, one of which is HLA binding —
the binding workstream's subject.

That makes the separation argument structural rather than a matter of taste:

- A NetCTLpan score cannot attribute performance to processing, because binding
  signal is mixed into it by construction.
- The two workstreams would no longer be independent. Any later comparison of a
  Curenetics processing method against this baseline would be partly a
  comparison against a binding predictor.

**The trade-off being accepted.** NetCTLpan models more of the presentation
pathway, so it is the better tool for predicting *epitopes*. It is the wrong tool
for measuring *cleavage prediction*, which is what this arm is responsible for.
If the project later needs an integrated-pathway baseline, NetCTLpan is the right
choice for that question and should be benchmarked separately, against
epitope-level labels rather than cleavage-site labels.

---

## 2. The 20S network rather than C-term

NetChop 3.1 ships two networks. The distribution's `etc/` directory contains
exactly `20S.3.0` and `Cterm.3.0`.

Three independent arguments select 20S, and a fourth confirms it empirically.

### 2.1 Data-type compatibility

| | Trained on | CALIS 2015 is |
|---|---|---|
| **20S-3.0** | proteasome degradation data for enolase and beta-casein | in vitro proteasome digestion |
| **Cterm-3.0** | C-terminal cleavage sites of ~1260 MHC class I ligands | not ligand C-termini |

CALIS Table S1 records which internal bonds were cut when peptides were digested
in vitro and analysed by mass spectrometry. That is 20S's task. C-term answers a
different question — where naturally presented ligands end in vivo — and the
CALIS 2015 article itself reports that in vitro digestion patterns are best
captured by methods trained on in vitro data while in vivo proteolysis is best
captured by C-term.

### 2.2 Proteasome-type labelling

IEDB describes the 20S training data as deriving from the **constitutive**
proteasome. That corresponds directly to the CALIS `C` condition. No equivalent
correspondence exists for C-term, which is trained on ligand C-termini and has no
stated proteasome-type provenance.

Note what this does **not** license: there is no immunoproteasome network. Any
mapping of C-term to immunoproteasome would be unfounded. See
`week2_processing_input_target.md`.

### 2.3 Leakage exposure

This is the argument most relevant to the project's Week 1 standard.

- **20S** was trained on enolase and beta-casein. CALIS is HIV-1 GAG and TAT.
  Disjoint source proteins, so no CALIS site can have been a training example.
- **C-term** was trained on a ligand set that is not published. Kesmir et al.
  2002 state the ligands were drawn from SYFPEITHI and MHCPEP, both of which
  contain HIV ligands; HIV GAG epitopes are among the best-characterised MHC-I
  ligands in existence. IEDB instead describes the set as ligands from "188 human
  proteins", which would exclude HIV. These cannot both be correct.

So benchmarking on C-term would mean reporting a number with unquantifiable
contamination risk against an HIV evaluation set. Benchmarking on 20S permits a
verifiable disjointness claim. Full assessment in
`results/tables/netchop_calis_overlap_assessment.csv`.

### 2.4 Empirical confirmation

Both networks were run so the decision could be checked rather than assumed
(`results/tables/netchop_cleavage_benchmark_by_proteasome_type.csv`):

| Network | Proteasome | AUC | 95% CI |
|---|---|---|---|
| **20S** | constitutive | **0.681** | 0.609 – 0.749 |
| **20S** | immunoproteasome | **0.726** | 0.642 – 0.795 |
| Cterm | constitutive | 0.636 | 0.580 – 0.696 |
| Cterm | immunoproteasome | 0.720 | 0.661 – 0.778 |

20S scores higher on both conditions, as the data-type argument predicts. The
intervals overlap, so this is supporting evidence, not proof — the decision rests
on 2.1–2.3, and 2.4 shows nothing contradicts it.

---

## 3. The `-v` flag is documented backwards upstream

**This is the single most dangerous detail in this workstream.**

The authoritative mapping, from the binary's own man page (`netchop.1`) and
confirmed in its runtime banner:

```
-v value:  value = 0 use Cterm-3.0.  value = 1 use 20S-3.0.  Default 0.
```

IEDB's bundled `predict.py` contradicts this. Its validation error reads:

```
pred_method must be method id: '0', '1', corresponding method name: '20S', 'C term'.
```

which parses as `0` → `20S`. Line 33 of `netchop_3_1_executable/__init__.py`
passes that value straight through as `-v`, so `pred_method='0'` runs **C-term**.

Consequences for anyone following the wrapper's message:

1. They select `0` believing they chose 20S and silently get C-term.
2. C-term is the wrong network for in vitro digestion data.
3. C-term is the network carrying HIV leakage exposure.

All three failures are invisible in the output unless the banner line is read.
The runner (`scripts/run_netchop_benchmark.py`) therefore pins the flag
explicitly, records the banner-reported network in
`data/processed/cleavage/netchop_run_config.json`, and the container selftest
asserts that `-v 0` and `-v 1` produce different scores — so a silently ignored
flag fails loudly rather than producing plausible wrong numbers.

---

## 4. Reproducibility and licensing

**Container.** The shipped `netChop` is an ELF 32-bit i386 statically linked
Linux executable. It cannot run natively on Windows. The run is therefore
containerised (`docker/netchop/`), pinned to `debian:bookworm-20241111-slim`.

**Selftest.** The vendor ships `test/test.fsa` with `test/test.out`, but they do
not correspond: the input is Influenza A M2 (`BHB191648`, `MSLLTEVETPI…`) while
the expected output holds predictions for `gi|3333147` (`M`, `A`, `G`…), dated
2006. The fixture is unusable, so `docker/netchop/entrypoint.sh` verifies seven
properties instead — that predictions are produced, that the two networks differ,
that each banner names the network requested, that output is deterministic, that
scores lie in [0,1], and that the S/. flag agrees with the threshold.

**Licensing.** The NetChop distribution is licensed software (NetMHC 3.0 terms
plus Non-Profit OSL 3.0). It is **not** committed to this repository, which is
public; `.gitignore` excludes the tarball and any extracted binaries.
Redistribution would breach the licence.

Two items require a decision outside this workstream:

1. The vendor readme states the academic download site is for academic users and
   that "Other users are requested to contact software@cbs.dtu.dk." Non-Profit
   OSL is scoped by the *entity*, not by the research intent, so a commercial
   entity conducting research is not automatically covered.
2. The same licensing family covers NetMHCpan-4.1, whose **training data is
   already committed to this public repository** under
   `data/reference/netmhcpan41/`. That may be redistribution of licensed material
   and should be reviewed independently of this decision.

---

## Summary

| Question | Answer | Basis |
|---|---|---|
| Which tool? | NetChop 3.1 | NetCTLpan contains NetChop plus binding; cannot isolate processing |
| Which network? | 20S-3.0 | data-type match, constitutive labelling, verifiable HIV-disjointness |
| Which flag? | `-v 1` | binary man page and banner; IEDB's wrapper documents this backwards |
| Reproducible how? | pinned Docker image + 7-check selftest | binary is 32-bit Linux; vendor fixture is broken |
