# SparQ-HLA-Prediction

A research project investigating **HLA class I peptide-binding prediction** and **antigen processing** using classical machine-learning models and exploratory hybrid quantum machine-learning approaches.

The project is organised into two related but distinct research arms:

1. **HLA class I binding prediction** — predict peptide–HLA binding using experimentally measured binding data.
2. **Antigen processing** — investigate proteasomal cleavage and naturally presented HLA class I ligands.

The planned modelling sequence is: reference benchmark → in-house classical baseline → matched-capacity classical control → exploratory hybrid quantum model. Quantum-model results will be interpreted against matched-capacity controls, not just against a baseline with fewer parameters.

---

## Current Project Status

**Week 1 data preparation is complete.** The core binding, elution and cleavage datasets have been assembled, and sequence-level held-out datasets have been checked against the corresponding NetMHCpan-4.1 training references.

| Week 1 deliverable | Status | Evidence / output |
|---|---|---|
| Allele-stratified IEDB HLA class I binding dataset | Complete | `data/processed/binding/` |
| Measured-binder counts per allele | Complete | `results/tables/measured_binder_counts_by_allele.csv` and `final_binder_counts_by_allele.csv` |
| Cleavage and elution data assembled | Assembled; cleavage caveat recorded | `data/processed/cleavage/`, `data/processed/elution/` |
| Held-out overlap assessment against NetMHCpan-4.1 BA and EL training files | Complete for these references | `results/tables/netmhcpan_reference_disjoint_holdout_summary.csv` |

### Key Week 1 dataset figures

**IEDB binding**

- 38,395 binding observations were downloaded and processed.
- The allele summary covers 123 groups: 95 primary two-chain groups and 28 class-level groups.
- Peptide–allele pair classification and conflicts are retained in separate processing outputs. Conflicting observations are documented in `iedb_2chain_conflicts.csv`.

**IEDB eluted ligands**

The IEDB elution query returned **1,486,479 records** across HLA-A, HLA-B and HLA-C:

| Locus | Records |
|---|---:|
| HLA-A | 571,749 |
| HLA-B | 771,896 |
| HLA-C | 142,834 |
| **Total** | **1,486,479** |

All 151 compressed download chunks were read successfully; all 1,486,479 records parsed, with zero validation issues. The assembled observations are in `data/processed/elution/iedb_hla_i_elution_observations.csv.gz`.

**CALIS cleavage data**

- Supplementary Table S1 provides 16 HIV-1 GAG/TAT peptides, each represented under immunoproteasome (`I`) and constitutive proteasome (`C`) conditions.
- The extraction contains 32 peptide-condition patterns and 736 candidate-site labels.
- The table transcription passed structural checks with zero transcription errors.

**Important CALIS caveat:** the cleavage annotations counted directly from Supplementary Table S1 (I=140; C=139; 94 shared positions) do not reconcile with the totals reported in the article's main text. This discrepancy is recorded in `results/tables/calis_2015_cleavage_source_reconciliation.csv`. Keep the source-reconciliation report with the data. These labels should not be treated as fully reconciled ground truth until the discrepancy is resolved.

### NetMHCpan-4.1 reference-disjoint holdouts

NetMHCpan-4.1 BA and EL training files (`c000`–`c004`) were audited separately against the corresponding IEDB holdouts. A conservative sequence-level policy was then applied: **any holdout peptide sequence found in the corresponding NetMHCpan training files was excluded, regardless of whether the allele association was an exact or an ambiguous match.**

| Metric | Binding BA holdout | Elution EL holdout |
|---|---:|---:|
| Initial holdout rows | 6,502 | 297,040 |
| Rows excluded by sequence overlap | 3,203 | 205,225 |
| **Reference-disjoint rows retained** | **3,299** | **91,815** |
| Unique peptides retained | 2,198 | 54,579 |
| Residual shared peptide sequences | **0** | **0** |

The final reference-disjoint files are:

- `data/processed/binding/iedb_2chain_holdout_netmhcpan_disjoint.csv`
- `data/processed/elution/iedb_hla_i_elution_holdout_netmhcpan_disjoint.csv.gz`

The ordinary sequence-grouped split files are also retained. For the NetMHCpan-independent benchmark, use the **reference-disjoint** files above, not the original unfiltered holdouts.

**Scope of this audit:** the zero-overlap result applies to the corresponding NetMHCpan-4.1 BA and EL training files included in this repository. It does not establish independence from every possible external predictor or training corpus. The antigen-processing reference (NetChop or NetCTLpan) still requires its own final selection, data-compatibility review and overlap assessment.

---

## Week 2 Workload

The Week 2 objective in the project plan is to establish reference performance on the held-out data and report results by relevant allele or processing stratum rather than relying only on pooled metrics.

The two workstreams below are designed to run **in parallel**. Each owner has a separate dataset, reference-tool task and deliverable; neither workstream needs the other to finish first.

| Owner | Independent workstream | Main tasks | Friday deliverable |
|---|---|---|---|
| **Wasiu** | HLA binding reference benchmark | Run NetMHCpan-4.1 in BA mode against the reference-disjoint binding holdout; harmonise allele naming; collect predictions; calculate performance by allele/allele stratum; document unsupported alleles and failed predictions. | Per-allele NetMHCpan BA performance table, aggregate summary, run configuration and prediction outputs. |
| **Totan** | Antigen-processing reference benchmark | Finalise the processing reference choice (NetChop is currently preferred for a cleaner separation from the binding arm; NetCTLpan remains an alternative); define a biologically compatible input/target; review CALIS source limitations and assess reference-training overlap; run and summarise the agreed processing benchmark independently. | Processing-reference decision and rationale, input/target definition, overlap/compatibility note and stratified benchmark results where the available data support them. |

**Coordination point:** agree the processing reference and evaluation target at the Week 2 checkpoint. This is a decision for the processing workstream; it does not block Wasiu's NetMHCpan BA benchmark.

---

## Repository Structure

```text
SparQ-HLA-Prediction/
├── data/
│   ├── raw/
│   │   └── processing/
│   │       ├── cleavage/
│   │       │   └── calis_2015_supp_table_s1.pdf
│   │       └── elution/
│   │           ├── HLA_A_2_chain/
│   │           ├── HLA_B_2_chain/
│   │           ├── HLA_C_2_chain/
│   │           └── manifest.json
│   ├── reference/
│   │   └── netmhcpan41/
│   │       └── NetMHCpan_train/
│   │           ├── allelelist
│   │           ├── c000_ba ... c004_ba
│   │           └── c000_el ... c004_el
│   └── processed/
│       ├── binding/
│       ├── cleavage/
│       └── elution/
├── notebooks/
├── scripts/
├── src/
├── results/
│   └── tables/
└── docs/
```

The tree shows the principal project areas; some directories may contain additional files or be expanded as modelling work progresses.

### Important processed files

| Area | File | Purpose |
|---|---|---|
| Binding | `iedb_binding_observations.csv` | Processed IEDB binding observations |
| Binding | `iedb_2chain_training_candidates.csv` | Two-chain peptide–allele training candidates |
| Binding | `iedb_2chain_pair_status.csv` | Pair-level classification/status |
| Binding | `iedb_2chain_conflicts.csv` | Conflicting pair evidence |
| Binding | `measured_binder_counts_by_allele.csv` | Measured binder counts by allele |
| Elution | `iedb_hla_i_elution_observations.csv.gz` | Full validated elution observation set |
| Elution | `iedb_hla_i_elution_train.csv.gz` | Sequence-grouped training partition |
| Elution | `iedb_hla_i_elution_holdout_netmhcpan_disjoint.csv.gz` | NetMHCpan-EL sequence-disjoint holdout |
| Cleavage | `calis_2015_cleavage_patterns.csv` | Peptide-condition cleavage patterns |
| Cleavage | `calis_2015_cleavage_sites.csv` | Extracted cleavage-site annotations |
| Cleavage | `calis_2015_cleavage_site_labels.csv` | Candidate-site labels |
| Reports | `netmhcpan_reference_disjoint_holdout_summary.csv` | Final reference-disjoint overlap summary |
| Reports | `calis_2015_cleavage_source_reconciliation.csv` | CALIS source discrepancy documentation |

---

## Environment Setup

The Week 1 scripts were run locally with Python 3.13 in a virtual environment on Windows PowerShell.

From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install requests
```

`requests` is required by the IEDB API download scripts. The processing and audit scripts also use Python's standard library. If a later modelling component introduces additional dependencies, add them to a project dependency file and document the tested versions.

Do not commit `.venv/`, credentials, API tokens or other secrets.

---

## Reproducing the Week 1 Data Preparation

Run commands from the repository root. The raw IEDB/CALIS inputs and extracted NetMHCpan reference files are stored in the repository, so **download commands are not required** to reproduce processing from the committed inputs.

### 1. IEDB binding

```powershell
python scripts/validate_iedb_binding.py
python scripts/prepare_iedb_binding.py
python scripts/finalise_iedb_binding.py
python scripts/summarise_pair_conflicts.py
```

The binding acquisition and resolution-count scripts are available if the source data need to be refreshed:

```powershell
python scripts/count_iedb_resolution.py
python scripts/download_iedb_binding.py
```

### 2. IEDB elution

```powershell
python scripts/validate_and_prepare_iedb_elution.py
```

This validates the downloaded chunk files against the extraction manifest and creates the processed compressed observation file.

To refresh the source download, use:

```powershell
python scripts/download_iedb_elution.py
```

The small pilot inspection script is for query inspection only; it is not the full extraction or validation workflow.

### 3. CALIS cleavage data

```powershell
python scripts/prepare_calis_cleavage.py
```

Review both of the following after execution:

- `results/tables/calis_2015_cleavage_validation_summary.csv`
- `results/tables/calis_2015_cleavage_source_reconciliation.csv`

A structural validation pass does not resolve the published cleavage-count discrepancy.

### 4. NetMHCpan overlap audit and holdouts

```powershell
python scripts/inspect_netmhcpan_training.py
python scripts/create_holdouts_and_audit_netmhcpan.py
python scripts/create_netmhcpan_disjoint_holdouts.py
```

The first audit records exact, possible multi-allelic and sequence-only overlaps. The final script applies the conservative sequence-disjoint exclusion policy and writes the benchmark holdouts used for the NetMHCpan-independent evaluation.

---

## Data Sources and Provenance

### IEDB

- Source: [Immune Epitope Database](https://www.iedb.org/)
- Query API: `https://query-api.iedb.org/mhc_search`
- Binding and elution records were queried for linear peptides and HLA class I.
- Elution extraction covered HLA-A, HLA-B and HLA-C, with two-chain allele resolution and mass-spectrometry assay filtering.
- The elution manifest is stored at `data/raw/processing/elution/manifest.json`; it records the endpoint, extraction timestamp, filters, page size, chunk size and per-locus completion offsets.

### NetMHCpan-4.1

- Reference: [NetMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
- Training data are stored under `data/reference/netmhcpan41/NetMHCpan_train/`.
- BA and EL training files for configurations `c000` through `c004` and the `allelelist` mapping are included.
- The downloaded archive `NetMHCpan_train.tar.gz` is not tracked because the extracted training files are already stored in the repository. Re-download the archive from the official source if the extracted reference files need to be restored.

### CALIS et al. (2015)

- Source article: [DOI: 10.1007/s00251-014-0815-0](https://doi.org/10.1007/s00251-014-0815-0)
- Source data: Supplementary Table S1, stored at `data/raw/processing/cleavage/calis_2015_supp_table_s1.pdf`.
- The processed CSV files are derived from the cleavage bars shown in the supplementary table. The reconciliation report records the difference between the table transcription and totals stated in the article.

**Source-use note:** source datasets and tools may have separate access, citation and redistribution conditions. Check the relevant source terms before redistributing data or software outside this repository. This repository stores the NetMHCpan training reference files used in the audit; it does not include the NetMHCpan executable/software.

---

## Validation and Scientific Limitations

1. **Sequence-level leakage control:** the reference-disjoint holdouts were filtered against the corresponding NetMHCpan BA and EL training sequences. The final summaries report zero residual shared peptide sequences for those comparisons.
2. **Allele-specific versus sequence-only overlap:** exact peptide–allele matches, possible multi-allelic associations and sequence-only matches are distinct categories. The final conservative exclusion removes a sequence even when an exact allele assignment cannot be confirmed.
3. **CALIS source discrepancy:** the supplementary table counts (I=140; C=139; 94 shared positions) do not match the article's main-text counts. This needs source-level resolution before the labels are treated as definitive.
4. **Processing reference is a separate decision:** NetMHCpan BA/EL overlap results do not establish that a cleavage or antigen-processing reference has been audited. NetChop is currently preferred for a cleavage-focused comparison; NetCTLpan remains under consideration because it integrates more of the antigen-presentation pathway.
5. **Different data types answer different questions:** purified binding assays, eluted ligands and in-vitro cleavage annotations are not interchangeable labels. They must be evaluated with metrics and biological interpretations appropriate to each task.
6. **Allele-stratified reporting:** aggregate metrics alone can hide differences in allele support and performance. Report support counts and per-allele results alongside aggregate summaries.

---

## Project Workflow

1. Data collection and validation
2. Reference benchmarking
3. Classical model development
4. Matched-capacity controls
5. Hybrid quantum modelling
6. Evaluation and reporting

The current phase is reference benchmarking. Once both reference baselines and their evaluation protocols are established, the project proceeds to in-house classical models, then matched-capacity controls and exploratory hybrid quantum models.

---

## Scientific References

1. [IEDB — Immune Epitope Database](https://www.iedb.org/)
2. [IEDB MHC Class I Help](https://tools.iedb.org/mhci/help/)
3. [NetMHCpan-4.1 Technical Documentation](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
4. [NetMHCpan-4.1 Research Paper](https://doi.org/10.1093/nar/gkaa379)
5. [CALIS et al. (2015), source article](https://doi.org/10.1007/s00251-014-0815-0)
6. [NetCTLpan: Pan-specific MHC Class I Pathway Epitope Predictions](https://doi.org/10.1007/s00251-010-0441-4)

---

# Biological Background and Research Context

## 1. Introduction: What Are We Trying to Predict?

The human immune system continuously monitors cells for signs of infection, abnormal protein production and other cellular changes.

One important mechanism through which this monitoring occurs is the Human Leukocyte Antigen (HLA) system.

HLA molecules bind short fragments of proteins, known as peptides, and display them on the surface of cells. These displayed peptides can be examined by T cells, particularly CD8+ T cells in the case of HLA class I.

This project investigates two computational problems associated with this process:

1. HLA class I peptide binding prediction.
2. Antigen processing prediction.

The objective is to develop classical machine learning models and investigate whether hybrid quantum machine learning approaches offer any measurable advantage over appropriately matched classical controls.

## 2. Understanding HLA Class I

HLA stands for Human Leukocyte Antigen. It is the human version of the Major Histocompatibility Complex (MHC).

There are three major classical HLA class I gene groups:

- HLA-A
- HLA-B
- HLA-C

These genes are highly polymorphic, meaning that many different versions (alleles) exist across the human population.

Examples of HLA class I alleles include:

- HLA-A*02:01
- HLA-A*03:01
- HLA-B*07:02
- HLA-B*27:05
- HLA-C*07:01

Each allele has a distinctive peptide-binding groove. Consequently, different HLA alleles have different preferences for the amino acid sequences of peptides they bind.

This is an important consideration for our modelling work: a peptide that binds strongly to one HLA allele may bind weakly or not at all to another.

HLA class I molecules typically bind short peptides, commonly around 8–12 amino acids long, with 9-amino-acid peptides being particularly frequent.

Specific amino acid positions, known as anchor positions, contribute strongly to the interaction between a peptide and its HLA molecule.

These allele-specific sequence preferences are commonly referred to as peptide-binding motifs.

## 3. The HLA Class I Antigen Processing and Presentation Pathway

To understand our second research arm, it is important to understand how peptides become available for HLA class I presentation.

The process can be simplified into the following stages:

### Stage 1: Protein degradation

Proteins inside the cell, including viral proteins and normal cellular proteins, can be broken down into smaller fragments by the proteasome.

The proteasome is a large protein complex responsible for degrading proteins into shorter peptide fragments.

### Stage 2: Peptide transport

Some generated peptides are transported from the cytoplasm into the endoplasmic reticulum (ER).

This transport is facilitated by a protein complex called Transporter Associated with Antigen Processing (TAP), consisting of TAP1 and TAP2.

### Stage 3: Peptide trimming

Once inside the ER, some peptides may undergo additional processing by enzymes such as ERAP1 and ERAP2.

This trimming helps generate peptides of suitable lengths for HLA class I binding.

### Stage 4: HLA class I peptide loading

The processed peptides are loaded onto HLA class I molecules.

The HLA class I molecule associates with beta-2-microglobulin (ß2m) and other components of the peptide-loading machinery during assembly.

### Stage 5: Cell surface presentation

Successfully loaded HLA class I molecules travel through the secretory pathway to the cell surface.

The peptide-HLA complex is then available for immune surveillance by T cells.

CD8+ T cells can recognise specific peptide-HLA class I combinations through their T-cell receptors.

Importantly, successful peptide binding does not automatically mean that a peptide will be recognised by a T cell.

## 4. Our Two Research Arms

Although HLA binding and antigen processing are related, they represent different biological questions.

### Research Arm 1: HLA Class I Binding Prediction

**Biological question:**

Given a peptide sequence and an HLA class I allele, how likely is the peptide to bind to that HLA molecule, and what is its measured or predicted binding strength?

Inputs:
- Amino acid peptide sequence.
- HLA class I allele.

Expected outputs:
- Predicted binding score or affinity.
- Binder/non-binder classification where appropriate.

For this arm, experimentally measured HLA binding data will be obtained from the Immune Epitope Database (IEDB).

The reference benchmark is NetMHCpan-4.1, using its binding affinity (BA) prediction mode.

The in-house classical model will learn the relationship between peptide sequence, HLA allele and experimentally measured binding outcomes.

An important requirement is to evaluate results separately across HLA alleles rather than relying exclusively on pooled performance.

### Research Arm 2: Antigen Processing Prediction

**Biological question:**

How do intracellular proteins produce peptide fragments that can subsequently become available for HLA class I presentation?

This research arm focuses on the upstream processing of proteins into peptides and evidence of naturally presented peptides.

Relevant biological components include:

- Proteasomal cleavage: how proteins are broken into peptide fragments.
- TAP transport: movement of peptides into the endoplasmic reticulum.
- Peptide trimming: additional processing of peptide fragments.
- HLA eluted ligands: peptides experimentally identified from HLA molecules, commonly through mass spectrometry.

The planned dataset will include cleavage and elution information.

The processing reference is still subject to final selection:

- NetChop focuses on proteasomal cleavage prediction.
- NetCTLpan combines proteasomal cleavage, TAP transport and HLA class I binding predictions.

NetChop is currently preferred for a cleaner separation from the HLA binding research arm.

The precise processing prediction target and evaluation dataset will be documented once the datasets and reference-tool scope have been finalised.

## 5. Understanding the Biological Data

The main biological data types in this project are:

### A. Peptide sequences

Peptides are represented using the one-letter amino acid alphabet.

For example:

SIINFEKL

This represents a peptide containing eight amino acids.

The sequence is not simply an arbitrary string. Each character represents a specific amino acid with chemical and structural properties.

### B. HLA alleles

An HLA allele identifies a particular version of an HLA molecule.

For example:

HLA-A*02:01

This is important because the peptide-binding preferences depend on the HLA allele.

### C. Binding affinity measurements

Experimental binding assays provide quantitative information about the interaction between a peptide and an HLA molecule.

Common measurements include:

- IC50: concentration associated with 50% inhibition in a binding assay.
- KD: equilibrium dissociation constant.

For IC50 and KD measurements, lower values generally indicate stronger binding.

However, different assay types and units must be handled carefully during data preparation.

### D. Eluted ligands

Eluted ligands are peptides experimentally recovered from HLA molecules, commonly using mass spectrometry.

They provide evidence of peptides that have naturally reached and associated with HLA molecules.

Eluted ligand observations are not identical to purified-protein binding affinity measurements. They reflect additional biological steps involved in natural antigen presentation.

## 6. How the Computational Models Fit Into the Biology

The project translates biological sequences and measurements into computational learning problems.

For the binding arm, the model will learn relationships between:

- Peptide sequence.
- HLA allele identity.
- Experimentally measured binding outcomes.

For the processing arm, the model will investigate relationships between protein-derived sequences, cleavage evidence and peptide presentation observations.

These relationships can be represented mathematically and explored using classical machine learning techniques.

The subsequent hybrid quantum arm will investigate whether incorporating quantum computational layers offers any measurable improvement.

To make the comparison meaningful, each quantum layer will be compared with a classical layer having a matched parameter count in the same architectural position.

This is necessary to distinguish any observed difference associated with the quantum approach from improvements that could simply result from additional model capacity.

## 7. Important Scientific Distinctions

The following concepts should not be treated as interchangeable:

| Concept | Meaning |
|---|---|
| Antigen | A molecule or protein that can be recognised by the immune system. |
| Peptide | A short chain of amino acids, often generated from a protein. |
| HLA binding | Interaction between a peptide and an HLA molecule. |
| Antigen processing | Biological steps involved in generating and preparing peptides for presentation. |
| Antigen presentation | Display of peptide-HLA complexes on the cell surface. |
| Immunogenicity | Ability of an antigen or peptide to trigger an immune response. |
| Epitope | A specific molecular region recognised by an immune receptor or antibody. |

**Important:** This project investigates HLA binding and antigen processing. It does not directly establish whether a peptide will trigger a T-cell response or predict clinical immunotherapy outcomes.

## 8. Why Allele-Stratified Evaluation Matters

HLA alleles differ in their peptide-binding preferences.

A model that performs well on one allele may perform poorly on another.

Consequently, reporting only an overall accuracy or performance score may conceal important differences.

Our evaluation will therefore report performance separately across HLA allele groups, alongside relevant aggregate results.

The project will also examine dataset overlap with the training data of the reference tools to minimise information leakage and avoid misleading benchmark results.

## 9. Key Terminology

| Term | Definition |
|---|---|
| HLA | Human Leukocyte Antigen |
| MHC | Major Histocompatibility Complex |
| TAP | Transporter Associated with Antigen Processing |
| ER | Endoplasmic Reticulum |
| ERAP | Endoplasmic Reticulum Aminopeptidase |
| ß2m | Beta-2-microglobulin |
| IEDB | Immune Epitope Database |
| BA | Binding Affinity |
| EL | Eluted Ligand |
| IC50 | Half-maximal inhibitory concentration |
| CD8+ T cell | Immune cell involved in recognising peptide-HLA class I complexes |
| Proteasome | Cellular protein-degrading complex |
| Held-out dataset | Data reserved for evaluation rather than model training |

## 10. Biological and Technical References

1. NCBI Bookshelf - Major Histocompatibility Complex: Antigen Processing and Presentation.
   https://www.ncbi.nlm.nih.gov/books/NBK459467/

2. Immune Epitope Database (IEDB) - MHC Class I Help.
   https://tools.iedb.org/mhci/help/

3. NetMHCpan-4.1 - Technical Documentation.
   https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

4. NetMHCpan-4.1 Research Paper.
   https://doi.org/10.1093/nar/gkaa379

5. NetCTLpan: Pan-specific MHC Class I Pathway Epitope Predictions.
   https://doi.org/10.1007/s00251-010-0441-4
