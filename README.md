# SparQ-HLA-Prediction

## Project Overview
This project investigates HLA class I binding prediction and antigen processing using classical and exploratory hybrid quantum machine learning approaches.

## Research Arms
1. HLA Class I Binding Prediction
2. Antigen Processing

## Project Workflow
1. Data collection and validation
2. Reference benchmarking
3. Classical model development
4. Matched-capacity controls
5. Hybrid quantum modelling
6. Evaluation and reporting

## Repository Structure
- data/: Dataset storage
- notebooks/: Exploratory analysis
- scripts/: Data processing scripts
- src/: Model implementation
- results/: Experimental outputs
- docs/: Research documentation

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

