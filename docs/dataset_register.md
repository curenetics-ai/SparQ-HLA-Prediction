# Dataset Register

## Dataset 1: IEDB HLA Class I Binding

* Research arm: HLA class I binding prediction
* Source: Immune Epitope Database (IEDB)
* Website: https://iedb.org/
* Access date: 2026-11-18
* Dataset type: Experimentally measured MHC binding assays
* Planned output: TSV
* Raw storage: data/raw/binding/
* Query API URL: Pending
* Total records: Pending
* Number of unique HLA alleles: Pending
* Dataset version/date: Pending

### Initial Query Criteria

* Epitope: Linear peptide
* Assay category: MHC Ligand
* Assay type: Binding
* MHC restriction: Class I
* Outcomes: Positive and negative
* Peptide length: Unrestricted initially

### Notes

Binding affinity measurements and assay metadata
will be retained for subsequent quality control.

Eluted ligand observations must not be mixed into
this binding dataset.

The final human HLA allele selection will be
documented during data cleaning.

## Dataset 2: Antigen Processing

* Sources: Pending
* Dataset type: Proteasomal cleavage and HLA eluted ligands
* Raw storage: data/raw/processing/
* Download status: Pending



\### Preliminary IEDB Record Counts



Extraction date: 2026-11-19



| HLA locus | Preliminary records |

|---|---:|

| HLA-A | 43,124 |

| HLA-B | 52,231 |

| HLA-C | 8,514 |

| Total | 103,869 |



Count method: PostgREST planned count.



Status: Preliminary estimates, pending full extraction

and validation.



These counts represent assay records, not unique

peptides or confirmed positive binders.

