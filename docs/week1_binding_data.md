\# Week 1: HLA Class I Binding Dataset



\## 1. Dataset Source



\- Database: Immune Epitope Database (IEDB)

\- Website: https://iedb.org/

\- API: https://query-api.iedb.org/

\- Endpoint: mhc\_search

\- Extraction date: 19 November 2026

\- Data type: Experimentally measured HLA class I binding assays



\## 2. Extraction Criteria



The following filters were applied:



\- Peptide structure: Linear peptide

\- MHC class: Class I

\- HLA loci: HLA-A, HLA-B and HLA-C

\- Assay: Binding-related assays

\- Allele resolutions: 2 chain and class



The extraction used 500-record API pages and included

nested MHC export records.



\## 3. Dataset Size



| HLA locus | 2-chain | Class-level | Total |

|---|---:|---:|---:|

| HLA-A | 22,503 | 2,652 | 25,155 |

| HLA-B | 10,565 | 851 | 11,416 |

| HLA-C | 1,804 | 20 | 1,824 |

| Total | 34,872 | 3,523 | 38,395 |



\## 4. Raw Data Validation



\- Total records: 38,395

\- Unique record identifiers: 38,395

\- Nested assay records: 38,395

\- Empty assay exports: 0

\- Structural validation errors: 0



\## 5. Primary Allele-Specific Dataset



The two-chain resolution dataset contains:



\- 34,872 experimental records

\- 95 distinct HLA allele groups

\- 32,654 unique peptide-allele pairs



Pair-level classification:



| Classification | Count |

|---|---:|

| Unambiguous positive | 11,796 |

| Unambiguous negative | 20,648 |

| Conflicting | 210 |

| Unclassified | 0 |



Conflicting observations are retained separately

and excluded from initial single-label training candidates.



\## 6. Processing Decisions



\- Positive, Positive-Low, Positive-Intermediate

&#x20; and Positive-High are grouped as positive.

\- Negative observations are grouped as negative.

\- Original qualitative outcomes are preserved.

\- Conflicting peptide-allele pairs are flagged.

\- Quantitative values, units and measurement

&#x20; qualifiers are preserved.

\- Different measurement types are not assumed

&#x20; to be interchangeable.



\## 7. Generated Files



\### Raw data

data/raw/binding/



\### Processed data

\- iedb\_binding\_observations.csv

\- iedb\_2chain\_pair\_status.csv

\- iedb\_2chain\_training\_candidates.csv



\### Results

\- iedb\_raw\_validation.csv

\- iedb\_2chain\_conflicts.csv

\- iedb\_2chain\_pair\_summary\_by\_allele.csv

\- final\_binder\_counts\_by\_allele.csv



\## 8. Reproducibility



Scripts used for retrieval and preparation are

maintained in the scripts/ directory.



The original API records are retained separately

from derived datasets.



\## 9. Data Attribution



Source: Immune Epitope Database (IEDB).



Please acknowledge IEDB and cite the original

experimental publications associated with the

individual observations.



Website: https://iedb.org/

