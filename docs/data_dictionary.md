\# Data Dictionary



\## Primary observation dataset



File: iedb\_binding\_observations.csv



| Column | Description |

|---|---|

| elution\_id | IEDB search record identifier |

| assay\_id | Experimental assay identifier |

| peptide | Amino acid sequence |

| peptide\_length | Number of amino acids |

| hla\_allele | HLA allele designation |

| hla\_locus | HLA-A, HLA-B or HLA-C |

| mhc\_class | MHC class |

| allele\_resolution | IEDB allele resolution |

| qualitative\_outcome | Original experimental outcome |

| binary\_outcome | Grouped positive/negative outcome |

| assay\_response | Type of response measured |

| assay\_method | Experimental assay methodology |

| quantitative\_measurement | Reported numerical measurement |

| measurement\_units | Original measurement units |

| measurement\_inequality | Measurement inequality |

| reference\_id | IEDB publication reference identifier |

| pmid | PubMed identifier |



\## Pair-level dataset



File: iedb\_2chain\_pair\_status.csv



Each row represents a unique peptide-HLA allele

combination.



Important fields:



\- hla\_allele

\- peptide

\- positive\_assays

\- negative\_assays

\- other\_assays

\- classification

\- assay\_ids

\- pmids

\- evidence



\## Training candidate dataset



File: iedb\_2chain\_training\_candidates.csv



Contains unambiguous positive and negative

peptide-allele pairs.



binary\_label:

\- 1 = Positive

\- 0 = Negative



Conflicting observations are excluded.



\## Important distinction



An assay observation is not necessarily a unique

peptide-allele pair.



Multiple experiments can investigate the same

peptide and HLA allele.



Original measurements and experimental evidence

should be preserved during downstream analysis.

