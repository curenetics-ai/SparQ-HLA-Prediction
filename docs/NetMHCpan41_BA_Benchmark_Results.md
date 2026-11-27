# NetMHCpan 4.1b Binding-Affinity Benchmark Results

## Scope

This package consolidates the eight DTU NetMHCpan web-server PDF outputs and matches predictions to the exact peptide–allele target pairs in the corresponding batch `target_pairs.csv` files.

The PDFs report **NetMHCpan version 4.1b**, PEPTIDE input format, and both EL and BA outputs. The primary binding benchmark uses BA predictions.

## Result summary

| Item | Result |
|---|---:|
| DTU PDF jobs processed | 8 |
| PDF pages parsed | 515 |
| Allele-job sections in PDFs | 138 |
| Distinct HLA alleles represented | 78 |
| Full server cross-product rows parsed | 47,286 |
| Reference-disjoint holdout rows before length eligibility | 3,299 |
| Rows eligible and submitted | 3,285 |
| Unique peptides in scored targets | 2,184 |
| Target pairs matched to predictions | 3,285 |
| Target prediction coverage | 100.0% |
| Failed/missing eligible target predictions | 0 |
| Unsupported target alleles with no returned section | 0 |

Fourteen rows from the original 3,299-row holdout were not submitted because the peptide was shorter than the DTU server's accepted minimum length. They are not failures and have not been deleted from the original holdout. See `netmhcpan41_ba_length_eligibility_summary.csv`.

## Performance results

The continuous prediction score is **negative `%Rank_BA`** (lower percentile rank indicates stronger predicted binding). ROC-AUC and average precision are calculated against the experimental `binary_label` supplied in the target-pair files.

| Locus | Alleles | Target pairs | Positive | Negative | ROC-AUC | Average precision |
|---|---:|---:|---:|---:|---:|---:|
| HLA-A | 23 | 2,233 | 1,301 | 932 | 0.8564 | 0.8873 |
| HLA-B | 39 | 745 | 378 | 367 | 0.8402 | 0.8492 |
| HLA-C | 16 | 307 | 158 | 149 | 0.8097 | 0.8456 |
| All loci (pooled) | 78 | 3,285 | 1,837 | 1,448 | 0.8456 | 0.8727 |

### BA rank cutoffs

The server headers state a strong-binder rank threshold of 0.5% and a weak-binder rank threshold of 2.0%. BA categories in the target-pair table are therefore derived from `%Rank_BA`:
- Strong: `%Rank_BA <= 0.5`
- Weak: `0.5 < %Rank_BA <= 2.0`
- Non-binder: `%Rank_BA > 2.0`

At the combined weak-or-strong cutoff (`%Rank_BA <= 2.0`), pooled sensitivity is 0.668, specificity is 0.850, and balanced accuracy is 0.759. The allele-level and locus-level files contain the corresponding metrics and confusion counts.

**Important server-output detail:** the PDF `BindLevel` field matched the EL rank categories on all 3,285 target rows, not the BA rank categories. Do not use `BindLevel` as the BA predicted class. Use the BA percentile-rank thresholds above.

## Allele mapping / nearest-neighbour information

All 78 requested alleles returned prediction sections. In 14 cases, the reported nearest-neighbour allele differed from the requested allele; 8 alleles had a non-zero distance to training data. Preserve these fields when interpreting the per-allele results; they are recorded in `netmhcpan41_ba_allele_nearest_neighbor_mapping.csv`.

A per-allele AUC is estimable for 58 of 78 alleles because both experimental classes are represented. For 20 alleles, only one class is present, so ROC-AUC and average precision are not estimable. 58 alleles have fewer than 20 target observations; treat their estimates as descriptive only. 20 alleles have at least 20 observations and 14 have at least 50.

## Files

- `data/raw/reference_predictions/netmhcpan41_ba/`: all eight original DTU result PDFs.
- `data/interim/netmhcpan41_ba_submissions/`: the eight batch target-pair manifests plus peptide and allele inputs.
- `data/processed/reference_predictions/netmhcpan41_ba_el_all_server_pairs.csv.gz`: all 47,286 extracted server output rows (including peptide-by-allele cross-product rows not used in the evaluation).
- `data/processed/binding/netmhcpan41_ba_target_pair_predictions.csv`: the exact 3,285 target rows matched to their server predictions and experimental labels.
- `results/tables/netmhcpan41_ba_per_allele_metrics.csv`: per-allele sample counts, coverage, AUC/AP, rank-threshold metrics and nearest-neighbour fields.
- `results/tables/netmhcpan41_ba_per_locus_metrics.csv`: HLA-A/B/C summaries.
- `results/tables/netmhcpan41_ba_overall_metrics.csv`: pooled supporting summary.
- `results/tables/netmhcpan41_ba_allele_nearest_neighbor_mapping.csv`: server mapping and distance fields.
- `results/tables/netmhcpan41_ba_per_job_qc.csv`: job-by-job extraction and target-match checks.
- `results/tables/netmhcpan41_ba_extraction_qc.csv`: overall extraction and match checks.
- `results/tables/netmhcpan41_ba_failed_predictions.csv`: empty report; no eligible target predictions failed or were missing.
- `results/tables/netmhcpan41_ba_unsupported_alleles.csv`: empty report; all submitted target alleles returned prediction sections.
- `results/tables/netmhcpan41_ba_benchmark_manifest.json`: setup and run provenance.

## Reproduction

From the repository root, after placing the eight PDFs and the eight batch input folders in the paths shown above:

```powershell
python -m pip install pymupdf scikit-learn
python scripts/extract_and_evaluate_netmhcpan41_pdfs.py
```

The script reads PDF text (no OCR), validates each batch's target pairs against its corresponding PDF, and regenerates the extracted prediction and metric tables.

## Interpretation limitations

1. The evaluation is for the **3,285 server-length-eligible reference-disjoint target pairs**, not all 3,299 rows in the preceding holdout.
2. Per-allele scores with small sample sizes are uncertain. AUC/AP are not reported where one experimental class is absent.
3. Fourteen requested alleles were evaluated through a different nearest-neighbour allele according to the DTU PDF header. Eight of the 14 also report non-zero distance to training data. Preserve this context in interpretation.
4. The hosted server model parameter files are not contained in these PDFs. The repository's `NetMHCpan_train` files are reference training datasets, not the executable predictor or installed model parameter package.
5. This result is a reference-predictor benchmark; it is not the performance of the in-house SparQ classical model or the hybrid quantum model.
