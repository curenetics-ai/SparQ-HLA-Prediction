# Week 2 GitHub update — NetMHCpan 4.1b BA benchmark

This update adds the completed Wasiu Week 2 binding-reference benchmark to the SparQ-HLA-Prediction repository.

## Included
- All eight original DTU NetMHCpan 4.1b result PDFs.
- Eight reproducible submission input folders with peptide lists, allele lists and exact target-pair manifests.
- PDF extraction and evaluation script.
- All 47,286 extracted server prediction rows (compressed CSV).
- Exact 3,285-row reference-disjoint target-pair prediction table.
- Per-allele, per-locus and pooled metric tables.
- Job-level and aggregate QC, unsupported/failed prediction reports, length eligibility report, allele nearest-neighbour mapping and run manifests.
- Benchmark results documentation and updated project README.

## Key result
3,285 / 3,285 eligible target peptide–allele pairs were matched (100% coverage). NetMHCpan 4.1b BA pooled ROC-AUC = 0.8456 and average precision = 0.8727. See `docs/NetMHCpan41_BA_Benchmark_Results.md` and `results/tables/` for full allele-stratified results and caveats.

## Re-run extraction and evaluation
From the repository root, with Python dependencies installed:

```powershell
python -m pip install pymupdf scikit-learn
python scripts/extract_and_evaluate_netmhcpan41_pdfs.py
```

The run uses the eight PDFs and target-pair manifests included in this update. It reads PDF text rather than OCR.

## Important interpretation notes
- Fourteen of the initial 3,299 reference-disjoint pairs were below the DTU server's accepted peptide-length minimum and remain documented separately.
- The evaluation is based on BA-specific `%Rank_BA`; the server `BindLevel` field follows EL rank categorisation and must not be used for BA labels.
- Fourteen requested alleles were evaluated using a different nearest-neighbour allele; eight reported non-zero distance to training data.
- AUC/AP cannot be estimated for 20 alleles because each has only one experimental class.
- These are NetMHCpan reference-predictor results, not the in-house SparQ classical or quantum model results.
