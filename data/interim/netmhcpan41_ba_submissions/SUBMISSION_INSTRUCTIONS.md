# NetMHCpan-4.1 BA submission instructions

Server: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

Only rows with peptide lengths from 8 to
20000 residues are included in the submission batches.
Shorter/longer rows are preserved in
`results/tables/netmhcpan41_length_ineligible_target_pairs.csv`.

For each job directory under `data/interim/netmhcpan41_ba_submissions/`:
1. Open the DTU NetMHCpan-4.1 server.
2. Select input type **PEPTIDE**.
3. Upload that job's `peptides.txt` (one peptide per line).
4. Copy the comma-separated allele names from `alleles.txt` into the allele-name field.
5. Enable **Include BA predictions**.
6. Set the %Rank filtering threshold to **-99** to print all predictions.
7. Keep strong/weak annotation thresholds at **0.5** and **2.0** percent.
8. Submit the job and save the complete raw result to the `save_raw_output_as`
   path recorded in `results/tables/netmhcpan41_submission_batch_manifest.csv`.
9. Retain the raw header because it records the predictor version and selected mode.

The web service allows at most 500 peptide sequences and
20 MHC alleles per submission. It requires peptides to be
8–20000 residues long.

The server predicts all peptide-by-allele combinations within a job. Each job's
`target_pairs.csv` is the authoritative list of benchmark pairs. During result
consolidation, retain only these exact target pairs.

For binding evaluation use the BA-specific percentile rank output (`%Rank_BA`)
as the primary score. Retain BA score and predicted affinity as secondary outputs.
Do not substitute EL score or EL rank for BA output.
