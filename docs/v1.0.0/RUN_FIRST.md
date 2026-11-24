# Run first

1. Download `ASD_TREG_CRISPR_processed_data_v1.0.0.zip` from Zenodo: [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659), then extract it.
2. From this repository, run `bash scripts/reproduce_from_processed/run_all_processed.sh path/to/ASD_TREG_CRISPR_processed_data_v1.0.0`.
3. Confirm `MINIMAL_REPRODUCTION_STATUS=PASS` and `EXPECTED_RESULTS_FAILURES=0`.

This verifies the processed-data workflow only; it is not a raw-source rebuild.
