# R7 statistical reproduction

From the repository root, use a Python environment configured with the pinned requirements. Point the command to the folder containing the extracted deposit's `data/` directory. Choose a new or empty output directory.

```sh
python scripts/r7_stats/reproduce.py \
  --data-dir /path/to/ASD_TREG_CRISPR_processed_data_v1.1.0 \
  --out-dir ./r7-rebuilt
```

Add `--verify-r` to verify the primary CPM ≥1 M2 fits with the independent base-R objective. The Python pipeline uses the frozen Jeffreys penalized likelihood and tests the same-design profiles. The two primary M2 P values share a Holm family. M3 estimates are secondary sensitivities. The expected-results fixture records numerical tolerances and both out-of-support and non-estimable outcomes.

The deposit supplies accepted canonical gene identities, screen-native FDR and log-fold changes, author-derived per-gene expression summaries and threshold votes, and frozen property mappings and LoF eligibility. Starting from those inputs, the entry point reconstructs detection states and stratum counts, runs all 15 Fisher backgrounds and 30 exact conditional comparisons, refits 33 property models and the R6 bridge, and reruns fixed-caliper matching and candidate direction checks. Gene-level detection and result tables are written into the requested output directory. See the project reproducibility documentation for primary-source retrieval and upstream preprocessing before harmonization.

The validation table `checks/EXPECTED_RESULTS_QA.tsv` compares the computed results to `config/r7_expected_results.tsv`. It uses exact comparison for counts and the listed tolerances for fitted or profiled values. Attempts that use sources outside the extracted deposit are rejected by the input resolver.
