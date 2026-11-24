# R7 evidence and figure-source verification

From the repository root, pass the unpacked Zenodo archive's versioned R7 data directory. The verifier uses the Python standard library.

```bash
python3 -I scripts/r7_evidence/verify_release_evidence.py \
  path/to/ASD_TREG_CRISPR_processed_data_v1.1.0/data/r7 \
  --report reproduced/r7_evidence_checks.tsv
```

The verifier reconciles 26 Figure 3d source records to the frozen enrichment, strict conditioning, common-case regression and historical bridge tables. It checks the complete-case margins shared by M0–M3 and locates Holm-corrected primary M2 values. The exact LoF104 sensitivity is recalculated from the released stratum margins by hypergeometric convolution. Scope-specific case counts use the released per-gene and case-support records.

For Supplementary Figure S6, it checks 245 score-1 memberships, 104 LoF-compatible genes, all 17 candidate rows, the measured and unmeasured screen contexts, source-native LFC signs, standardized KO direction, and the published experiments associated with each assay label. It independently groups endpoint records into experiments using the frozen design table and recovers each assay category's display labels. The 104 eligible genes include 100 inferred and four G2P evidence-labelled genes with mixed ASD and neurodevelopmental scope.

The verification starts from released processed data. Original artwork geometry and rebuilding those processed inputs from their underlying public or separately hosted resources require the upstream files listed in the release provenance. The checks concern the numeric figure sources and evidence scope. Results from the four published studies refer to their reported experiments and readouts. Their identifiers and source figure/table locators are provided in `function/VERIFIED_PUBLICATION_IDENTIFIERS.tsv` and the evidence tables.
