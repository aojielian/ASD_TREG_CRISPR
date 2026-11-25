# Experimental evidence and source verification

The [release workflow](../../docs/REPRODUCIBILITY.md) runs the processed-data verifier. It can also run directly with an extracted data root.

```sh
python3 -I scripts/experimental_evidence/verify_release_evidence.py DATA --report OUT/evidence.tsv
```

The verifier reconciles 26 historical source rows to enrichment, conditioning and model results. Its independent LoF exact test checks the frozen 104-member scope sensitivity. Candidate checks span 245 score-1 memberships, 17 candidates, source-native LFC and standardized KO directions, measurement availability, and published experiment sources and denominators. Original author workbooks and historical figure underlays have their own source terms and acquisition instructions.
