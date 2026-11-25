# Functional comparisons

The [release guide](../../docs/REPRODUCIBILITY.md) gives the main replay command. This branch can also run separately with an extracted data root and a new output directory.

```sh
bash scripts/functional_comparisons/run.sh DATA OUT
```

The analysis reconstructs directional loss-of-function enrichment and the shared FOXP3–TCR proliferation comparison. The [pinned inputs](frozen_inputs.tsv) give the source hashes. [Methods](METHODS.md) define the exact strata, discordant Firth sensitivity and correction families. The parity report under `checks/` compares all released tables with relative tolerance 1e-7 and absolute tolerance 1e-10 for platform roundoff. Independent exact and base-R Firth checks are also written under `checks/`.
