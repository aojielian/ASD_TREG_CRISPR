# Mouse RNA-seq replay

```sh
bash scripts/conditional_knockout/run.sh /path/to/extracted/processed/archive /path/to/new/output
```

The entry point fits the frozen DESeq2 models from the three original GEO count matrices, runs edgeR sensitivity fits, repeats independent negative-binomial checks, calculates signature coverage, per-library ranks, cohort models, corrections and CAMERA sensitivity, and integrates candidate-specific human and mouse endpoints. It compares all resulting tables to the released results and writes the statistical parity report under `07_qa/`. The prepared GEO sample and gene identifier maps, prespecified signature memberships and closed literature curation are versioned inputs. Counts are KO minus control. The two Foxp1 models include the genotype-by-state interaction and the Sin3a model adjusts for sex.

`inputs.tsv` provides the required provenance and SHA-256 digests. The four main GEO files are the three source gene-count matrices and the GSE121279 authors' reported result table. The optional checks below redo the GEO metadata and signature preparation. Retrieve the official files marked `external` in the `origin` column of `inputs.tsv` and place them at the indicated paths within `EXTERNAL_INPUT_ROOT`. The MGI report is an upstream snapshot that can change in place. The current report can be used when its SHA-256 matches the listed snapshot.

```sh
bash scripts/conditional_knockout/run.sh DATA OUT EXTERNAL_INPUT_ROOT --rebuild-samples --rebuild-signatures --rebuild-literature
```

The optional flags can run separately. GEO SOFT files are used only when reconstructing sample mappings. The signature reconstruction reads exact MSigDB 2025.1 reference files and the frozen ImmGen expression reference. The literature reconstruction joins the frozen citations and curation to the original Chen workbooks, rechecks the author's protein bins and scores, and regenerates the literature audit, screening and summary tables. It uses the original bounded search results, including indirect and unavailable evidence. Direct RNA abundance, computed signature effects, human protein-bin screen entities and publication-derived assay summaries keep their original units and provenance.

The [release guide](../../docs/REPRODUCIBILITY.md) describes the integrated workflow and environment. Select an installed Python interpreter with `PYTHON` and an installed Rscript with `RSCRIPT` if needed. The R processes run with `--vanilla`. Output paths must be new or empty. In the tested environment, these commands use R 4.5.1, DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, org.Mm.eg.db 3.22.0, AnnotationDbi 1.72.0, and Python 3.13. The primary replay verifies the prepared maps and signature inputs by exact digest. Numeric comparison with other compatible runtime builds retains full categorical and missing-value agreement and records the maximum absolute difference when bytes differ.
