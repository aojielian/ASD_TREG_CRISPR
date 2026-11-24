# Reproduce the release

The complete processed-data workflow uses this repository and the [v1.1.0 data archive](https://doi.org/10.5281/zenodo.24487407). Use CPython 3.13 and base R with `Rscript` on the executable path. The environment was tested with CPython 3.13.2 and R 4.5.1. Use a local environment and the pinned dependency files.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -r environment/r7_figure_requirements.txt
bash scripts/smoke_test.sh ../ASD_TREG_CRISPR_processed_data_v1.1.0.zip release-check
```

Paths are relative to the repository root in the examples. Replace the data argument with a relative path to the archive or its extracted root. Specify a new or empty output directory. Verification checks the archival SHA-256 sidecar when it accompanies the ZIP, then the active data checksums, manifest inventory and table schemas. It runs the legacy results, the frozen R7 background and model comparisons, source evidence checks, and the ten vector PDF / 600 dpi TIFF renders. Base R independently checks the two primary M2 profiles and the legacy Loo conditional estimate. The results and tolerances are in `expected_results/EXPECTED_CORE_RESULTS.tsv`, `config/r7_expected_results.tsv` and `expected_results.tsv`. The R checks and Python model grid use the stated profile and convergence comparison rules.

The reproduced outputs include historical numeric checks under `legacy/`, adjusted coefficients and sensitivities under `r7/`, source evidence verification and figures. Comparison of statistical results uses the published tolerances. File digests identify the shipped source tables. PDF and TIFF are deterministic within the pinned rendering environment, and font or platform builds can affect byte-level similarity on other systems.

## Individual components

With an extracted archive, define a path relative to the repository root:

```sh
DATA=../ASD_TREG_CRISPR_processed_data_v1.1.0
python scripts/legacy/run_legacy_reproduction.py "$DATA" legacy-check --rscript Rscript
python scripts/r7_stats/reproduce.py --data-dir "$DATA" --out-dir r7-check --verify-r
python scripts/zhu/extract_zhu_evidence.py verify --derived-dir "$DATA/data/r8/zhu_derived" --candidate-master "$DATA/data/r7/mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv"
python figures/render_main.py "$DATA" rendered
python figures/render_supp.py "$DATA" rendered
```

The main renderer writes four `Main_Figure_*.pdf` / `.tif` pairs. The supplementary renderer writes six individual pairs under `individual_supplementary/` and a combined `Supplementary_Figures_FINAL.pdf`. The [figure source map](FIGURE_SOURCE_MAP.md) links each panel to the particular source rows and transformations. The separate historical figure scripts retain their original inputs and render earlier presentation versions.

## Optional reconstruction from provider sources

The processed workflow needs no additional scientific inputs. Original count normalization and Zhu derivative extraction require the exact independently obtained files described in the [source guide](SOURCES.md) and companion `SOURCE_ACCESS_AND_REUSE.md`. These files are outside the public dataset archive. The full upstream reconstruction also depends on the frozen versions of provider databases and the cited manual evidence-scope adjudication.

```sh
python scripts/upstream/reconstruct_expression.py "$DATA" \
  --itreg-counts provider/count_summary.txt \
  --peripheral-counts provider/RNAseq_salmon.merged.gene_counts.tsv
python scripts/zhu/extract_zhu_evidence.py extract \
  --source-dir provider/zhu-originals \
  --candidate-master "$DATA/data/r7/mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv" \
  --out-dir reconstructed-zhu
```

Expression reconstruction accepts optional `--hpa provider/rna_immune_cell.tsv` when the historical HPA file is available. The source scripts verify frozen file hashes before rebuilding their derivative fields. When a provider supplies a later snapshot, select the recorded version or a lawful hash-matching archive. Keep original third-party files subject to their separate reuse terms.
