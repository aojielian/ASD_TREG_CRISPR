# ASD_TREG_CRISPR

[![Zenodo data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24487407.svg)](https://doi.org/10.5281/zenodo.24487407)

This repository reproduces a functional-genomics analysis of high-confidence autism spectrum disorder (ASD) risk genes in published human FOXP3 perturbation screens. Version 1.1.0 preserves the v1.0.0 processed core and adds expression-detectable backgrounds, gene-property sensitivities, allelic-mechanism annotations, published experimental evidence, and source-based renderers for Figures 1–4 and Supplementary Figures S1–S6.

| Resource | Location |
| --- | --- |
| Code, analysis thresholds and numerical fixtures | This repository |
| Companion processed dataset | `ASD_TREG_CRISPR_processed_data_v1.1.0.zip` · [10.5281/zenodo.24487407](https://doi.org/10.5281/zenodo.24487407) |
| Collection DOI | [10.5281/zenodo.22704658](https://doi.org/10.5281/zenodo.22704658) |
| Published v1.0.0 dataset | [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659) |

## Reproduce the results and figures

Download `ASD_TREG_CRISPR_processed_data_v1.1.0.zip` from the [v1.1.0 Zenodo record](https://doi.org/10.5281/zenodo.24487407). From the repository root, with CPython 3.13 and base R available:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -r environment/r7_figure_requirements.txt
bash scripts/smoke_test.sh ../ASD_TREG_CRISPR_processed_data_v1.1.0.zip release-check
```

The workflow uses the companion ZIP or its extracted data directory. It verifies the data inventory, recomputes the original screen, specificity and program analyses, refits the M0–M3 grid, checks the fixed mechanism/evidence tables and re-renders all ten figures. Base R independently checks primary profiles and the historical conditional odds ratio. The optional raw-source reconstructions require files from the original providers. See the [reproduction guide](docs/REPRODUCIBILITY.md), [data dictionary](docs/DATA_DICTIONARY.md), [source guide](docs/SOURCES.md) and [figure source map](docs/FIGURE_SOURCE_MAP.md).

## Analysis scope

| Comparison | Chen human iTreg | Umhoefer human primary Treg |
| --- | ---: | ---: |
| Original measured score-1 genes, hits/tested | 14/242 | 5/33 |
| Expression-detectable score-1 genes at CPM ≥1, hits/tested | 14/174 | 5/18 |
| Same-complete-case M2 OR [95% penalized profile CI] | 1.332 [0.706, 2.345] | 3.838 [0.998, 13.097] |
| M2 two-screen Holm P | 0.360 | 0.101 |

M2 adjusts for regulator class, expression and LOEUF, with targeted-library membership also included for Chen. Neither M2 test reached the adjusted 5% threshold. The complete results retain the exploratory primary-Treg M3 estimate, the expression and coverage sensitivities, and condition-specific failures of model estimation. The [methods](docs/METHODS_AND_STATISTICS.md) specify the denominators and interpretation for each comparison.

The frozen score-1 list contains 245 genes. The 104-gene LoF-compatible sensitivity includes 100 inferred and four evidence-tagged classifications under the disease-scope rules. Published evidence records retain the assay, perturbation, cellular context and endpoint at experiment level. Zhu et al. (2026) add candidate-level GATA3 protein evidence in polarized conventional CD4 T cells. The FOXP3 results and adjustment model definitions remain unchanged. [Zhu et al., Cell](https://doi.org/10.1016/j.cell.2026.08.002)

The [MIT license](LICENSE) covers the repository code. Source-specific terms and attribution apply to upstream data. Cite the applicable data version alongside the code. See the [source guide](docs/SOURCES.md) and [changelog](CHANGELOG.md).
