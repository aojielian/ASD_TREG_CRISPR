# Figure reproducibility

Processed figure-source data are deposited in Zenodo: [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659).

The final figure source tables are in the companion analysis-ready snapshot under `06_figure_table_source_data/`. The renderer inputs are the same public tables used to generate the final figure artifacts.

| Figure | Public source table | Renderer |
|---|---|---|
| Figure 1 | `06_figure_table_source_data/main/Figure1_source_data.tsv` | `scripts/figures/main/render_main_figures.py` |
| Figure 2 | `06_figure_table_source_data/main/Figure2_source_data.tsv` | `scripts/figures/main/render_main_figures.py` |
| Figure 3 | `06_figure_table_source_data/main/Figure3_source_data.tsv` | `scripts/figures/main/render_main_figures.py` |
| Figure 4 | `06_figure_table_source_data/main/Figure4_source_data.tsv` | `scripts/figures/main/render_main_figures.py` |
| Supplementary Figures S1–S5 | `06_figure_table_source_data/supplementary/Supplementary_Figure_S1_source_data.tsv` through `Supplementary_Figure_S5_source_data.tsv` | `scripts/figures/supplementary/render_supplement.py` |

Display transformations are limited to sorting, percentages, `-log10(P)` labels/coordinates where specified, empirical cumulative distribution steps, and the symmetric-log axis used for signed-effect display. Figure 3 categorical squares encode real membership categories and are not continuous coordinates.

Figure 4 M1/M6/M8/M9 curves are empirical step ECDFs from the real Moreno feature rows, with no density fitting or smoothing. The M1–M10 short descriptors are source-derived annotations recorded in the Figure 4 table.

Supplementary Figure S4 uses the real 20 pairwise ASD-minus-sibling differences per program. The displayed interval is the conventional paired-t 95% interval, `mean +/- t_(0.975,19) * sample_sd / sqrt(20)`. It is a display interval only; inferential P values are the exact two-sided sign-flip tests, with BH10 q values retained separately.

The core numerical reproduction does not require these optional plotting dependencies.


The renderers read the same public source-table filenames listed above; no filename alias or link is required.
