# Figure reproducibility scripts

These scripts reproduce the final main figures and supplementary figures from the public figure source tables in the companion analysis-ready snapshot.

The core numerical reproduction remains separate and uses only Python standard-library modules. Figure rendering is optional and uses the packages listed in `environment/figure_requirements.txt`.

Example commands from the repository root:

```sh
python scripts/figures/main/render_main_figures.py \
  --input-dir path/to/extracted/ASD_TREG_CRISPR_processed_data_v1.0.0/06_figure_table_source_data/main \
  --output-dir outputs/main_figures

python scripts/figures/supplementary/render_supplement.py \
  --input-dir path/to/extracted/ASD_TREG_CRISPR_processed_data_v1.0.0/06_figure_table_source_data/supplementary \
  --output-dir outputs/supplementary_figures

python scripts/figures/supplementary/assemble_supplement.py \
  --output-dir outputs/supplementary_figures
```

The renderers read relative/argument-provided paths and do not require raw datasets.

The supplementary assembler reads the exact frozen public legend payload from `scripts/figures/supplementary/figure_legends.json` by default. An explicit `--legends` path may be supplied when reproducing in another directory.

Processed figure-source tables are distributed in the Zenodo v1.0.0 dataset: [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659).
