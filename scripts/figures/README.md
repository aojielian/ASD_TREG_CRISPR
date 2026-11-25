# Compatibility figure renderers

These scripts reproduce the earlier figure presentation for compatibility checks. Use the [reproduction guide](../../docs/REPRODUCIBILITY.md#figures-and-supplementary-tables) for the current manuscript figure set.

Use the `DATA`, `PYTHON`, and `WORK_ROOT` values initialized by the [archive verification workflow](../../docs/REPRODUCIBILITY.md#verify-and-extract). The following commands rebuild the required input paths from the current release and place the outputs in the external work directory. The [optional figure dependencies](../../environment/figure_requirements.txt) include the libraries for this compatibility renderer.

```bash
"$PYTHON" scripts/release/verify_release.py "$DATA" --legacy-view "$WORK_ROOT/figure-compatibility-inputs"

"$PYTHON" scripts/figures/main/render_main_figures.py \
  --input-dir "$WORK_ROOT/figure-compatibility-inputs/06_figure_table_source_data/main" \
  --output-dir "$WORK_ROOT/compatibility-main-figures"

"$PYTHON" scripts/figures/supplementary/render_supplement.py \
  --input-dir "$WORK_ROOT/figure-compatibility-inputs/06_figure_table_source_data/supplementary" \
  --output-dir "$WORK_ROOT/compatibility-supplementary-figures"

"$PYTHON" scripts/figures/supplementary/assemble_supplement.py \
  --output-dir "$WORK_ROOT/compatibility-supplementary-figures"
```

The supplementary assembler uses `figure_legends.json` in its script directory. The `--legends` option accepts an explicit path.
