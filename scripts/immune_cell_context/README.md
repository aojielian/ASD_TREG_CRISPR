# Immune-cell context

See the [release guide](../../docs/REPRODUCIBILITY.md) for the full v1.2.0 replay and environment.

```sh
bash scripts/immune_cell_context/run.sh DATA OUT [EXTERNAL_INPUT_ROOT]
```

With the deposited inputs this branch recomputes the human Monaco and mouse ImmGen atlas results, donor contrasts, expression breadth and sensitivity analyses, plus four historical screen-level comparisons. With exact upstream files it additionally rebuilds the metadata, the direct Umhoefer Treg–Tconv comparison and its sensitivities, and the Chen Perturb-icCITE-seq results. The [input manifest](config/inputs.tsv) identifies the external files and the SHA256s of the archive inputs. The parity table under `08_logs/` lists 42 comparisons. The archive-only mode recomputes 18 results and identifies the other 24 as `REQUIRES_EXTERNAL`. The independent Monaco and ImmGen checks are written under `01_audit/`. The full replay also writes independent Umhoefer and cross-screen numeric checks under `08_logs/`.

- [Monaco methods](config/monaco_methods.md)
- [ImmGen methods](config/immgen_methods.md)
- [Umhoefer methods](config/umhoefer_methods.md) and [table dictionary](config/umhoefer_data_dictionary.md)

The source-specific historical figure layout can be reconstructed with `python3 scripts/immune_cell_context/figures/render.py DATA OUT`. The final submission figures use the [release renderer](../release_figures/README.md).
