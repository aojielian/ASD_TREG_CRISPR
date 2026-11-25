# Reproduction guide

Use the v1.2.0 code with the processed-data archive identified by the [Zenodo v1.2.0 dataset DOI](https://doi.org/10.5281/zenodo.24532983). The commands below run from the repository root. Keep replay outputs in a new or empty directory outside the repository and data directories. The example creates a new sibling work directory on each run.

## Environment

The pinned workflow uses CPython 3.13.2, Bash and R 4.5.1 with `Rscript` on `PATH`. Python libraries are listed in [requirements.txt](../requirements.txt). Bioconductor dependencies are DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, AnnotationDbi 1.72.0 and org.Mm.eg.db 3.22.0. Optional signature-source reconstruction also uses msigdbr 25.1.1. The workbook builder uses the Python standard library.

Create a Python environment and install the pinned libraries.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

Set `PYTHON` to select an installed interpreter when needed. Install the R packages in an environment compatible with the pinned versions.

## Verify and extract

Use the downloaded archive's actual absolute path for `DATA_ZIP` and retain its published filename. The ZIP and its SHA-256 sidecar must be in the same directory. Keep the original sidecar when available. If the sidecar is unavailable, use the release digest below. The verifier checks this digest against the ZIP before extraction.

```bash
DATA_ZIP="${DATA_ZIP:-/absolute/path/to/ASD_TREG_CRISPR_processed_data_v1.2.0.zip}"
if [ -f "$DATA_ZIP" ] && [ ! -e "${DATA_ZIP}.sha256" ]; then
  (
    set -o noclobber
    printf '%s  %s\n' '06faee76f0a3d04b787712483516be2f8b796e4a43752c0bf9bfefd9014fc923' "${DATA_ZIP##*/}" > "${DATA_ZIP}.sha256"
  )
fi
```

```bash
export PYTHON="${PYTHON:-python3}"
WORK_ROOT="$(mktemp -d "$(cd .. && pwd -P)/ASD_TREG_reproduction.XXXXXX")"
"$PYTHON" scripts/release/verify_release.py "$DATA_ZIP" --extract-to "$WORK_ROOT/verified-data"
DATA="$WORK_ROOT/verified-data/ASD_TREG_CRISPR_processed_data_v1.2.0"
```

The verifier checks the detached archive checksum, ZIP paths and CRC, file checksums, TSV structure, scientific fixtures, and historical compatibility inputs. The workflow restores the necessary historical paths in temporary storage when running the compatibility checks.

## Statistical replay

Run one of these commands from the same shell. Use a separate, new or empty output directory for each run.

```bash
bash workflows/reproduce.sh "$DATA" "$WORK_ROOT/replay"
```

For official-source replay, set `EXTERNAL_INPUT_ROOT` to the local directory that holds the files described under [exact official-source inputs](#exact-official-source-inputs).

```bash
bash workflows/reproduce.sh "$DATA" "$WORK_ROOT/replay-official" "$EXTERNAL_INPUT_ROOT"
```

Archive-only replay covers the processed core, gene-property models, directional loss-of-function and direct FOXP3–TCR analyses, both immune-cell atlases, and mouse conditional-knockout analyses. The immune-cell context branch recomputes 18 tables with archived inputs. Its parity report marks the 24 source-dependent tables as `REQUIRES_EXTERNAL`. Official-source replay checks all 42 tables in that branch and reconstructs the mouse sample, reference-signature and literature evidence inputs. Both modes verify the distributed statistical fixtures.

The workflow writes results under `historical_core/`, `functional_comparisons/`, `immune_cell_context/`, and `conditional_knockout/` within the selected output directory. The [methods and statistics](METHODS_AND_STATISTICS.md) document the analyses and data boundaries.

## Figures and supplementary tables

```bash
bash scripts/release_figures/render.sh "$DATA" "$WORK_ROOT/rendered" --strict-byte-parity --submission-legends
"$PYTHON" scripts/release/build_supplementary_tables.py --source-dir "$DATA/tables/submission" --output "$WORK_ROOT/Supplementary_Tables_S1-S20.xlsx"
```

The figure replay recreates the four main PDFs and their 600 DPI TIFFs and assembles all 23 supplementary pages. Strict byte parity depends on the pinned rendering libraries and font environment, including the Matplotlib DejaVu Sans fonts and PyMuPDF layout fonts. The workbook builder verifies the typed worksheet inputs before assembly. [Figure and table source maps](FIGURE_SOURCE_MAP.md)

## Exact official-source inputs

The extracted archive's `data/manifests/EXTERNAL_REPLAY_INPUTS.tsv` lists approved provider URLs, checksums, file sizes and the paths relative to `EXTERNAL_INPUT_ROOT`. The official-source workflow uses the rows marked `required_for_full_replay`. The `provenance_only` rows record citation or acquisition context. Obtain the originals under their provider terms and verify the checksums before replay. Rolling annotations and updated metadata can differ from the recorded snapshots.

The source notes also describe optional reconstruction with the exact DDBJ E-GEAD-649 counts, HPA input and Zhu conventional CD4 inputs. [Sources and reuse](SOURCES.md) · [Zhu extractor](../scripts/zhu/README.md)
