# ASD_TREG_CRISPR

## Overview

This repository accompanies the ASD_TREG_CRISPR study, which examines how autism-associated genes relate to FOXP3 regulation and regulatory T-cell programs. The analyses integrate human perturbation screens, immune-cell expression references, and published mouse conditional-knockout data.

## Repository contents

- `workflows/` and `scripts/` contain the statistical analyses, validation tools, and figure and table builders.
- `data/manifests/` and `expected_results/` identify source inputs and record checksums for the data and reference results.
- `docs/` describes the methods, data fields, source attribution, and reproduction procedures.

## Data availability

The companion dataset for this release is identified by [Zenodo dataset DOI](https://doi.org/10.5281/zenodo.24532983). The archive contains the analytical inputs, statistical results, and source data for Figures 1–4, Supplementary Figures S1–S23, and Supplementary Tables S1–S20.

Use `ASD_TREG_CRISPR_processed_data_v1.2.0.zip` and keep its matching SHA-256 sidecar in the same directory. The [archive verification instructions](docs/REPRODUCIBILITY.md#verify-and-extract) give the exact checksum and explain how to create the sidecar if needed. Use the [source instructions](docs/REPRODUCIBILITY.md#exact-official-source-inputs) to obtain required third-party files separately.

## Software requirements

Use the documented Python, Bash, R and Bioconductor environment. Python package versions are pinned in [requirements.txt](requirements.txt), and R package versions are listed in the [reproduction guide](docs/REPRODUCIBILITY.md#environment).

## Reproducing the analyses

Run these commands from the repository root. Set `DATA_ZIP` to the absolute path of the downloaded file. Each run creates a new sibling directory for its outputs, preserving earlier results.

```bash
DATA_ZIP="${DATA_ZIP:-/absolute/path/to/ASD_TREG_CRISPR_processed_data_v1.2.0.zip}"
export PYTHON="${PYTHON:-python3}"
WORK_ROOT="$(mktemp -d "$(cd .. && pwd -P)/ASD_TREG_reproduction.XXXXXX")"
"$PYTHON" scripts/release/verify_release.py "$DATA_ZIP" --extract-to "$WORK_ROOT/verified-data"
DATA="$WORK_ROOT/verified-data/ASD_TREG_CRISPR_processed_data_v1.2.0"
bash workflows/reproduce.sh "$DATA" "$WORK_ROOT/replay"
```

The workflow above uses the archived inputs. See the [reproduction guide](docs/REPRODUCIBILITY.md#statistical-replay) for the scope of each replay mode.

## Reproducing figures and supplementary tables

In the same shell session, run:

```bash
bash scripts/release_figures/render.sh "$DATA" "$WORK_ROOT/rendered" --strict-byte-parity --submission-legends
"$PYTHON" scripts/release/build_supplementary_tables.py --source-dir "$DATA/tables/submission" --output "$WORK_ROOT/Supplementary_Tables_S1-S20.xlsx"
```

Strict figure parity depends on the documented software and font environment. The [source maps](docs/FIGURE_SOURCE_MAP.md) identify the data used by each panel and worksheet.

## Citation

Cite the software and its companion dataset using [CITATION.cff](CITATION.cff).

## License

The code is distributed under the [MIT License](LICENSE). The contributors' processed-data compilation is available under CC BY 4.0, with source-specific terms and attribution documented in the archive. [Sources and reuse](docs/SOURCES.md)
