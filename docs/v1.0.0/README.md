# ASD_TREG_CRISPR

[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22704659.svg)](https://doi.org/10.5281/zenodo.22704659)

Reproducible analysis and figure-generation code for a computational functional-genomics study of **high-confidence autism spectrum disorder (ASD) risk genes in human FOXP3/regulatory T-cell (Treg) regulatory programs**.

This repository contains the public **code, configuration, documentation, expected results, manifests, and figure-generation workflows**. The complete analysis-ready processed data are distributed separately through Zenodo:

**Zenodo processed-data release:**  
**DOI:** [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)

The GitHub repository and Zenodo dataset are designed to be used together.

---

## Overview

The released workflow supports processed-data reproduction of the principal analyses underlying the study, including:

- enrichment of high-confidence ASD genes in a human genome-wide FOXP3 CRISPR screen;
- independent support in a primary human Treg targeted CRISPR screen;
- direct comparison of FOXP3-regulatory phenotypes with a generic primary T-cell proliferation/TCR-response screen;
- regulator-class-conditioned and stricter joint-conditioned sensitivity analyses;
- cross-screen ASD-gene overlap and influence/recalibration analyses;
- a mouse primary-Treg screen used as a cross-species boundary analysis;
- source-defined FOXP3-associated expression-program analyses in purified ASD Treg transcriptomic data;
- paired whole-blood ASD–sibling program-level analyses;
- reproduction of the final main and supplementary figures from frozen public figure-source tables.

The public release is intended for **transparent reproduction from harmonized processed data**. It does not claim to reconstruct every historical raw-data ingest, source-download, or upstream canonicalization step.

---

## Repository / data organization

The public release is split into two components.

### GitHub — code and documentation

This repository contains:

```text
ASD_TREG_CRISPR/
├── config/                 # frozen public configuration
├── data/                   # documentation / lightweight public data notes
├── docs/                   # analysis map and figure reproducibility documentation
├── environment/            # optional figure-rendering requirements
├── expected_results/       # frozen expected outputs for verification
├── scripts/
│   ├── reproduce_from_processed/
│   └── figures/
├── CHANGELOG.md
├── DATA_AVAILABILITY.md
├── DATA_SOURCES.md
├── LICENSE
├── README.md
├── REPRODUCIBILITY.md
└── RUN_FIRST.md
```

### Zenodo — processed analysis-ready data

Download:

```text
ASD_TREG_CRISPR_processed_data_v1.0.0.zip
```

from:

[https://doi.org/10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)

The Zenodo snapshot contains the row-level processed tables required by the public reproduction scripts, including the analysis universes, conditioning tables, program memberships, patient-context tables, and final figure-source data.

---

## Quick start

### 1. Clone the repository

Using SSH:

```bash
git clone git@github.com:aojielian/ASD_TREG_CRISPR.git
cd ASD_TREG_CRISPR
```

or HTTPS:

```bash
git clone https://github.com/aojielian/ASD_TREG_CRISPR.git
cd ASD_TREG_CRISPR
```

### 2. Download the Zenodo processed-data archive

Download:

```text
ASD_TREG_CRISPR_processed_data_v1.0.0.zip
```

from:

[10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)

Extract it to a local directory.

For example:

```text
/path/to/ASD_TREG_CRISPR_processed_data_v1.0.0/
```

### 3. Run the minimal processed-data reproduction

From the GitHub repository root:

```bash
bash scripts/reproduce_from_processed/run_all_processed.sh \
    /path/to/ASD_TREG_CRISPR_processed_data_v1.0.0
```

The minimal numerical workflow uses **Python 3 standard-library modules only**.

### 4. Confirm successful reproduction

A successful run should end with:

```text
MINIMAL_REPRODUCTION_STATUS=PASS
EXPECTED_RESULTS_FAILURES=0
```

If both conditions are satisfied, the released processed-data workflow has reproduced the expected public numerical results.

---

## What the minimal workflow checks

The processed-data workflow recomputes and verifies, where applicable:

- contingency-table counts;
- one-sided Fisher exact tests;
- rank-based tests;
- direct FOXP3-versus-T-cell-proliferation phenotype-preference statistics;
- regulator-class-conditioned observed/expected enrichment;
- stricter joint-conditioned observed/expected enrichment;
- direct two-test Holm adjustment;
- complete-case coverage;
- cross-screen ASD score-1 gene partitions;
- whole-blood paired sign-flip tests;
- Moreno purified-Treg program-level statistics;
- released mouse boundary counts and Fisher P values.

Detailed implementation and scope notes are provided in:

- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)
- [`RUN_FIRST.md`](RUN_FIRST.md)
- [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)

---

## Reproducing the final figures

The final main and supplementary figures can be regenerated directly from the frozen public figure-source tables included in the Zenodo processed-data snapshot.

Figure rendering is optional and requires the packages listed in:

```text
environment/figure_requirements.txt
```

Current requirements:

```text
numpy
pandas
matplotlib
pypdf
reportlab
```

### Optional: create a clean Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/figure_requirements.txt
```

### Main Figures 1–4

```bash
python scripts/figures/main/render_main_figures.py \
  --input-dir /path/to/ASD_TREG_CRISPR_processed_data_v1.0.0/06_figure_table_source_data/main \
  --output-dir outputs/main_figures
```

### Supplementary Figures S1–S5

```bash
python scripts/figures/supplementary/render_supplement.py \
  --input-dir /path/to/ASD_TREG_CRISPR_processed_data_v1.0.0/06_figure_table_source_data/supplementary \
  --output-dir outputs/supplementary_figures
```

### Assemble the combined supplementary PDF

After rendering Supplementary Figures S1–S5:

```bash
python scripts/figures/supplementary/assemble_supplement.py \
  --output-dir outputs/supplementary_figures
```

The supplementary assembler uses the frozen public legend payload:

```text
scripts/figures/supplementary/figure_legends.json
```

by default.

Additional figure-reproduction documentation is available in:

- [`scripts/figures/README.md`](scripts/figures/README.md)
- [`docs/FIGURE_REPRODUCIBILITY.md`](docs/FIGURE_REPRODUCIBILITY.md)

---

## Data provenance

The study integrates processed information derived from published human and mouse CRISPR screens and transcriptomic resources.

The public release distinguishes:

- primary functional evidence;
- independent human primary-Treg support;
- phenotype-preference analyses;
- conditioned robustness analyses;
- cross-species boundary evidence;
- program-level patient transcriptomic context;
- descriptive/contextual analyses.

The detailed source inventory and public-data provenance are documented in:

[`DATA_SOURCES.md`](DATA_SOURCES.md)

The Zenodo archive contains harmonized, analysis-ready public tables rather than unnecessary redistribution of complete third-party raw datasets.

---

## Reproducibility scope

This repository supports **processed-data minimal reproduction**.

It is designed to allow an independent user to:

1. download the frozen processed-data snapshot;
2. rerun the principal released statistical calculations;
3. verify the expected results;
4. regenerate the final main and supplementary figures.

It does **not** claim to provide a complete from-scratch reconstruction of every original third-party source acquisition or historical upstream data-processing step.

For implementation-specific details, including retained reference quantities and the scope of recomputation, see:

[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)

---

## Expected-results verification

Frozen expected outputs are stored under:

```text
expected_results/
```

The minimal workflow compares newly reproduced values with these released expectations rather than relying only on script completion.

The two final status indicators are:

```text
MINIMAL_REPRODUCTION_STATUS=PASS
EXPECTED_RESULTS_FAILURES=0
```

These provide the primary public verification gate for the processed-data workflow.

---

## Data availability

The analysis-ready processed-data release is permanently archived at Zenodo:

**DOI:** [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)

GitHub intentionally does not duplicate the complete processed-data payload.

See:

[`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)

for the detailed release scope.

---

## Code availability

Public code:

[https://github.com/aojielian/ASD_TREG_CRISPR](https://github.com/aojielian/ASD_TREG_CRISPR)

Processed data:

[https://doi.org/10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)

---

## Version

Current public release:

```text
v1.0.0
```

See [`CHANGELOG.md`](CHANGELOG.md) for release notes.

---

## Citation

If you use the processed dataset, please cite the Zenodo record:

> ASD_TREG_CRISPR processed data, version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.22704659

A manuscript-specific citation can be added after formal publication of the associated article.

---

## License

Original code in this repository is distributed under the terms provided in [`LICENSE`](LICENSE).

Third-party source data remain subject to their original source terms and citations. The processed-data archive should be used together with the provenance information in [`DATA_SOURCES.md`](DATA_SOURCES.md) and the metadata attached to the Zenodo record.

---

## Key documentation

For first-time users:

- [`RUN_FIRST.md`](RUN_FIRST.md) — shortest reproduction route
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) — numerical reproduction scope and implementation details
- [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md) — what is distributed and where
- [`DATA_SOURCES.md`](DATA_SOURCES.md) — source datasets and provenance
- [`scripts/figures/README.md`](scripts/figures/README.md) — figure-rendering commands
- [`docs/FIGURE_REPRODUCIBILITY.md`](docs/FIGURE_REPRODUCIBILITY.md) — final figure-source mapping and reproduction notes

---

## Contact / issues

For reproducibility problems, please open a GitHub Issue and include:

- operating system;
- Python version;
- command executed;
- relevant error message;
- whether the numerical reproduction or figure-rendering workflow was being run.

Do not include private participant information or restricted source data in public issues.
