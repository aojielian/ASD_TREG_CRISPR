# ASD_TREG_CRISPR

[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22704659.svg)](https://doi.org/10.5281/zenodo.22704659)

This public release contains code, configuration, documentation, manifests, expected results, and a standard-library processed-data reproduction workflow. Complete row-level analysis-ready tables are distributed in the companion Zenodo dataset ([DOI: 10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)); GitHub does not duplicate them.

## Minimal reproduction

Download the processed-data snapshot from Zenodo ([DOI: 10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659)), extract it, place the GitHub repository and extracted snapshot in accessible locations, and run:

    bash scripts/reproduce_from_processed/run_all_processed.sh path/to/ASD_TREG_CRISPR_processed_data_v1.0.0

The workflow uses Python 3 standard library modules only. It recomputes the released counts, Fisher and rank-sum tests, direct phenotype preference, class and strict joint conditioning, Loo boundary counts/P, whole-blood paired tests, and Moreno program tests, then compares all expected rows. Class and joint O/E, exact P, direct two-test Holm P, complete-case coverage, and the common ASD score-1 gene partition are checked explicitly.

The public workflow supports processed-data minimal reproduction. It does not claim a complete rebuild from raw files or every original third-party source-ingest step. No manuscript or patient identifiers are included. See `DATA_AVAILABILITY.md`, `DATA_SOURCES.md`, and the data manifests for scope and provenance.
