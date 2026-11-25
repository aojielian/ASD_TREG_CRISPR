# Data dictionary

Read `DATA_DICTIONARY.md` in the extracted archive for detailed table fields and units. Data paths in this repository's maps are relative to that archive root.

| Canonical location | Contents |
| --- | --- |
| `data/processed/` | Harmonized core gene, stratum and program tables, including anonymous whole-blood pair differences. |
| `data/crispr_screens/` | Original screen-level adjudications and comparison source projections. |
| `data/gene_property_adjustment/` | Detection scenarios, gene covariates, Firth profile models, sensitivity and matching tables, and scoped allelic mechanism. |
| `data/functional_comparisons/` | Directional loss of function, FOXP3–TCR, and matched Treg–Tconv estimands, with their distinct universes, strata, null distributions and checks. |
| `data/immune_cell_expression/` | Human Monaco TPM and mouse ImmGen TMM CPM results, processed matrices, sample and gene mappings, breadth, sensitivity and matching. |
| `data/conditional_knockout/` | Published GEO gene counts, mouse library and gene mappings, primary and sensitivity differential expression, prespecified transcripts and CAMERA results. |
| `data/experimental_evidence/` | Published assays, Perturb-icCITE-seq, conventional CD4 Zhu tables, bounded literature search and candidate-level integration. |
| `data/manifests/` | Source and reuse ledgers, canonical file checksums, official external inputs and temporary compatibility bindings. |
| `figure_sources/` | Per-panel numeric inputs, supplementary sources and fixed vector components. |
| `tables/submission/` | Public per-worksheet TSVs, column dictionary, index and portable typed cell data. |

Missing or unmeasured outcomes retain their recorded `NA`, blank or status field. DNA-level labels, screen-defined directions, FOXP3 protein abundance and mouse Foxp3 RNA abundance refer to their own data types and biological contexts. [Methods](METHODS_AND_STATISTICS.md) · [Figure and table source maps](FIGURE_SOURCE_MAP.md)
