# Figure source map

Source paths are relative to the extracted v1.1.0 data directory. The renderers use released data rows and deterministic display transformations. Main figure tables retain the published compatible paths. The supplementary source directory points back to the original table and SHA-256 where it carries a frozen copy.

| Figure and panel | Data input | Recorded source and transformation |
| --- | --- | --- |
| Figure 1a | `06_figure_table_source_data/main/Figure1_source_data.tsv` | All 18,966 `screen_gene` rows, using native Chen LFC and −log10 of the smaller native one-sided P. Marker groups use score-1, hit and standardized direction flags. |
| Figure 1b,c | Same table | The enrichment rows contain 14/242 and 243/18,724. The 14 points in c are the score-1 hit subset of the scatter, sorted by native LFC. |
| Figure 2a–c | `06_figure_table_source_data/main/Figure2_source_data.tsv` | Primary Treg Fisher margins, 33 `common_asd_gene` rows, and the five genes with a hit in either experiment. Panels b and c respectively summarize and subset the common measured rows. S1 retains the full matrix. |
| Figure 3a,b | `06_figure_table_source_data/main/Figure3_source_data.tsv` | Measured Chen–Shifrut intersection and counts of the four joint hit categories from the 18,615 shared rows, with 242 score-1 genes. |
| Figure 3c | Same table | Observed difference from FOXP3-only and TCR-only hit counts. Exact discrete null over sampling without replacement across the plus, minus and neutral categories. The visible score range is a subset of the full exact support. |
| Figure 3d | `data/r7/figure_sources/Figure_3d_source.tsv` | Eight same-scenario M0–M3 OR/CI records, and the two M2 Holm P values. The forest uses a log odds-ratio scale and preserves the source intervals. |
| Figure 4a,b | `06_figure_table_source_data/main/Figure4_source_data.tsv` | Frozen M1–M10 membership coverage and eligible-program feature-level distribution comparison, with eligibility and BH10 documented. |
| Figure 4c | Same table | Step ECDFs for programs M1, M6, M8 and M9 versus other features. Uses all 5,851 finite Moreno effects, including repeated or unmapped feature symbols, with a symmetric-log display. |
| S1 | `figure_sources/supp/S1_source.tsv` | 33 jointly measured score-1 genes with original Chen and Umhoefer Treg / Tconv signed effects and hit status. |
| S2 | `figure_sources/supp/S2_source.tsv` | Measured, recalibrated and both-hit-removal enrichments, class-conditioned O/E, and full leave-one-gene-out P series. |
| S3 | `figure_sources/supp/S3_source.tsv`, `S3_conditioning.tsv`, `S3_regression.tsv` and `S3_model_sensitivities.tsv` in the same directory | Measured regulator strata, exact conditional analyses by reference, coverage gates and fixed complete-case models. Historical and exploratory rows are distinct. |
| S4 | `figure_sources/supp/S4_source.tsv` | Twenty anonymous paired differences per program, exact two-sided sign-flip P and BH10 q, paired-t display intervals and descriptive Yan list overlap. |
| S5 | `figure_sources/supp/S5_source.tsv` | Moreno feature-level program and DEG overlap, including the 0/14 discovery-gene outcome. |
| S6 | `figure_sources/supp/S6_membership.tsv` and `figure_sources/supp/S6_candidates.tsv` | G2P mechanism category counts, LoF-compatible sensitivity, 17 standardized candidate directions and assay availability. The distinct conventional CD4 Zhu context is described in the S6 legend and S12 data, with its numerical source under `data/r8/zhu_derived/`. |

`figures/render_main.py` renders Figures 1–4. `figures/render_supp.py` renders S1–S6 and combines those pages. See the [reproduction guide](REPRODUCIBILITY.md) for commands. The `06_figure_table_source_data/FIGURE_SOURCE_DATA_MAP.tsv` index describes the earlier display renderers. The v1.1.0 sources and guide here apply to the final renderers.
