# Data dictionary

All paths below are relative to the extracted processed-data root. Numeric columns use decimal points. Empty numeric cells indicate an unavailable value or a non-estimable test, as recorded by the applicable status field. A case gene is in the frozen SFARI 2026 Q2 score-1 group. A control gene in enrichment analyses is an eligible screened gene outside that group. These labels are gene-level, with no patient risk prediction intended.

| Path | Unit, key and interpretation |
| --- | --- |
| `data/processed/` | Legacy v1.0.0 analysis-ready gene, stratum and program results. Definitions in `data/manifests/v1.0.0/DATA_DICTIONARY.tsv` and `COLUMN_DICTIONARY.tsv`. The whole-blood pair table retains anonymous ordinal, program, paired difference and sign. |
| `06_figure_table_source_data/` | Legacy figure-source rows and source map. The S4 table retains aggregate data and anonymous paired differences. |
| `data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz` | Unique `(screen, approved_symbol)` pairs. Keeps measured genes plus explicitly absent frozen ASD genes. Includes source symbols/LFC/FDR orientation, harmonized membership, detectability by threshold, mean log expression, freeze-mapped gene constraint, span and guide counts. |
| `data/r7/inputs/R6_GENE_COVARIATES.tsv` | One gene per core screen. Historical model covariates, including the original standardized expression and LOEUF terms used in the bridge. |
| `data/r7/inputs/R1_CONDITIONING_STRATA.tsv` | Unique scenario–screen–stratum records with total `N`, author hit count `H`, ASD count `A` and ASD hits `K`. Retains non-informative strata. |
| `data/r7/inputs/R1_MATCHING_QC.tsv` | Historical class/strict balance, comparator support and matching notes. `record_type` identifies row schema. |
| `data/r7/results/R1_FROZEN_SENSITIVITY_GRID.tsv` | Frozen Fisher, exact conditional and historical sensitivity results. `effect_type` distinguishes OR and O/E; `null_lo, null_hi` are conditional-null hit-count bounds, separate from an OR interval. |
| `data/r7/results/LOEUF_ADJUSTED_MODELS.tsv` | Unique `screen,reference,threshold,units,model` rows. M0–M3 use the same complete genes within each scenario. 33 rows include the Chen no-library comparison. Blank estimates for unstable models remain absent. |
| `data/r7/results/LOEUF_MODEL_COEFFICIENTS.tsv` | Model term coefficients, expected-information standard error and VIF where defined. Coefficients use log odds. |
| `data/r7/results/R6_TO_R2_MODEL_BRIDGE.tsv` | One row for each screen and historical or bridge model condition. New penalized LR P and legacy Wald diagnostics are named separately. |
| `data/r7/results/MEASURED_TO_DETECTABLE_FLOW.tsv` | Mutually described flow stages or exclusions within screen and threshold, counted in genes, ASD genes and hits. `ASSAYED`, `REFERENCE_MISSING`, `BELOW_THRESHOLD`, `DETECTABLE`, strict complete and informative are distinct. |
| `data/r7/results/PROPERTY_BALANCE_BEFORE_AFTER.tsv` | Screen, universe/phase and covariate SMDs. `weight_rule` records unit-gene or strata weights. |
| `data/r7/results/REGULATOR_CLASS_BREAKDOWN.tsv` | Counts and Fisher statistics by `screen,reference,threshold,regulator_class,background`. |
| `data/r7/results/MATCHING_FEASIBILITY.tsv` and `MATCH_CASE_SUPPORT.tsv` | Fixed-caliper matching feasibility, paired gene identities, indicator margins and predeclared coverage limits. Matched P is absent when coverage or balance fails. |
| `data/r7/results/DIRECTION_ORIGINAL_TO_MANUSCRIPT.tsv` | Native positive and negative LFC, measured/hit states and standardized KO direction for 17 candidates across three screens. Historical filename retained for traceability. The corresponding full manuscript is outside this archive. |
| `data/r7/mechanism/` | Gene-level annotation master (245 genes), gene–disease records (294 rows including explicit missing records), 17 candidate crosswalk, 51 screen context rows and scoped LoF count and change ledgers. |
| `data/r7/function/` | 14 grouped published experiment contexts plus 208 linked source-evidence or coverage rows; PMIDs, DOIs, figure/table citations, endpoints and distinctions between expression, phenotype and actual suppression. |
| `data/r7/figure_sources/` | Final Figure 3d inputs, S6 candidate/readout and membership inputs, legends, and version-linked input indexes. |
| `data/r7/provenance/` | Versioned source registration and file/column schema for projected inputs and released results. |

## Core fields

- **Identity and provenance.** `approved_symbol` is the freeze-approved HGNC symbol. `raw_primary_id` is the symbol or identifier in the source screen. `source_row` is the original workbook position, never a sample ID. `measured` and `author_hit` have separate meanings. Unmeasured rows carry no result.
- **Expression and counts.** `primary_mean_log_expression` is the mean of log2(CPM + 1) per reference library. Mean CPM is stored separately. `*_group*_n_ge_*` counts libraries meeting that cutoff. Each Chen iTreg NTC guide requires two of three libraries. Each hTreg/hTconv reference requires three of four libraries. `hpa_mean_log_expression` is log2(nTPM + 1) from the appropriate HPA aggregate. Different expression units are never merged.
- **Gene attributes.** `LOEUF_v2_1_1` is the gnomAD v2.1.1 upper observed/expected pLoF bound. `GENE_SPAN_BP` is genomic span in base pairs. `num` is guides per gene from the author table. TF and EpiFactors labels represent source catalogue membership. `library_member` refers to the targeted library.
- **Direction.** The author `pos|lfc` and `neg|lfc` fields keep the source contrast. The interpreted `KO_REDUCES_FOXP3` and `KO_INCREASES_FOXP3` labels refer to the screened protein endpoint and context. They do not predict the effect of a patient allele.
- **Models.** M0 uses only the ASD label. M1 adds catalogue flags, plus library membership in Chen. M2 adds continuous expression and LOEUF. M3 adds genomic span and guide count. Reported ORs are exp(ASD log-odds coefficient); P and CI are based on the Firth/Jeffreys penalized profile unless separately identified. The primary M2 Holm adjustment is over the two CPM=1 core tests.
- **Mechanism and evidence.** `molecular_mechanism_support` distinguishes the source database's inferred and evidence-tagged status. `lof_eligible` is the final scoped analysis flag; its evidence includes non-ASD neurodevelopmental disorders. Explicit placeholders preserve missing coverage without interpreting it as LoF or GoF. Experiment groups link to individual published endpoint or coverage entries.

## Final figure and supplementary-table sources

| Path | Unit, key and interpretation |
| --- | --- |
| `06_figure_table_source_data/main/Figure1_source_data.tsv` | `record_type` partitions gene, enrichment and direction rows. The 18,966 measured genes use native screen LFC and minimum tail P for plotted strength. |
| `06_figure_table_source_data/main/Figure2_source_data.tsv` | Enrichment records and 33 uniquely identified common genes. Hit flags distinguish the two cell contexts. |
| `06_figure_table_source_data/main/Figure3_source_data.tsv` | Gene rows for the 18,615 shared measured universe and source summary records. `D` is +1 for FOXP3-only, −1 for proliferation-only, 0 otherwise. |
| `06_figure_table_source_data/main/Figure4_source_data.tsv` | Program summaries and 5,851 finite source feature records. Multiple features can map to one symbol and unmapped features remain in the denominator. Program eligibility counts matched symbols. |
| `figure_sources/main/` | Source index for the final Figures 1–4 and numeric support tables where applicable. Existing compatible source paths remain in use. |
| `figure_sources/supp/SOURCE_MAP.tsv` | One frozen copy per source table with the relative original companion path and SHA-256. S1 carries the 33-gene screen matrix, S2 carries sensitivity/omission rows, S3 covers conditional, model and reference results, S4 carries 20 anonymous pair differences per program and list overlap, S5 carries feature-level program/DEG results, and S6 covers mechanism and candidate-assay context. |
| `tables/TABLE_INDEX.tsv` | One row per substantive supplementary worksheet projection. `Table` is the stable S1–S12 identifier, `Worksheet` is the current sub-sheet name and `Public TSV` is relative to `tables/`. |
| `tables/COLUMN_SCHEMA.tsv` | `file, ordinal` identify each released worksheet column. Includes emitted headings and role. The three S9 TSVs retain their multi-section data and original worksheet row numbers. |
| `tables/` | Original analysis and evidence rows at spreadsheet precision. S12 retains 245 mechanism master rows and 17 candidate rows, with focused Zhu candidate coverage in S12f. Source-specific endpoint values retain their own P and correction labels. Local-only source pointer columns are omitted from the public equivalents. |

## Zhu conventional CD4 candidate evidence

| Path within `data/r8/zhu_derived/` | Unit, key and interpretation |
| --- | --- |
| `ZHU_CANDIDATE_17_SUMMARY.tsv` | One record per frozen candidate symbol. Guide records carry source quality flags. DE coverage refers to a reported perturbation-condition summary, with lack of a row represented as missing assessment. |
| `ZHU_CANDIDATE_17_BY_CONDITION.tsv` | One record per `gene, culture_condition` across 17 genes and Rest, Stim8hr and Stim48hr. Of 51 entries, 46 are observed DE records. Counts distinguish original cells, guides and DE genes. Cells summarize coverage and are not inferential replicates. |
| `ZHU_ARRAYED_VALIDATION_CANDIDATES.tsv` | One record per available selected candidate/condition/modality/endpoint. `source_effect` uses `effect_measure`. P and adjusted P are original source statistics. Flow matching counts use distinct donor labels within assay and batch. Bulk RNA donor counts describe a separate assay. |
| `ZHU_SOURCE_FILES.tsv` | One record per original source filename with bytes, SHA-256, direct source link, role and redistribution disposition. |
| `ZHU_BIBLIOGRAPHY.tsv` | One publication with DOI, PMID, article and processed-dataset locators, terms and cell context. |

The detailed derivative dictionary records field-level definitions and exact source columns. A missing guide or DE result, a measured nonhit, and a reported effect whose adjusted P exceeds a threshold are distinct. The ZMYND8 result pertains to conventional T-cell Th1 GATA3 median fluorescence intensity. It remains separate from all primary-Treg FOXP3 outcomes.

## Data encoding

TSVs are UTF-8, tab separated, with a header row. Fields containing delimiters, quotation marks or line breaks are CSV-quoted under the tab-delimited dialect. Use the relevant status field when interpreting a blank statistic. Raw one-sided enrichment P, penalized profile P, Holm adjustments, BH q and source-reported adjusted P are distinct columns or documented by the row's analysis label. The current manifest covers released bytes, rather than expected tolerances for recomputed numeric output.

See the [source guide](SOURCES.md) and the companion data archive’s Zhu per-column dictionary.
