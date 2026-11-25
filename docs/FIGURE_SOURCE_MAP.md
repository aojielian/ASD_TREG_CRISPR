# Figure and table source maps

The figure map lists the current source paths relative to the processed-data root, SHA256 digests and rendering commands for each main panel and Supplementary Figures S1–S23. It identifies reused vector components and the source data for re-rendered panels. [Figure source map](FIGURE_SOURCE_MAP.tsv)

The table map links each of the 76 scientific worksheets in Supplementary Tables S1–S20 to its canonical TSV, original source context and column dictionary. Typed cell values and spreadsheet metadata needed to reconstruct the workbook are under `tables/submission/`. [Supplementary table source map](SUPPLEMENTARY_TABLE_SOURCE_MAP.tsv)

Publication figures and legends preserve original panel numbering and scientific notation. The workbook index retains reference labels from the released tables where needed for byte-for-byte workbook reproduction. The table source map and accompanying dictionary resolve those labels to the current paths. [Data dictionary](DATA_DICTIONARY.md) · [Build commands](REPRODUCIBILITY.md#figures-and-supplementary-tables)
