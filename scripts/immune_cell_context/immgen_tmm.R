#!/usr/bin/env Rscript
# Normalize the complete vM25 count matrix before selecting atlas comparisons.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
    stop("Usage: Rscript immgen_tmm.R DOWNLOAD_DIR OUT_DIR")
}
download_dir <- arguments[[1L]]
out_dir <- arguments[[2L]]
audit_dir <- file.path(out_dir, "01_audit")
result_dir <- file.path(out_dir, "03_results")
count_file <- file.path(download_dir, "GSE109125",
                        "GSE109125_Gene_count_table_GENCODE_vM25.csv.gz")

read_tab <- function(path) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
               na.strings = "NA", quote = "", comment.char = "")
}

write_tab <- function(data, path, compressed = FALSE) {
    connection <- if (compressed) gzfile(path, "wt", compression = 6) else file(path, "wt")
    on.exit(close(connection))
    write.table(data, file = connection, sep = "\t", row.names = FALSE,
                quote = FALSE, na = "NA")
}

soft <- read_tab(file.path(audit_dir, "R10_IMMGEN_SAMPLE_CROSSWALK.tsv"))
group_data <- read_tab(file.path(audit_dir, "R10_IMMGEN_GROUP_CROSSWALK.tsv"))
panel <- read_tab(file.path(audit_dir, "R10_IMMGEN_GENE_CROSSWALK.tsv"))
raw <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(names(raw)[[1L]] == "GeneSymbol", !anyDuplicated(raw[[1L]]),
          !anyDuplicated(names(raw)))

counts <- as.matrix(raw[-1L])
rownames(counts) <- raw[[1L]]
rm(raw)
stopifnot(is.numeric(counts), all(is.finite(counts)), all(counts >= 0),
          all(counts == floor(counts)), all(nzchar(rownames(counts))))
match_gsm <- match(colnames(counts), soft$matrix_library)
stopifnot(!anyNA(match_gsm), !anyDuplicated(match_gsm))
library_meta <- soft[match_gsm, ]
depth <- colSums(counts)
stopifnot(all(depth > 0))

# edgeR uses all input libraries. The all-zero rows have no effect on TMM.
y <- edgeR::DGEList(counts = counts)
y <- edgeR::normLibSizes(y, method = "TMM")
factors <- y$samples$norm.factors
stopifnot(all(is.finite(factors)), all(factors > 0),
          abs(sum(log(factors))) < 1e-12, all(y$samples$lib.size == depth))
tmm_cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
chosen_rows <- c(1L, nrow(counts) %/% 2L, nrow(counts))
expected <- sweep(counts[chosen_rows, , drop = FALSE], 2L, depth * factors / 1e6, "/")
stopifnot(isTRUE(all.equal(tmm_cpm[chosen_rows, , drop = FALSE], expected)))

primary_groups <- group_data[group_data$include_primary == "yes", ]
primary <- library_meta$include_primary == "yes"
median_cpm <- vapply(primary_groups$cell_type, function(group) {
    idx <- which(library_meta$cell_type == group & primary)
    stopifnot(length(idx) == primary_groups$n_libraries[primary_groups$cell_type == group])
    if (length(idx) == 1L) {
        return(tmm_cpm[, idx])
    }
    matrixStats::rowMedians(tmm_cpm, cols = idx)
}, numeric(nrow(counts)))
stopifnot(ncol(median_cpm) == nrow(primary_groups), all(is.finite(median_cpm)))
colnames(median_cpm) <- primary_groups$cell_type
median_export <- data.frame(mouse_gene = rownames(counts), median_cpm, check.names = FALSE)
write_tab(median_export, file.path(result_dir, "R10_IMMGEN_CELLTYPE_CPM.tsv.gz"), compressed = TRUE)

measured <- panel[panel$measured == "yes", ]
mouse_genes <- measured$mouse_gene
target_counts <- counts[mouse_genes, , drop = FALSE]
target_normal <- tmm_cpm[mouse_genes, , drop = FALSE]
ordinary_cpm <- sweep(target_counts, 2L, depth / 1e6, "/")
target_values <- data.frame(
    human_gene = rep(measured$human_gene, times = ncol(counts)),
    mouse_gene = rep(mouse_genes, times = ncol(counts)),
    panel_role = rep(measured$panel_role, times = ncol(counts)),
    matrix_library = rep(colnames(counts), each = nrow(measured)),
    gsm = rep(library_meta$gsm, each = nrow(measured)),
    cell_type = rep(library_meta$cell_type, each = nrow(measured)),
    include_primary = rep(library_meta$include_primary, each = nrow(measured)),
    matrix_match = rep(library_meta$matrix_match, each = nrow(measured)),
    raw_gene_count = as.vector(target_counts),
    tmm_cpm = as.vector(target_normal),
    ordinary_cpm = as.vector(ordinary_cpm)
)
write_tab(target_values, file.path(result_dir, "R10_IMMGEN_PANEL_LIBRARY_CPM.tsv.gz"), compressed = TRUE)

qc <- data.frame(
    dataset = "GSE109125",
    matrix_library = colnames(counts),
    gsm = library_meta$gsm,
    cell_type = library_meta$cell_type,
    include_primary = library_meta$include_primary,
    raw_library_size = as.numeric(depth),
    tmm_factor = factors,
    tmm_effective_library_size = as.numeric(depth * factors),
    n_nonzero_genes = colSums(counts > 0),
    stringsAsFactors = FALSE
)
write_tab(qc, file.path(audit_dir, "R10_IMMGEN_LIBRARY_QC.tsv"))
schema <- data.frame(
    dataset = "GSE109125",
    source = "GSE109125_Gene_count_table_GENCODE_vM25.csv.gz",
    feature_id = "GENCODE vM25 GeneSymbol",
    matrix_grain = "Gene by RNA-seq library",
    source_unit = "nonnegative raw integer gene counts",
    result_unit = "TMM CPM",
    loci = nrow(counts),
    libraries = ncol(counts),
    min_count = min(counts),
    max_count = max(counts),
    n_missing_count = sum(is.na(counts)),
    n_all_zero_genes = sum(rowSums(counts) == 0),
    n_duplicate_genes = sum(duplicated(rownames(counts))),
    library_size_min = min(depth),
    library_size_max = max(depth),
    tmm_factor_min = min(factors),
    tmm_factor_max = max(factors),
    stringsAsFactors = FALSE
)
write_tab(schema, file.path(audit_dir, "R10_IMMGEN_SCHEMA.tsv"))
cat("Normalized", nrow(counts), "genes and", ncol(counts), "libraries by TMM in",
    "edgeR", as.character(packageVersion("edgeR")), "\n")
cat("Group medians:", ncol(median_cpm), "populations;", sum(primary), "baseline libraries\n")
cat("TMM factor range:", format(min(factors), digits = 5L), "to",
    format(max(factors), digits = 5L), "\n")
print(sessionInfo())
