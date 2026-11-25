args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: fit_differential_expression.R RAW_FILES STAGE")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(script_dir, "rnaseq_io.R"))
raw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
stage <- normalizePath(args[[2L]], mustWork = TRUE)
audit <- file.path(stage, "01_audit")
results_dir <- file.path(stage, "03_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
set.seed(20261125)
map_all <- r11_read(file.path(audit, "R11_ID_MAP.tsv.gz"))
sample_all <- r11_read(file.path(audit, "R11_SAMPLE_MAP.tsv"))
write_result <- function(x, name) r11_write(x, file.path(results_dir, name))
design_labels <- c(GSE121279 = "~ state + genotype + state:genotype",
                   GSE263830 = "~ sex + genotype", GSE121251 = "~ state + genotype + state:genotype")
track_genes <- c("Foxp3", "Foxp1", "Sin3a", r11_markers, "Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
sample_columns <- c("dataset", "GSM", "sample_id", "state", "genotype", "sex", "replicate")
width <- qnorm(0.975)

get_de <- function(dds, map, model, cohort, contrast, weights, state) {
    disable_if <- list(object = dds, contrast = weights, independentFiltering = FALSE, parallel = FALSE)
    r <- as.data.frame(do.call(DESeq2::results, disable_if))
    uncensored <- as.data.frame(do.call(DESeq2::results, c(disable_if, list(cooksCutoff = FALSE))))
    ids <- match(rownames(r), map$gene_id)
    stopifnot(!anyNA(ids))
    id_extras <- map[ids, c("gene_id", "gene_symbol", "symbol_unique", "mapping_status", "samples_with_count_ge10")]
    n_control <- if (contrast == "S3") sum(model$samples$genotype == "control") else
                         sum(model$samples$genotype == "control" & model$samples$state == state)
    n_KO <- if (contrast == "S3") sum(model$samples$genotype == "KO") else
                         sum(model$samples$genotype == "KO" & model$samples$state == state)
    reason <- ifelse(is.finite(r$pvalue), "",
        ifelse(is.finite(uncensored$pvalue), "DESeq2 Cook's distance P masking",
               ifelse(r$baseMean == 0, "All counts zero", "Non-estimable Wald P")))
    stopifnot(isTRUE(all.equal(r$padj, p.adjust(r$pvalue, method = "BH"), check.attributes = FALSE)))
    data.frame(contrast = contrast, dataset = cohort, state = state, target = r11_targets[[cohort]],
               id_extras, n_control = n_control, n_KO = n_KO, base_mean = r$baseMean,
               log2FC = r$log2FoldChange, SE = r$lfcSE,
               CI_low = r$log2FoldChange - width * r$lfcSE,
               CI_high = r$log2FoldChange + width * r$lfcSE,
               wald_z = r$stat, pvalue = r$pvalue, padj = r$padj,
               pvalue_without_cooks_mask = uncensored$pvalue, pvalue_NA_reason = reason,
               beta_converged = as.data.frame(S4Vectors::mcols(dds))$betaConv,
               dispersion = DESeq2::dispersions(dds))
}

endpoint_rows <- function(symbols, full, map) {
    pieces <- lapply(symbols, function(symbol) {
        candidates <- which(map$gene_symbol == symbol & !is.na(map$gene_symbol))
        eligible <- candidates[map$symbol_unique[candidates]]
        j <- if (length(eligible) == 1L) match(map$gene_id[eligible], full$gene_id) else NA_integer_
        if (is.na(j)) {
            one <- full[1L, ]
            clear <- names(one)[!names(one) %in% c("contrast", "dataset", "state", "target", "n_control", "n_KO")]
            one[clear] <- NA
            one$gene_symbol <- symbol
            if (length(eligible) == 1L) {
                one$gene_id <- map$gene_id[eligible]
                one$samples_with_count_ge10 <- map$samples_with_count_ge10[eligible]
            }
            one$status <- if (length(eligible) == 1L) "Below locked count filter" else
                if (length(candidates)) "Ambiguous gene symbol in count matrix" else "Gene absent or unmapped"
        } else {
            one <- full[j, ]
            one$status <- if (is.finite(one$pvalue)) "Estimable" else one$pvalue_NA_reason
        }
        one
    })
    do.call(rbind, pieces)
}

save_samples <- function(dds, map, raw_counts, sample, cohort) {
    all_counts <- sweep(raw_counts, 2L, DESeq2::sizeFactors(dds), "/")
    do.call(rbind, lapply(track_genes, function(symbol) {
        eligible <- which(map$gene_symbol == symbol & map$symbol_unique)
        raw_value <- norm_value <- rep(NA_real_, nrow(sample))
        id <- NA_character_
        status <- "Source count absent or ambiguous"
        if (length(eligible) == 1L) {
            id <- map$gene_id[eligible]
            raw_value <- raw_counts[id, ]
            norm_value <- all_counts[id, ]
            status <- if (map$pass_count_filter[eligible]) "Passed locked count filter" else "Below locked count filter"
        }
        row <- sample[, sample_columns, drop = FALSE]
        data.frame(row, gene_symbol = symbol, gene_id = id, raw_count = raw_value,
            normalized_count = norm_value, log2_norm_count_plus1 = log2(norm_value + 1),
            status = status, normalization = "Study-specific DESeq2 median-ratio size factor")
    }))
}

all_de <- list()
all_sample_counts <- list()
sensitivity <- list()
edge_genes <- list()
qc <- list()
qc_pca <- list()
qc_pairs <- list()
model_audit <- list()
comparisons <- list()
interaction <- list()
for (cohort in r11_cohorts) {
    message("Running ", cohort)
    raw <- r11_read_counts(raw_dir, cohort)$counts
    sample <- sample_all[sample_all$dataset == cohort & sample_all$include, ]
    sample <- sample[match(colnames(raw), sample$sample_id), ]
    map <- map_all[map_all$dataset == cohort, ]
    map <- map[match(rownames(raw), map$gene_id), ]
    stopifnot(identical(sample$sample_id, colnames(raw)), identical(rownames(raw), map$gene_id))
    keep <- rowSums(raw >= 10L) >= r11_min_samples[[cohort]]
    stopifnot(identical(keep, setNames(map$pass_count_filter, map$gene_id)))
    counts <- raw[keep, , drop = FALSE]
    design <- r11_design(sample, cohort)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = design$samples, design = design$formula)
    dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", betaPrior = FALSE,
                         minReplicatesForReplace = Inf, parallel = FALSE)
    stopifnot(length(DESeq2::resultsNames(dds)) == ncol(design$design),
              all(is.finite(DESeq2::sizeFactors(dds))), all(DESeq2::sizeFactors(dds) > 0))
    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design$design, robust = TRUE)
    edge <- edgeR::glmQLFit(y, design$design, robust = TRUE)
    meta <- as.data.frame(S4Vectors::mcols(dds))
    sf <- DESeq2::sizeFactors(dds)
    dds_fc <- DESeq2::counts(dds, normalized = TRUE)
    model_audit[[cohort]] <- data.frame(dataset = cohort, formula = design_labels[[cohort]],
        coefficient_order = paste(DESeq2::resultsNames(dds), collapse = " | "),
        matrix_columns = paste(colnames(design$design), collapse = " | "),
        model_rank = qr(design$design)$rank, n_libraries = ncol(counts), n_source_gene_rows = nrow(raw),
        n_count_filtered_gene_rows = nrow(counts), n_beta_nonconverged = sum(!meta$betaConv, na.rm = TRUE),
        n_dispersion_outliers = sum(meta$dispOutlier, na.rm = TRUE),
        fit_type_requested = "parametric", fit_type_used = attr(DESeq2::dispersionFunction(dds), "fitType"),
        DESeq2 = as.character(packageVersion("DESeq2")), edgeR = as.character(packageVersion("edgeR")),
        outlier_replacement = "Disabled by minReplicatesForReplace=Inf",
        cooks_p_mask = "DESeq2 default; reason in contrast table",
        independent_filter = FALSE, control = "control", sign = "KO minus control")
    all_sample_counts[[cohort]] <- save_samples(dds, map, raw, design$samples, cohort)

    # Normalized sample QC is descriptive and retains every independent library.
    log_cpm <- log2(edgeR::cpm(y, normalized.lib.sizes = TRUE) + 0.5)
    pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
    explained <- pca$sdev^2 / sum(pca$sdev^2)
    pc <- pca$x[, 1:2, drop = FALSE]
    # Canonicalize arbitrary eigenvector sign without changing coordinates' relative geometry.
    for (j in 1:2) {
        biggest <- which.max(abs(pca$rotation[, j]))
        pc[, j] <- pc[, j] * sign(pca$rotation[biggest, j])
    }
    qc_pca[[cohort]] <- data.frame(design$samples[, sample_columns], PC1 = pc[, 1], PC2 = pc[, 2],
        pc1_variance_fraction = explained[1], pc2_variance_fraction = explained[2],
        method = "Centred log2(TMM CPM+0.5), all count-filtered genes, no per-gene scaling")
    corr <- cor(log_cpm)
    dist <- as.matrix(dist(t(log_cpm)))
    pairs <- expand.grid(sample_i = colnames(log_cpm), sample_j = colnames(log_cpm), stringsAsFactors = FALSE)
    qc_pairs[[cohort]] <- data.frame(dataset = cohort, pairs,
        pearson_correlation = as.vector(corr), euclidean_distance = as.vector(dist))
    q_counts <- t(apply(raw, 2L, quantile, probs = c(.25, .50, .75, .95, .99), names = FALSE))
    qc[[cohort]] <- data.frame(design$samples[, sample_columns],
        raw_library_size = colSums(raw), filtered_library_size = colSums(counts),
        size_factor = sf, TMM_factor = y$samples$norm.factors,
        effective_TMM_library_size = y$samples$lib.size * y$samples$norm.factors,
        zero_fraction = colMeans(raw == 0L), fraction_gene_rows_nonzero = colMeans(raw > 0L),
        raw_count_q25 = q_counts[, 1], raw_count_median = q_counts[, 2],
        raw_count_q75 = q_counts[, 3], raw_count_q95 = q_counts[, 4], raw_count_q99 = q_counts[, 5],
        status = "All assay-matched libraries retained")

    for (contrast in names(design$contrasts)) {
        weights <- design$contrasts[[contrast]]
        state <- design$states[[contrast]]
        full <- get_de(dds, map, design, cohort, contrast, weights, state)
        write_result(full, paste0("R11_DE_", contrast, ".tsv.gz"))
        all_de[[contrast]] <- endpoint_rows(c("Foxp3", r11_targets[[cohort]], r11_markers), full, map)
        et <- edgeR::topTags(edgeR::glmQLFTest(edge, contrast = weights), n = Inf, sort.by = "none")$table
        j <- match(full$gene_id, rownames(et))
        stopifnot(!anyNA(j))
        ee <- data.frame(contrast = contrast, dataset = cohort, gene_id = full$gene_id,
             gene_symbol = full$gene_symbol, symbol_unique = full$symbol_unique,
             DESeq2_log2FC = full$log2FC, DESeq2_P = full$pvalue, DESeq2_BH = full$padj,
             edgeR_log2FC = et$logFC[j], edgeR_logCPM = et$logCPM[j], edgeR_QL_F = et$F[j],
             edgeR_P = et$PValue[j], edgeR_BH = et$FDR[j])
        common <- is.finite(ee$DESeq2_log2FC) & is.finite(ee$edgeR_log2FC)
        ee$direction_agrees <- ifelse(common, sign(ee$DESeq2_log2FC) == sign(ee$edgeR_log2FC), NA)
        edge_genes[[contrast]] <- ee
        both_fdr <- !is.na(ee$DESeq2_BH) & !is.na(ee$edgeR_BH) & ee$DESeq2_BH < .05 & ee$edgeR_BH < .05
        fox <- ee[which(ee$symbol_unique & ee$gene_symbol == "Foxp3"), ]
        stopifnot(nrow(fox) == 1L)
        sensitivity[[contrast]] <- data.frame(
            contrast = contrast, dataset = cohort, state = state, target = r11_targets[[cohort]],
            count_filtered_gene_rows = nrow(full), n_DESeq2_finite_P = sum(is.finite(full$pvalue)),
            n_cooks_masked = sum(full$pvalue_NA_reason == "DESeq2 Cook's distance P masking", na.rm = TRUE),
            n_other_na_p = sum(!is.finite(full$pvalue) & full$pvalue_NA_reason != "DESeq2 Cook's distance P masking"),
            n_finite_both_LFC = sum(common),
            n_direction_agree = sum(ee$direction_agrees, na.rm = TRUE),
            direction_agreement = mean(ee$direction_agrees, na.rm = TRUE),
            spearman_log2FC = cor(ee$DESeq2_log2FC[common], ee$edgeR_log2FC[common], method = "spearman"),
            DESeq2_BH_below_0p05 = sum(ee$DESeq2_BH < .05, na.rm = TRUE),
            edgeR_BH_below_0p05 = sum(ee$edgeR_BH < .05, na.rm = TRUE),
            both_BH_below_0p05 = sum(both_fdr),
            both_BH_below_0p05_same_direction = sum(both_fdr & ee$direction_agrees, na.rm = TRUE),
            Foxp3_DESeq2_log2FC = fox$DESeq2_log2FC, Foxp3_DESeq2_P = fox$DESeq2_P,
            Foxp3_DESeq2_BH = fox$DESeq2_BH, Foxp3_edgeR_log2FC = fox$edgeR_log2FC,
            Foxp3_edgeR_P = fox$edgeR_P, Foxp3_edgeR_BH = fox$edgeR_BH,
            Foxp3_direction_agrees = fox$direction_agrees)
        if (cohort == "GSE121279") comparisons[[contrast]] <- full
    }
    if (cohort != "GSE263830") {
        code <- if (cohort == "GSE121279") "F1_INTERACTION" else "F2_INTERACTION"
        full <- get_de(dds, map, design, cohort, code, c(0, 0, 0, 1), "interaction difference of KO effects")
        full$n_control <- NA_integer_
        full$n_KO <- NA_integer_
        interaction[[cohort]] <- endpoint_rows(c("Foxp3", "Foxp1", r11_markers), full, map)
    }
}

result_order <- r11_contrast_ids
endpoints <- do.call(rbind, all_de[result_order])
fox <- endpoints[endpoints$gene_symbol == "Foxp3", ]
fox$contrast_family <- ifelse(fox$contrast %in% c("F1N", "F1A", "S3"), "Core Holm (3)", "Auxiliary Holm (2)")
fox$adjusted_p <- ave(fox$pvalue, fox$contrast_family,
    FUN = function(x) p.adjust(x, method = "holm", n = length(x)))
markers <- endpoints[endpoints$gene_symbol %in% r11_markers, ]
markers$contrast_family <- "All measurable fixed marker x contrast tests (BH)"
markers$adjusted_p <- p.adjust(markers$pvalue, method = "BH")
target <- endpoints[endpoints$gene_symbol == endpoints$target, ]
target$contrast_family <- "Descriptive target-RNA QC"
target$adjusted_p <- NA_real_
write_result(fox, "R11_FOXP3_ENDPOINT.tsv")
write_result(markers, "R11_MARKER_EFFECTS.tsv")
write_result(target, "R11_TARGET_QC.tsv")
write_result(do.call(rbind, all_sample_counts), "R11_SAMPLE_EXPRESSION.tsv")
write_result(do.call(rbind, sensitivity[result_order]), "R11_DE_METHOD_SENSITIVITY.tsv")
write_result(do.call(rbind, edge_genes[result_order]), "R11_DE_METHOD_SENSITIVITY_GENES.tsv.gz")
write_result(do.call(rbind, qc), "R11_QC_SAMPLES.tsv")
write_result(do.call(rbind, qc_pca), "R11_QC_SAMPLE_PCA.tsv")
write_result(do.call(rbind, qc_pairs), "R11_QC_SAMPLE_DISTANCES.tsv.gz")
write_result(do.call(rbind, model_audit), "R11_DE_MODEL_DIAGNOSTICS.tsv")
write_result(do.call(rbind, interaction), "R11_DE_INTERACTION_MARKERS.tsv")

author_path <- file.path(raw_dir, "GSE121279", "GSE121279_rnaseq-foxp1fl-over-wt-results.txt.gz")
if (file.exists(author_path)) {
    author <- read.delim(author_path, check.names = FALSE, quote = "", comment.char = "",
                         strip.white = TRUE, na.strings = c("NA", ""))
    stopifnot(all(vapply(author[-1L], is.numeric, logical(1))))
    stopifnot(!anyDuplicated(author$gene), setequal(author$gene, rownames(r11_read_counts(raw_dir, "GSE121279")$counts)))
    author_list <- lapply(names(comparisons), function(contrast) {
        ours <- comparisons[[contrast]]
        state <- if (contrast == "F1N") "naive" else "act"
        raw_author <- author[match(ours$gene_id, author$gene), ]
        stopifnot(identical(raw_author$gene, ours$gene_id))
        original_effect <- raw_author[[paste0("log2FC.", state)]]
        original_q <- raw_author[[paste0("padj.", state)]]
        original_base <- raw_author[[paste0("base.", state)]]
        agree <- ifelse(is.finite(original_effect) & is.finite(ours$log2FC),
                       sign(original_effect) == sign(ours$log2FC), NA)
        data.frame(contrast = contrast, gene_id = ours$gene_id, gene_symbol = ours$gene_symbol,
            r11_log2FC = ours$log2FC, r11_pvalue = ours$pvalue, r11_BH = ours$padj,
            author_log2FC = original_effect, author_FDR = original_q, author_baseMean = original_base,
            direction_identical = agree, effect_difference_r11_minus_author = ours$log2FC - original_effect)
    })
    audit_author <- do.call(rbind, author_list)
    write_result(audit_author, "R11_AUTHOR_RESULTS_GENE_PAIRS.tsv.gz")
    overview <- do.call(rbind, lapply(author_list, function(x) {
        shared <- !is.na(x$direction_identical)
        both_1pct <- !is.na(x$author_FDR) & !is.na(x$r11_BH) & x$author_FDR < .01 & x$r11_BH < .01
        one <- x[x$gene_id == "Foxp3", ]
        original_naive <- x$contrast[1L] == "F1N"
        official_table_n <- sum(author[[if (original_naive) "padj.naive" else "padj.act"]] < .01, na.rm = TRUE)
        paper_text_n <- if (original_naive) 1091L else 2353L
        data.frame(contrast = x$contrast[1L], author_gene_rows = nrow(author),
            author_native_FDR_below_0p01 = official_table_n,
            paper_reported_FDR_below_0p01 = paper_text_n,
            table_minus_paper_reported_count = official_table_n - paper_text_n,
            published_source_count_note = "GEO table values are rounded; its count and the article Methods count differ, reason unresolved",
            count_filtered_gene_rows = nrow(x), both_finite_effect = sum(shared),
            spearman = cor(x$r11_log2FC[shared], x$author_log2FC[shared], method = "spearman"),
            pearson = cor(x$r11_log2FC[shared], x$author_log2FC[shared]),
            sign_agreement_n = sum(x$direction_identical, na.rm = TRUE),
            sign_agreement_fraction = mean(x$direction_identical, na.rm = TRUE), strict_FDR_overlap_0p01 = sum(both_1pct),
            strict_FDR_overlap_same_direction_n = sum(both_1pct & x$direction_identical, na.rm = TRUE),
            author_FDR_below_0p01_in_shared = sum(x$author_FDR < .01, na.rm = TRUE),
            r11_BH_below_0p01 = sum(x$r11_BH < .01, na.rm = TRUE),
            r11_BH_below_0p05 = sum(x$r11_BH < .05, na.rm = TRUE),
            Foxp3_author_log2FC = one$author_log2FC, Foxp3_author_FDR = one$author_FDR,
            Foxp3_r11_log2FC = one$r11_log2FC, Foxp3_r11_BH = one$r11_BH,
            join = "Exact author gene symbol, same sign KO−WT",
            source_model = "Pairwise DESeq2 1.16.1, article Methods Bioinformatic Analyses",
            comparability = "R11 interaction model uses joint normalization, dispersion estimation and a locked count filter; source table is rounded",
            PMCID = "PMC7534899", status = "Compared")
    }))
    write_result(overview, "R11_AUTHOR_RESULTS_COMPARISON.tsv")
} else {
    write_result(data.frame(contrast = c("F1N", "F1A"), status = "Author GEO result file unavailable"),
                 "R11_AUTHOR_RESULTS_COMPARISON.tsv")
}
writeLines(capture.output(sessionInfo()), file.path(stage, "08_logs/R11_RNASEQ_R_SESSION.txt"))
print(fox[, c("contrast", "gene_symbol", "n_control", "n_KO", "log2FC", "CI_low", "CI_high", "pvalue", "padj", "adjusted_p")],
      row.names = FALSE, digits = 8)
