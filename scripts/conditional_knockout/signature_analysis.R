# Score independently frozen Treg and Hallmark references in each RNA-seq cohort.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript signature_analysis.R GEO_DOWNLOAD_FILES OUT")
download <- normalizePath(args[1])
out <- normalizePath(args[2])
audit <- file.path(out, "01_audit")
result_dir <- file.path(out, "03_results")
read_tsv <- function(path, ...) read.delim(path, sep = "\t", check.names = FALSE, quote = "", ...)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_gz <- function(x, path) {
    target <- gzfile(path, "wt", compression = 9L)
    on.exit(close(target))
    write_tsv(x, target)
}
filter_min <- c(GSE121279 = 3, GSE263830 = 3, GSE121251 = 2)
file_names <- c(
    GSE121279 = "GSE121279_rnaseq-counts-raw.txt.gz",
    GSE263830 = "GSE263830_Sin3aTreg_gene_count.txt.gz",
    GSE121251 = "GSE121251_Counts.mRNA.foxp1.txt.gz"
)
contrast_list <- list(
    GSE121279 = c("F1N", "F1A"),
    GSE263830 = "S3",
    GSE121251 = c("F2R", "F2A")
)
all_contrasts <- unlist(contrast_list, use.names = FALSE)
source_member <- read_tsv(gzfile(file.path(result_dir, "R11_SIGNATURE_MEMBERS.tsv.gz")))
defs <- read_tsv(file.path(audit, "R11_SIGNATURE_DEFINITIONS.tsv"))
meta <- read_tsv(file.path(audit, "R11_SAMPLE_MAP.tsv"), colClasses = c(GSM = "character"))
known_design <- read_tsv(file.path(audit, "R11_RNASEQ_DESIGN.tsv"))
id_map <- read_tsv(gzfile(file.path(audit, "R11_ID_MAP.tsv.gz")))
active <- source_member[source_member$included, ]
stopifnot(nrow(defs) == 10L, length(unique(defs$program)) == 8L,
          !anyDuplicated(active[c("program", "component", "mouse_gene")]))
program_order <- unique(defs$program)
components <- unique(defs[c("program", "component")])

read_count_matrix <- function(accession, samples, ids) {
    path <- file.path(download, accession, file_names[[accession]])
    x <- if (accession == "GSE121279") {
        read_tsv(gzfile(path), row.names = 1)
    } else {
        raw <- read_tsv(gzfile(path))
        row.names(raw) <- raw[[1]]
        raw[samples$sample_id]
    }
    stopifnot(setequal(names(x), samples$sample_id))
    x <- as.matrix(x[, samples$sample_id, drop = FALSE])
    is_count <- is.finite(x) & x >= 0 & x == trunc(x)
    stopifnot(all(is_count), !anyDuplicated(rownames(x)),
              all(ids$gene_id %in% rownames(x)))
    x <- x[ids$gene_id, , drop = FALSE]
    filter <- rowSums(x >= 10) >= filter_min[[accession]]
    stopifnot(identical(unname(filter), ids$pass_count_filter))
    x[filter, , drop = FALSE]
}

contrast_statistics <- function(y, design, vector) {
    if (all(is.na(y))) return(rep(NA_real_, 5L))
    fit <- lm.fit(design, y)
    dof <- nrow(design) - fit$rank
    deviation <- sqrt(sum(fit$residuals^2) / dof)
    error <- sqrt(drop(t(vector) %*% solve(crossprod(design), vector))) * deviation
    est <- drop(crossprod(vector, fit$coefficients))
    limits <- est + c(-1, 1) * qt(0.975, dof) * error
    p <- if (error > 0) 2 * pt(-abs(est / error), dof) else NA_real_
    c(est, error, limits, p)
}
select_members <- function(dataset, ids) {
    lookup <- match(active$mouse_gene, ids$gene_symbol)
    unique <- ids$symbol_unique[lookup]
    filtered <- ids$pass_count_filter[lookup]
    status <- ifelse(is.na(lookup), "absent_or_unmapped",
                     ifelse(is.na(unique) | !unique, "ambiguous_gene_symbol",
                            ifelse(filtered, "measured", "count_filter_failed")))
    data.frame(dataset = dataset, program = active$program, component = active$component,
               mouse_gene = active$mouse_gene, gene_id = ifelse(status == "measured", ids$gene_id[lookup], NA_character_),
               status = status)
}

analyse_dataset <- function(dataset) {
    samples <- meta[meta$dataset == dataset & meta$include, c("dataset", "GSM", "sample_id", "state", "genotype", "sex")]
    ids <- id_map[id_map$dataset == dataset, ]
    stopifnot(!anyDuplicated(samples$sample_id), !anyDuplicated(ids$gene_id))
    storage <- read_count_matrix(dataset, samples, ids)
    samples$genotype <- factor(samples$genotype, levels = c("control", "KO"))
    if (dataset == "GSE263830") {
        samples$sex <- factor(samples$sex, levels = c("female", "male"))
        fit_design <- model.matrix(~ sex + genotype, samples)
        wanted <- c("intercept", "sex_male", "genotype_KO")
    } else {
        lower <- if (dataset == "GSE121279") "naive" else "resting"
        samples$state <- factor(samples$state, levels = c(lower, "activated"))
        fit_design <- model.matrix(~ state + genotype + state:genotype, samples)
        wanted <- c("intercept", "state_activated", "genotype_KO", "interaction")
    }
    recorded <- known_design[known_design$dataset == dataset, ]
    prior_mat <- as.matrix(recorded[match(samples$sample_id, recorded$sample_id), wanted])
    stopifnot(identical(dim(fit_design), dim(prior_mat)), all(fit_design == prior_mat),
              qr(fit_design)$rank == ncol(fit_design))
    contrasts <- lapply(contrast_list[[dataset]], function(contrast) {
        vector <- numeric(ncol(fit_design))
        vector[match("genotypeKO", colnames(fit_design))] <- 1
        if (endsWith(contrast, "A")) vector[match("stateactivated:genotypeKO", colnames(fit_design))] <- 1
        vector
    })
    names(contrasts) <- contrast_list[[dataset]]
    # Only a library-scale factor changes within a sample; normalization leaves within-sample rank order unchanged.
    edge <- edgeR::calcNormFactors(edgeR::DGEList(storage), method = "TMM")
    log_cpm <- log2(edgeR::cpm(edge) + 0.5)
    percentile <- apply(log_cpm, 2, function(column) (rank(column, ties.method = "average") - 0.5) / length(column))
    rownames(percentile) <- rownames(storage)
    weights <- limma::voom(edge, design = fit_design, plot = FALSE)
    measurement <- select_members(dataset, ids)
    index <- list()
    measure_n <- list()
    cover_rows <- list()
    for (i in seq_len(nrow(defs))) {
        d <- defs[i, ]
        candidates <- measurement[measurement$program == d$program & measurement$component == d$component, ]
        measured_id <- candidates$gene_id[candidates$status == "measured"]
        stopifnot(all(measured_id %in% rownames(storage)), !anyDuplicated(measured_id),
                  nrow(candidates) == d$frozen_member_n)
        key <- paste(d$program, d$component, sep = ":")
        measure_n[[key]] <- length(measured_id)
        index[[key]] <- match(measured_id, rownames(storage))
        frac <- length(measured_id) / d$frozen_member_n
        cover_rows[[i]] <- data.frame(dataset = dataset, program = d$program, component = d$component,
            source_set = d$source_set, starting_member_n = d$starting_member_n, frozen_member_n = d$frozen_member_n,
            measured_n = length(measured_id), measured_fraction = frac, absent_or_unmapped_n = sum(candidates$status == "absent_or_unmapped"),
            ambiguous_symbol_n = sum(candidates$status == "ambiguous_gene_symbol"),
            count_filter_failed_n = sum(candidates$status == "count_filter_failed"),
            rank_universe_n = nrow(storage), min_samples_count_ge10 = filter_min[[dataset]],
            passes_component_coverage = if (d$program == "Hill_LN") length(measured_id) >= 20 && frac >= 0.5 else length(measured_id) >= 1)
    }
    coverage <- do.call(rbind, cover_rows)
    coverage$passes_score_coverage <- vapply(coverage$program, function(program) {
        all(coverage$passes_component_coverage[coverage$program == program])
    }, logical(1L))
    normalized <- data.frame(dataset = dataset, sample_id = samples$sample_id, GSM = samples$GSM,
        filtered_library_size = edge$samples$lib.size, TMM_factor = edge$samples$norm.factors,
        effective_library_size = edge$samples$lib.size * edge$samples$norm.factors,
        rank_universe_n = nrow(storage))
    scores <- list()
    for (p in program_order) {
        b <- coverage[coverage$program == p, ]
        has_arms <- nrow(b) == 2L
        if (has_arms) {
            c_high <- paste(p, "Treg_high", sep = ":")
            c_low <- paste(p, "Tconv_high", sep = ":")
            high_score <- if (length(index[[c_high]]) > 0) colMeans(percentile[index[[c_high]], , drop = FALSE]) else rep(NA_real_, nrow(samples))
            low_score <- if (length(index[[c_low]]) > 0) colMeans(percentile[index[[c_low]], , drop = FALSE]) else rep(NA_real_, nrow(samples))
            score <- if (all(b$passes_score_coverage)) high_score - low_score else rep(NA_real_, nrow(samples))
            scores[[p]] <- rbind(
                data.frame(samples, program = p, component = "score", score = score, status = ifelse(is.na(score), "insufficient_reference_coverage", "estimated")),
                data.frame(samples, program = p, component = "Treg_high", score = high_score, status = ifelse(is.na(high_score), "no_measured_members", "component_descriptive")),
                data.frame(samples, program = p, component = "Tconv_high", score = low_score, status = ifelse(is.na(low_score), "no_measured_members", "component_descriptive"))
            )
        } else {
            mi <- index[[paste(p, "member_mean", sep = ":")]]
            score <- if (length(mi)) colMeans(percentile[mi, , drop = FALSE]) else rep(NA_real_, nrow(samples))
            scores[[p]] <- data.frame(samples, program = p, component = "score", score = score,
                                      status = ifelse(is.na(score), "no_measured_members", "estimated"))
        }
    }
    scores <- do.call(rbind, scores)
    results <- list()
    camera_rows <- list()
    for (label in names(contrasts)) {
        vector <- contrasts[[label]]
        state <- if (dataset == "GSE263830") rep(TRUE, nrow(samples)) else {
            samples$state == if (endsWith(label, "A")) "activated" else levels(samples$state)[1]
        }
        fits <- lapply(seq_len(nrow(components)), function(i) {
            p <- components$program[i]
            component <- components$component[i]
            if (component == "member_mean") component <- "score"
            ss <- scores[scores$program == p & scores$component == component, ]
            stopifnot(identical(ss$sample_id, samples$sample_id))
            if (component != "score") {
                stats <- contrast_statistics(ss$score, fit_design, vector)
                stats[5] <- NA_real_
            } else {
                stats <- contrast_statistics(ss$score, fit_design, vector)
            }
            ss$component <- component
            b <- coverage[coverage$program == p, ]
            n_each <- setNames(b$measured_n, b$component)
            data.frame(dataset = dataset, contrast = label, program = p, component = component,
                estimate = stats[1], SE = stats[2], CI_low = stats[3], CI_high = stats[4],
                df = nrow(fit_design) - ncol(fit_design), P = stats[5], adjusted_P = NA_real_,
                adjustment = "none", adjustment_family = "component_descriptive",
                n_total = nrow(samples), n_KO = sum(state & samples$genotype == "KO"),
                n_control = sum(state & samples$genotype == "control"),
                n_reference_high = unname(if (nrow(b) == 1L) n_each["member_mean"] else n_each["Treg_high"]),
                n_reference_low = unname(if (nrow(b) == 1L) NA_integer_ else n_each["Tconv_high"]),
                rank_universe_n = nrow(storage), adjustment_family_size = NA_integer_,
                status = if (is.na(stats[1])) unique(ss$status) else if (component == "score") "estimated" else "component_descriptive",
                unit = "within-sample midrank percentile difference, KO minus control"
            )
        })
        main_rows <- lapply(program_order[program_order %in% c("Hill_LN", "ImmGen_spleen")], function(p) {
            ss <- scores[scores$program == p & scores$component == "score", ]
            base <- fits[[match(p, components$program)]]
            stats <- contrast_statistics(ss$score, fit_design, vector)
            base$component <- "score"
            base[c("estimate", "SE", "CI_low", "CI_high", "P")] <- as.list(stats)
            base$status <- if (is.na(stats[1])) unique(ss$status) else "estimated"
            base
        })
        cstats <- do.call(rbind, c(main_rows, fits))
        results[[label]] <- cstats
        run_keys <- coverage$passes_component_coverage & coverage$measured_n >= 2
        paper <- limma::camera(weights, index[run_keys], design = fit_design, contrast = vector,
                               inter.gene.cor = NA_real_, allow.neg.cor = FALSE, sort = FALSE)
        rows <- vector("list", nrow(defs))
        for (i in seq_len(nrow(defs))) {
            p <- defs$program[i]
            component <- defs$component[i]
            key <- paste(p, component, sep = ":")
            camera <- paper[key, ]
            score <- cstats[cstats$program == p & cstats$component == if (component == "member_mean") "score" else component, ]
            direction <- if (is.na(score$estimate)) NA_character_ else if (score$estimate == 0) "No change" else if (score$estimate > 0) "Up" else "Down"
            rows[[i]] <- data.frame(dataset = dataset, contrast = label, program = p, component = component,
                source_set = defs$source_set[i], measured_n = measure_n[[key]], tested_gene_universe_n = nrow(storage),
                n_total = nrow(samples), n_KO = score$n_KO, n_control = score$n_control,
                CAMERA_direction = if (run_keys[i]) camera$Direction else NA_character_,
                CAMERA_P = if (run_keys[i]) camera$PValue else NA_real_,
                adjusted_P = NA_real_, adjustment_family = NA_character_, adjustment_family_size = NA_integer_,
                mean_inter_gene_correlation = if (run_keys[i]) camera$Correlation else NA_real_,
                score_component_direction = direction,
                score_component_effect = score$estimate,
                directional_agreement = if (run_keys[i]) direction == camera$Direction else NA,
                rotations = 0L, method = "limma-voom CAMERA; analytic competitive, directional",
                status = if (run_keys[i]) "estimated" else "insufficient_reference_coverage"
            )
        }
        camera_rows[[label]] <- do.call(rbind, rows)
    }
    list(scores = scores, coverage = coverage, stats = do.call(rbind, results),
         camera = do.call(rbind, camera_rows), normalized = normalized, measurement = measurement)
}

analyses <- lapply(names(filter_min), analyse_dataset)
bind <- function(key) do.call(rbind, lapply(analyses, "[[", key))
effects <- bind("stats")
effects <- effects[order(match(effects$contrast, all_contrasts), match(effects$program, program_order),
                         match(effects$component, c("score", "Treg_high", "Tconv_high"))), ]
score_rows <- effects$component == "score"
family_spec <- list(
    Hill_primary = list(method = "holm", rows = score_rows & effects$program == "Hill_LN" & effects$contrast %in% c("F1N", "F1A", "S3")),
    Hill_independent_Foxp1 = list(method = "holm", rows = score_rows & effects$program == "Hill_LN" & effects$contrast %in% c("F2R", "F2A")),
    ImmGen_five_contrasts = list(method = "BH", rows = score_rows & effects$program == "ImmGen_spleen"),
    six_Hallmark_by_five_contrasts = list(method = "BH", rows = score_rows & startsWith(effects$program, "HALLMARK_"))
)
for (family_name in names(family_spec)) {
    f <- family_spec[[family_name]]
    effects$adjusted_P[f$rows] <- p.adjust(effects$P[f$rows], method = f$method, n = sum(f$rows))
    effects$adjustment[f$rows] <- if (f$method == "holm") "Holm" else f$method
    effects$adjustment_family[f$rows] <- family_name
    effects$adjustment_family_size[f$rows] <- sum(f$rows)
}
stopifnot(sum(score_rows) == 40, all(effects$CI_low[!is.na(effects$CI_low)] <= effects$CI_high[!is.na(effects$CI_high)]))
tests <- bind("camera")
tests$adjustment_family <- ifelse(tests$program == "Hill_LN", "CAMERA_Hill_components",
    ifelse(tests$program == "ImmGen_spleen", "CAMERA_ImmGen_components", "CAMERA_six_Hallmarks"))
for (family_name in unique(tests$adjustment_family)) {
    i <- tests$adjustment_family == family_name
    tests$adjusted_P[i] <- p.adjust(tests$CAMERA_P[i], method = "BH", n = sum(i))
    tests$adjustment_family_size[i] <- sum(i)
}
tests <- tests[order(match(tests$contrast, all_contrasts), match(tests$program, program_order),
                     match(tests$component, c("Treg_high", "Tconv_high", "member_mean"))), ]
write_tsv(effects, file.path(result_dir, "R11_SIGNATURE_CONTRASTS.tsv"))
write_tsv(tests, file.path(result_dir, "R11_SIGNATURE_ROBUSTNESS.tsv"))
write_tsv(bind("scores"), file.path(result_dir, "R11_SIGNATURE_SAMPLE_SCORES.tsv"))
write_tsv(bind("coverage"), file.path(result_dir, "R11_SIGNATURE_COVERAGE.tsv"))
write_tsv(bind("normalized"), file.path(result_dir, "R11_SIGNATURE_TMM_FACTORS.tsv"))
write_gz(bind("measurement"), file.path(result_dir, "R11_SIGNATURE_MEMBER_MEASUREMENT.tsv.gz"))
print(effects[score_rows, c("contrast", "program", "estimate", "CI_low", "CI_high", "P", "adjusted_P", "status")], row.names = FALSE)
cat("\n", R.version.string, "\nedgeR ", as.character(packageVersion("edgeR")),
    "\nlimma ", as.character(packageVersion("limma")), "\n", sep = "")
