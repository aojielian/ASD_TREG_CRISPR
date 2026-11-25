r11_cohorts <- c("GSE121279", "GSE263830", "GSE121251")
r11_counts_name <- c(GSE121279 = "GSE121279_rnaseq-counts-raw.txt.gz",
                     GSE263830 = "GSE263830_Sin3aTreg_gene_count.txt.gz",
                     GSE121251 = "GSE121251_Counts.mRNA.foxp1.txt.gz")
r11_min_samples <- c(GSE121279 = 3L, GSE263830 = 3L, GSE121251 = 2L)
r11_markers <- c("Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18", "Tigit", "Entpd1", "Lrrc32",
                  "Il10", "Tgfb1", "Satb1", "Tcf7", "Tet1", "Tet3")
r11_contrast_ids <- c("F1N", "F1A", "S3", "F2R", "F2A")
r11_reference <- c(GSE121279 = "naive", GSE263830 = "bulk", GSE121251 = "resting")
r11_targets <- c(GSE121279 = "Foxp1", GSE263830 = "Sin3a", GSE121251 = "Foxp1")

r11_read <- function(path) {
    read.delim(path, check.names = FALSE, quote = "", na.strings = c("NA", ""),
               comment.char = "", stringsAsFactors = FALSE)
}

r11_write <- function(data, path) {
    if (grepl("\\.gz$", path)) {
        connection <- gzfile(path, "wt")
        on.exit(close(connection))
        path <- connection
    }
    write.table(data, path, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)
}

r11_read_counts <- function(raw_dir, cohort) {
    path <- file.path(raw_dir, cohort, r11_counts_name[[cohort]])
    connection <- gzfile(path, "rt")
    on.exit(close(connection))
    header <- strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1L]]
    # Konopacki's first physical data field is a row symbol. The header only lists libraries.
    names <- if (cohort == "GSE121279") c("gene_id", header) else header
    data <- read.delim(connection, header = FALSE, col.names = names, check.names = FALSE,
                       stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = "")
    id <- data[[1L]]
    if (cohort == "GSE121251") {
        author_symbol <- data[[2L]]
        counts <- as.matrix(data[-c(1L, 2L)])
    } else {
        author_symbol <- if (cohort == "GSE121279") id else rep(NA_character_, length(id))
        counts <- as.matrix(data[-1L])
    }
    if (anyNA(id) || anyDuplicated(id) || anyDuplicated(colnames(counts)) ||
        !is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) ||
        any(counts %% 1 != 0) || any(counts > .Machine$integer.max)) {
        stop("Raw gene/count schema check failed for ", cohort)
    }
    if (cohort != "GSE121279" && !all(grepl("^ENSMUSG[0-9]+(\\.[0-9]+)?$", id))) {
        stop("Unexpected Ensembl stable gene ID schema in ", cohort)
    }
    rownames(counts) <- id
    storage.mode(counts) <- "integer"
    list(counts = counts, author_symbol = author_symbol, original_header = header)
}

r11_design <- function(samples, cohort) {
    samples$genotype <- factor(samples$genotype, levels = c("control", "KO"))
    if (cohort == "GSE263830") {
        samples$sex <- factor(samples$sex, levels = c("female", "male"))
        formula <- ~ sex + genotype
    } else {
        samples$state <- factor(samples$state, levels = c(r11_reference[[cohort]], "activated"))
        formula <- ~ state + genotype + state:genotype
    }
    rownames(samples) <- samples$sample_id
    design <- model.matrix(formula, samples)
    stopifnot(nrow(design) == nrow(samples), qr(design)$rank == ncol(design))
    if (cohort == "GSE263830") {
        contrasts <- list(S3 = c(0, 0, 1))
        states <- c(S3 = "bulk")
    } else {
        ids <- if (cohort == "GSE121279") c("F1N", "F1A") else c("F2R", "F2A")
        contrasts <- setNames(list(c(0, 0, 1, 0), c(0, 0, 1, 1)), ids)
        states <- setNames(levels(samples$state), ids)
    }
    stopifnot("genotypeKO" %in% colnames(design),
              (cohort == "GSE263830" && identical(colnames(design), c("(Intercept)", "sexmale", "genotypeKO"))) ||
              (cohort != "GSE263830" && identical(colnames(design),
                c("(Intercept)", "stateactivated", "genotypeKO", "stateactivated:genotypeKO"))))
    list(samples = samples, formula = formula, design = design, contrasts = contrasts, states = states)
}
