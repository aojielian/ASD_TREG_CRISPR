args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_rnaseq_metadata.R RAW_FILES STAGE")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(script_dir, "rnaseq_io.R"))
raw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
stage <- normalizePath(args[[2L]], mustWork = TRUE)
audit <- file.path(stage, "01_audit")
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

metadata_info <- data.frame(
    dataset = r11_cohorts,
    author = c("Konopacki et al.", "Christensen et al.", "Ren et al. (H. Wang senior author)"),
    pmid = c("30643266", "39156895", "31125332"),
    doi = c("10.1038/s41590-018-0291-z", "10.3389/fimmu.2024.1444937", "10.1371/journal.pbio.3000270"),
    bioproject = c("PRJNA496510", "PRJNA1099388", "PRJNA496482"),
    paper_url = c("https://pmc.ncbi.nlm.nih.gov/articles/PMC7534899/",
                  "https://pmc.ncbi.nlm.nih.gov/articles/PMC11327135/",
                  "https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000270"),
    source_sections = c("Results: Foxp1-dependent regulation of gene expression in Treg cells; Methods: Mice, RNA Sequencing, RNA-seq; GEO sample extract and growth protocols",
                        "Methods: Analysis of RNA; Fig. 4 caption; GEO sample treatment and extract protocols",
                        "Methods: RNA-Seq analysis; Supplementary Fig. S1; GEO sample treatment and extract protocols"),
    genome = c("GRCm38/mm10", "GRCm38/mm10", "GRCm38/mm10"),
    gene_model = c("Ensembl GRCm38.83; Methods", "Gene-model release not stated in GEO/paper", "Gene-model release not stated in GEO/paper"),
    assay_unit = c("Biological library; number of mice represented in each library not reported",
                   "Pooled biological replicate/library; number of mice per pool not reported",
                   "Biological library; number of mice represented in each library not reported"),
    pooling = c("Spleen/LN material pooled; no per-mouse, cross-state or cross-genotype pool link",
                "Six pooled biological replicates per genotype, three female and three male; no overlapping-pool identifiers",
                "Spleen and peripheral LN source; per-library mouse pooling and cross-state linkage not reported"),
    age = c("8–10 weeks is source default; RNA library-specific ages unreported",
             "RNA library-specific ages unreported",
             "6–8 weeks in SOFT treatment protocol"),
    check = TRUE
)
r11_write(metadata_info, file.path(audit, "R11_RNASEQ_STUDY_METADATA.tsv"))

r11_read_soft <- function(cohort) {
    path <- file.path(raw_dir, cohort, paste0(cohort, "_family.soft.gz"))
    connection <- gzfile(path, "rt")
    on.exit(close(connection))
    lines <- readLines(connection, encoding = "UTF-8")
    starts <- grep("^\\^SAMPLE = GSM", lines)
    stops <- c(starts[-1L] - 1L, length(lines))
    series <- grep("^!Series_relation = BioProject:", lines, value = TRUE)
    bioproject <- regmatches(series, regexpr("PRJNA[0-9]+", series))
    stopifnot(length(bioproject) == 1L && bioproject ==
        metadata_info$bioproject[match(cohort, metadata_info$dataset)])
    pieces <- vector("list", length(starts))
    for (i in seq_along(starts)) {
        lines_here <- lines[starts[[i]]:stops[[i]]]
        attr <- function(field) {
            prefix <- paste0("!Sample_", field, " = ")
            matches <- lines_here[startsWith(lines_here, prefix)]
            substring(matches, nchar(prefix) + 1L)
        }
        characteristic <- function(name) {
            prefix <- paste0(tolower(name), ": ")
            all <- attr("characteristics_ch1")
            ans <- substring(all[startsWith(tolower(all), prefix)], nchar(prefix) + 1L)
            if (length(ans)) paste(ans, collapse = " | ") else NA_character_
        }
        title <- attr("title")
        description <- attr("description")
        relation <- attr("relation")
        biosample <- grep("^BioSample", relation, value = TRUE)
        sra <- grep("^SRA", relation, value = TRUE)
        stopifnot(length(title) == 1L, length(biosample) == 1L, length(sra) == 1L)
        gsm <- sub("^\\^SAMPLE = ", "", lines[[starts[[i]]]])
        include <- cohort != "GSE121279" || "RNA-seq" %in% description
        sample_id <- NA_character_
        replicate <- NA_integer_
        state <- genotype <- sex <- gate <- NA_character_
        exact_key <- "Excluded non-RNA-seq assay"
        if (include && cohort == "GSE121279") {
            desc_id <- grep("^rna_(wt|fl)_(ntr|atr)[1-3]$", description, value = TRUE)
            stopifnot(length(desc_id) == 1L)
            sample_id <- sub("^rna_(wt|fl)_(ntr|atr)([1-3])$", "\\1_\\2_\\3", desc_id)
            replicate <- as.integer(sub(".*_([1-3])$", "\\1", sample_id))
            genotype <- if (startsWith(sample_id, "wt_")) "control" else "KO"
            state <- if (grepl("_ntr_", sample_id)) "naive" else "activated"
            gate <- if (state == "naive") "CD62Lhi YFP+" else "CD44hi YFP+"
            stopifnot(grepl(if (state == "naive") "naïve" else "activated", title, fixed = TRUE),
                      grepl(if (genotype == "control") "Foxp1+/+" else "Foxp1f/f",
                            characteristic("genotype"), fixed = TRUE),
                      endsWith(title, paste0("rep", replicate)))
            sex <- "female"
            exact_key <- "SOFT rna description transformed to count header; title, activation gate and genotype agree"
        }
        if (cohort == "GSE263830") {
            stopifnot(length(description) == 1L)
            sample_id <- if (description == "SIn3aKO_F_2") "Sin3aKO_F_2" else description
            stopifnot(grepl("^(CTR|Sin3aKO)_[FM]_[1-3]$", sample_id))
            genotype <- if (startsWith(sample_id, "CTR")) "control" else "KO"
            state <- "bulk"
            sex <- tolower(characteristic("sex"))
            replicate <- as.integer(sub(".*_([1-3])$", "\\1", sample_id))
            gate <- "CD25+Foxp3+ post CD4+CD25+ enrichment"
            stopifnot(grepl(paste0(", ", sex, ", biol rep ", replicate), title, fixed = TRUE),
                      if (sex == "female") grepl("_F_", sample_id) else grepl("_M_", sample_id),
                      if (genotype == "KO") grepl("Sin3a-/-", characteristic("genotype"), fixed = TRUE)
                          else characteristic("genotype") == "FoxP3YFPcre")
            exact_key <- if (description == sample_id) "Exact SOFT description and count header; title confirms sex/genotype/rep"
                         else "One letter-case typo SIn3aKO_F_2; independently exact title/genotype/sex/rep"
        }
        if (cohort == "GSE121251") {
            stopifnot(grepl("^(WT|Foxp1 cKO) [ar]Treg [1-2]$", title))
            sample_id <- gsub(" ", "_", sub("^Foxp1 ", "", title))
            genotype <- if (startsWith(title, "WT")) "control" else "KO"
            state <- if (grepl("rTreg", title, fixed = TRUE)) "resting" else "activated"
            replicate <- as.integer(sub(".*_([1-2])$", "\\1", sample_id))
            gate <- attr("source_name_ch1")
            stopifnot(grepl(if (state == "resting") "CD44lowCD62Lhigh" else "CD44highCD62Llow",
                            gate, fixed = TRUE),
                      if (genotype == "control") characteristic("genotype") == "Foxp3-YFP-Cre"
                          else grepl("Foxp1fl/fl", characteristic("genotype"), fixed = TRUE))
            exact_key <- "SOFT title transformed to underscore count header; source gate and genotype agree"
        }
        assay <- if (include) "RNA-seq" else "ChIP-seq"
        stopifnot(include || "Foxp1 ChIP-seq" %in% description || "Foxp3 ChIP-seq" %in% description ||
                  grepl("ChIP-seq", title, fixed = TRUE))
        source <- metadata_info[match(cohort, metadata_info$dataset), ]
        pieces[[i]] <- data.frame(dataset = cohort, GSM = gsm, sample_id = sample_id,
            include = include, state = state, genotype = genotype, sex = sex, replicate = replicate,
            assay = assay, soft_title = title, soft_description = paste(description, collapse = " | "),
            genotype_author = characteristic("genotype"), activation_author = characteristic("activation status"),
            gate = gate, source_ch1 = paste(attr("source_name_ch1"), collapse = " | "),
            tissue = if (include) switch(cohort, GSE121279 = "Spleen and lymph nodes",
                GSE263830 = "Secondary lymphoid organs: spleen and lymph nodes",
                GSE121251 = "Spleen and peripheral lymph nodes") else NA_character_,
            mouse_age = if (include) source$age else NA_character_,
            sex_evidence = if (!include) NA_character_ else switch(cohort,
                GSE121279 = "Paper Results specify female Foxp3-YFP-Cre/+ RNA-seq animals",
                GSE263830 = "SOFT sex and paper Fig. 4 caption agree",
                GSE121251 = "No library-level sex reported"),
            experimental_unit = if (include) source$assay_unit else NA_character_,
            pooling_and_pairing = if (include) source$pooling else NA_character_,
            batch = NA_character_, bioproject = bioproject,
            biosample = sub("^BioSample: .*[/]", "", biosample),
            sra_experiment = sub("^.*[?]term=", "", sra), match_evidence = exact_key,
            PMID = source$pmid, DOI = source$doi,
            paper_url = source$paper_url, soft_url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm))
    }
    do.call(rbind, pieces)
}

r11_get_map <- function(data, cohort) {
    id <- rownames(data$counts)
    stable <- if (cohort == "GSE121279") rep(NA_character_, length(id)) else sub("\\.[0-9]+$", "", id)
    if (cohort != "GSE121279") stopifnot(!anyDuplicated(stable))
    symbol <- data$author_symbol
    source <- rep(switch(cohort, GSE121279 = "Author gene symbols; Ensembl GRCm38.83",
                         GSE263830 = "org.Mm.eg.db current stable-ENSMUSG bridge",
                         GSE121251 = "Author GeneName linked to author Ensembl ID"), length(id))
    same_id_symbol_count <- as.integer(!is.na(symbol) & nzchar(symbol))
    symbols_for_id <- symbol
    if (cohort == "GSE263830") {
        org <- AnnotationDbi::select(org.Mm.eg.db, keys = stable, keytype = "ENSEMBL",
                                     columns = "SYMBOL")
        org <- unique(org[!is.na(org$SYMBOL), c("ENSEMBL", "SYMBOL")])
        linked <- split(org$SYMBOL, org$ENSEMBL)
        indexed <- unname(linked[stable])
        same_id_symbol_count <- lengths(indexed)
        symbols_for_id <- vapply(indexed, function(x) if (length(x)) paste(x, collapse = "|")
                                                             else NA_character_, character(1))
        symbol <- ifelse(same_id_symbol_count == 1L, symbols_for_id, NA_character_)
    }
    symbol <- ifelse(is.na(symbol) | !nzchar(symbol), NA_character_, symbol)
    stable_symbol_frequency <- table(symbol)
    symbol_unique <- !is.na(symbol) & as.integer(stable_symbol_frequency[symbol]) == 1L
    symbol_unique[is.na(symbol_unique)] <- FALSE
    n10 <- rowSums(data$counts >= 10L)
    minimum <- r11_min_samples[[cohort]]
    result <- data.frame(dataset = cohort, gene_id = id, stable_ensembl_id = stable,
        author_gene_symbol = data$author_symbol, gene_symbol = symbol,
        symbol_unique = symbol_unique, n_source_symbols = same_id_symbol_count, source_symbols = symbols_for_id,
        mapping_source = source,
        mapping_version = switch(cohort, GSE121279 = "Ensembl GRCm38.83 (original annotation)",
             GSE263830 = paste0("org.Mm.eg.db ", packageVersion("org.Mm.eg.db"),
                       "; AnnotationDbi ", packageVersion("AnnotationDbi")),
             GSE121251 = "Author mm10 GeneName (release unreported)"),
        mapping_status = ifelse(is.na(symbol), ifelse(same_id_symbol_count > 1L, "one ID to multiple symbols", "unmapped symbol"),
                               ifelse(symbol_unique, "unique symbol in count matrix", "symbol occurs for several count IDs")),
        nonzero_samples = rowSums(data$counts > 0L), samples_with_count_ge10 = n10,
        min_samples_required = minimum, pass_count_filter = n10 >= minimum)
    stopifnot(!anyDuplicated(result$gene_id))
    result
}

sample_rows <- list()
maps <- list()
schemas <- list()
design_rows <- list()
for (cohort in r11_cohorts) {
    data <- r11_read_counts(raw_dir, cohort)
    samples <- r11_read_soft(cohort)
    usable <- samples[samples$include, ]
    stopifnot(!anyDuplicated(usable$sample_id), setequal(colnames(data$counts), usable$sample_id))
    usable <- usable[match(colnames(data$counts), usable$sample_id), ]
    stopifnot(identical(colnames(data$counts), usable$sample_id))
    if (cohort == "GSE121279") stopifnot(sum(!samples$include) == 17L)
    levels <- if(cohort == "GSE263830") with(usable, table(genotype, sex)) else with(usable, table(genotype, state))
    stopifnot(all(levels == r11_min_samples[[cohort]]), length(levels) == 4L)
    design <- r11_design(usable, cohort)
    sample_rows[[cohort]] <- rbind(usable, samples[!samples$include, ])
    map <- r11_get_map(data, cohort)
    maps[[cohort]] <- map
    binary_ids <- paste(data$original_header, collapse = "|")
    matrix <- data$counts
    schemas[[cohort]] <- data.frame(dataset = cohort, source_file = r11_counts_name[[cohort]],
        assay = "raw integer gene count", original_header = binary_ids, count_matrix_rows = nrow(matrix),
        count_matrix_libraries = ncol(matrix), GSM_total = nrow(samples), RNA_GSMs = nrow(usable),
        excluded_non_RNA_GSMs = sum(!samples$include), ID_type = switch(cohort,
             GSE121279 = "Exact author gene symbol", GSE263830 = "Unversioned mouse Ensembl gene ID",
             GSE121251 = "Unversioned mouse Ensembl gene ID in column labelled RefSeq"),
        genome = metadata_info$genome[match(cohort, metadata_info$dataset)],
        gene_annotation = metadata_info$gene_model[match(cohort, metadata_info$dataset)],
        duplicate_gene_ids = anyDuplicated(rownames(matrix)), duplicate_samples = anyDuplicated(colnames(matrix)),
        missing_count_values = sum(is.na(matrix)), min_count = min(matrix), max_count = max(matrix),
        zero_count_values = sum(matrix == 0L), all_zero_genes = sum(rowSums(matrix) == 0),
        tested_design_count_filter_n = sum(map$pass_count_filter), unambiguous_symbol_rows = sum(map$symbol_unique),
        rows_with_nonunique_symbol = sum(!is.na(map$gene_symbol) & !map$symbol_unique),
        rows_without_unique_retrievable_symbol = sum(is.na(map$gene_symbol)),
        sample_min_depth = min(colSums(matrix)), sample_max_depth = max(colSums(matrix)))
    design_rows[[cohort]] <- data.frame(dataset = cohort, sample_id = usable$sample_id,
         intercept = design$design[, 1], state_activated = if (cohort == "GSE263830") NA else design$design[, 2],
         sex_male = if (cohort == "GSE263830") design$design[, 2] else NA,
         genotype_KO = design$design[, 3], interaction = if (cohort == "GSE263830") NA else design$design[, 4],
         rank = qr(design$design)$rank, columns = ncol(design$design),
         samples = nrow(design$design), residual_df = nrow(design$design) - ncol(design$design))
}
sample_rows <- do.call(rbind, sample_rows)
stopifnot(!anyDuplicated(sample_rows$GSM), !anyDuplicated(sample_rows$biosample),
          !anyDuplicated(sample_rows$sra_experiment))
r11_write(sample_rows, file.path(audit, "R11_SAMPLE_MAP.tsv"))
r11_write(do.call(rbind, maps), file.path(audit, "R11_ID_MAP.tsv.gz"))
r11_write(do.call(rbind, schemas), file.path(audit, "R11_SCHEMA.tsv"))
r11_write(do.call(rbind, design_rows), file.path(audit, "R11_RNASEQ_DESIGN.tsv"))
r11_write(AnnotationDbi::metadata(org.Mm.eg.db), file.path(audit, "R11_MOUSE_ORGDB_METADATA.tsv"))
