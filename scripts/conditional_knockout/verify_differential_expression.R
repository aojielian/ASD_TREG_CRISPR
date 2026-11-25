args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: verify_differential_expression.R RAW_FILES STAGE")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(script_dir, "rnaseq_io.R"))
raw_dir <- normalizePath(args[1L], mustWork = TRUE)
stage <- normalizePath(args[2L], mustWork = TRUE)
result <- function(name) r11_read(file.path(stage, "03_results", name))
sample_map <- r11_read(file.path(stage, "01_audit", "R11_SAMPLE_MAP.tsv"))
fx <- result("R11_FOXP3_ENDPOINT.tsv")
qc <- result("R11_QC_SAMPLES.tsv")
sample_expression <- result("R11_SAMPLE_EXPRESSION.tsv")
markers <- result("R11_MARKER_EFFECTS.tsv")
qa <- list()
records <- list()
check <- function(item, passed, detail) {
    qa[[length(qa) + 1L]] <<- data.frame(check = item, status = if (isTRUE(passed)) "PASS" else "FAIL",
                                         evidence = detail)
}
approx <- function(a, b, rtol = 2e-5, atol = 2e-7) {
    is.finite(a) & is.finite(b) & (abs(a - b) <= atol + rtol * abs(b))
}
for (i in seq_len(nrow(fx))) {
    f <- fx[i, ]
    cohort <- f$dataset
    raw <- r11_read_counts(raw_dir, cohort)$counts
    source <- sample_map[sample_map$dataset == cohort & sample_map$include, ]
    source <- source[match(colnames(raw), source$sample_id), ]
    qc0 <- qc[qc$dataset == cohort, ]
    qc0 <- qc0[match(source$sample_id, qc0$sample_id), ]
    # Rebuild the design manually as 0/1 variables from the raw sample author labels.
    ko <- as.integer(source$genotype == "KO")
    activated <- as.integer(source$state == "activated")
    male <- as.integer(source$sex == "male")
    sex_model <- f$contrast == "S3"
    x <- if (sex_model) cbind(intercept = 1, male = male, ko = ko) else
        cbind(intercept = 1, activated = activated, ko = ko, activated_x_ko = activated * ko)
    v <- if (sex_model) c(0, 0, 1) else c(0, 0, 1, as.integer(f$state == "activated"))
    y <- as.numeric(raw[f$gene_id, ])
    # Fixed-dispersion negative binomial Fisher scoring is an independent fit.
    family <- MASS::negative.binomial(theta = 1 / f$dispersion, link = "log")
    g <- glm.fit(x, y, offset = log(qc0$size_factor), family = family,
                 control = glm.control(epsilon = 1e-12, maxit = 100L))
    mu <- g$fitted.values
    information <- crossprod(x, (mu / (1 + f$dispersion * mu)) * x)
    effect <- drop(v %*% g$coefficients) / log(2)
    se <- sqrt(drop(v %*% solve(information, v))) / log(2)
    p <- 2 * pnorm(-abs(effect / se))
    interval <- effect + c(-1, 1) * qnorm(.975) * se
    e <- sample_expression[sample_expression$dataset == cohort & sample_expression$gene_symbol == "Foxp3", ]
    e <- e[match(source$sample_id, e$sample_id), ]
    target <- if (sex_model) rep(TRUE, nrow(source)) else source$state == f$state
    check(paste(f$contrast, "NB coefficients, Wald SE, CI and two-sided P"),
          g$converged && all(approx(c(effect, se, interval, p), c(f$log2FC, f$SE, f$CI_low, f$CI_high, f$pvalue))),
          "Independent fixed-dispersion unpenalized NB fit and Fisher-information covariance")
    check(paste(f$contrast, "KO−control estimate sign / sample counts"),
          sign(effect) == sign(f$log2FC) && sum(target & ko == 1L) == f$n_KO &&
          sum(target & ko == 0L) == f$n_control &&
          all(approx(y, e$raw_count, 0, 0)) && all(approx(y / qc0$size_factor, e$normalized_count, 1e-13, 1e-13)),
          "One-to-one count columns, independent KO/state/sex design, offset-normalized counts and n")
    d <- result(paste0("R11_DE_", f$contrast, ".tsv.gz"))
    finite <- is.finite(d$pvalue)
    ord <- order(d$pvalue[finite])
    values <- d$pvalue[finite][ord]
    n <- length(values)
    independent_bh <- pmin(1, rev(cummin(rev(values * n / seq_len(n)))))
    check(paste(f$contrast, "All-background BH / count rule"), nrow(d) == sum(rowSums(raw >= 10L) >= r11_min_samples[[cohort]]) &&
          all(approx(d$padj[finite][ord], independent_bh, 1e-13, 1e-14)) &&
          all(is.na(d$padj[!finite])), "Independent BH rank calculation and raw count-filter size")
    records[[f$contrast]] <- data.frame(contrast = f$contrast, gene_id = f$gene_id, dispersion = f$dispersion,
        independent_log2FC = effect, independent_SE = se, independent_CI_low = interval[1L],
        independent_CI_high = interval[2L], independent_P = p,
        DESeq2_log2FC = f$log2FC, DESeq2_SE = f$SE, DESeq2_P = f$pvalue,
        max_abs_difference = max(abs(c(effect, se, interval, p) -
                                      c(f$log2FC, f$SE, f$CI_low, f$CI_high, f$pvalue))))
}
for (family in unique(fx$contrast_family)) {
    a <- fx[fx$contrast_family == family, ]
    finite <- which(is.finite(a$pvalue))
    j <- finite[order(a$pvalue[finite])]
    independent_holm <- pmin(1, cummax((nrow(a) - seq_along(j) + 1) * a$pvalue[j]))
    check(paste(family, "independent Holm"),
        all(approx(independent_holm, a$adjusted_p[j], 1e-13, 1e-14)) &&
        all(is.na(a$adjusted_p[-j])) , "Fixed family sizes incorporating unavailable endpoints")
}
j <- which(is.finite(markers$pvalue))
j <- j[order(markers$pvalue[j])]
marker_bh <- pmin(1, rev(cummin(rev(markers$pvalue[j] * length(j) / seq_along(j)))))
check("Fixed marker family sizes and adjustment", nrow(markers) == 65L &&
      all(approx(marker_bh, markers$adjusted_p[j], 1e-13, 1e-14)),
      "13 fixed markers × five contrasts including unavailable endpoint rows")
included <- sample_map[sample_map$include, ]
check("No repeated GEO/BioSample/SRA source across studies", nrow(included) == 32L &&
    !anyDuplicated(included$GSM) && !anyDuplicated(included$biosample) &&
    !anyDuplicated(included$sra_experiment) && length(unique(included$bioproject)) == 3L,
    "No identical public sample IDs; individual/pool reuse cannot be excluded without identifiers")
dir.create(file.path(stage, "07_qa"), showWarnings = FALSE)
r11_write(do.call(rbind, qa), file.path(stage, "07_qa", "R11_RNASEQ_INDEPENDENT_QA.tsv"))
r11_write(do.call(rbind, records), file.path(stage, "07_qa", "R11_FOXP3_NB_RECOMPUTATION.tsv"))
print(do.call(rbind, qa)[, c("check", "status")], row.names = FALSE)
if (any(vapply(qa, function(x) x$status != "PASS", logical(1)))) quit(status = 1L)
