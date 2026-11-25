#!/usr/bin/env Rscript
# Independent numerical check of the two frozen R9 discordant Firth models.
# Uses weighted QR and finite-difference PORT optimization in base R.
# No add-on package fit is implied.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
    stop("Usage: Rscript verify_discordant_firth.R DATA_ROOT OUTPUT_DIR [PRODUCTION_RESULTS_DIR]")
}
data_root <- args[[1L]]
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(path) {
    read.delim(path, check.names = FALSE, quote = "", comment.char = "",
               na.strings = c("", "NA"), stringsAsFactors = FALSE)
}
boolean <- function(x) tolower(as.character(x)) %in% c("true", "1")
scale_population <- function(x) {
    d <- x - mean(x)
    d / sqrt(mean(d ^ 2))
}

background <- read_tsv(file.path(data_root, "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz"))
chen <- background[background$screen == "CHEN" & boolean(background$measured), ]
common <- read_tsv(file.path(data_root, "data/processed/chen_shifrut_common_universe.tsv"))
chen_key <- toupper(chen$approved_symbol)
common_key <- toupper(common$approved_symbol)
stopifnot(length(chen_key) == 18966L, !anyDuplicated(chen_key),
          length(common_key) == 18615L, !anyDuplicated(common_key))
indices <- match(common_key, chen_key)
stopifnot(!anyNA(indices))
joined <- chen[indices, ]
score1 <- as.numeric(boolean(common$SFARI_SCORE1))
fo <- as.numeric(boolean(common$CHEN_HIT_BOOL))
tcr <- as.numeric(boolean(common$SHIFRUT_HIT))
D <- fo - tcr
stopifnot(all(D == common$D), sum(score1) == 242,
          all(score1 == boolean(joined$sfari_score1)),
          all(fo == boolean(joined$author_hit)),
          all(common$CHEN_DIRECTION == joined$FOXP3_REGULATOR_DIRECTION),
          sum(common$approved_symbol != joined$approved_symbol) == 192L)
cc <- joined$class_label %in% c("NEITHER_LISTED", "TF_ONLY", "EPIFACTOR_ONLY", "BOTH_LISTED") &
      tolower(as.character(joined$library_member)) %in% c("true", "false", "0", "1") &
      is.finite(joined$HEALTHY_CD4_MEAN_LOGCPM) & is.finite(joined$LOEUF_v2_1_1)
take <- cc & D != 0
symbol <- joined$approved_symbol[take]
selected <- joined[take, ]
asd <- score1[take]
y <- as.numeric(D[take] == 1)
stopifnot(length(y) == 930L, sum(asd) == 19L, sum(y * asd) == 11L, sum(y) == 176L)
stopifnot(all(selected$TF_LISTED %in% c(0, 1)), all(selected$EPIFACTORS_GENES_LISTED %in% c(0, 1)))
nuisance <- cbind(
    TF = selected$TF_LISTED,
    EpiFactors = selected$EPIFACTORS_GENES_LISTED,
    expression = selected$HEALTHY_CD4_MEAN_LOGCPM,
    LOEUF = selected$LOEUF_v2_1_1,
    targeted_library = as.numeric(boolean(selected$library_member))
)
droppable <- apply(nuisance, 2L, function(x) length(unique(x)) == 1)
base_design <- cbind(intercept = 1, asd = asd)
design <- cbind(base_design, nuisance[, !droppable, drop = FALSE])
for (name in intersect(c("expression", "LOEUF"), colnames(design))) {
    design[, name] <- scale_population(design[, name])
}
write.table(data.frame(approved_symbol = symbol, score1 = asd, FOXP3_only = y),
            file.path(output_dir, "INDEPENDENT_FIRTH_DISCORDANT_MEMBERS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

finite_difference <- function(f, par) {
    vapply(seq_along(par), function(i) {
        delta <- rep(0, length(par))
        delta[i] <- 2e-5
        (f(par + delta) - f(par - delta)) / (2 * delta[i])
    }, numeric(1L))
}

minimize <- function(f, start) {
    # PORT optimization and a finite-difference Newton refinement have no
    # dependency on the analytic gradient or SciPy optimizer used in production.
    fit <- nlminb(start, f, control = list(iter.max = 400L, eval.max = 2000L, rel.tol = 1e-13))
    best <- list(par = fit$par, objective = fit$objective, convergence = fit$convergence,
                 optimizer = "nlminb_numeric")
    gradient <- function(p) finite_difference(f, p)
    score <- gradient(best$par)
    for (iteration in seq_len(8L)) {
        residual <- max(abs(score))
        if (residual <= 1e-6) break
        H <- optimHess(best$par, f, gr = gradient)
        step <- tryCatch(solve(H, score), error = function(e) rep(Inf, length(score)))
        if (!all(is.finite(step))) break
        accepted <- FALSE
        for (fraction in c(1, .5, .25, .125)) {
            next_par <- best$par - fraction * step
            next_value <- f(next_par)
            next_score <- gradient(next_par)
            if (is.finite(next_value) && next_value <= best$objective + 1e-10 &&
                max(abs(next_score)) < residual) {
                best$par <- next_par
                best$objective <- next_value
                best$optimizer <- "nlminb_numeric_and_Newton"
                score <- next_score
                accepted <- TRUE
                break
            }
        }
        if (!accepted) break
    }
    best$residual <- max(abs(score))
    best
}

fit <- function(x, model) {
    qr_rank <- qr(x)$rank
    stopifnot(qr_rank == ncol(x))
    negative_ll <- function(theta) {
        eta <- c(x %*% theta)
        w_half <- sqrt(plogis(eta) * plogis(-eta))
        r <- qr.R(qr(x * w_half, LAPACK = TRUE))
        diagonals <- abs(diag(r))
        if (any(diagonals <= 0) || any(!is.finite(diagonals))) return(1e100)
        loglik <- sum(y * eta - pmax(eta, 0) - log1p(exp(-abs(eta))))
        -loglik - sum(log(diagonals))
    }
    initial <- rep(0, ncol(x))
    initial[1L] <- qlogis((sum(y) + .5) / (length(y) + 1))
    optimum <- minimize(negative_ll, initial)
    initial[2L] <- 1
    shifted <- minimize(negative_ll, initial)
    if (max(optimum$residual, shifted$residual) > 2e-5 ||
        abs(optimum$objective - shifted$objective) > 1e-8 ||
        max(abs(optimum$par - shifted$par)) > 2e-5) {
        stop(sprintf("%s did not reach an independently verified full fit", model))
    }
    pmin <- optimum$objective
    mid <- optimum$par[2L]
    largest_profile_score <- 0
    profile_lr <- function(beta) {
        constrained <- function(free) negative_ll(append(free, beta, after = 1))
        prof <- minimize(constrained, optimum$par[-2L])
        largest_profile_score <<- max(largest_profile_score, prof$residual)
        if (prof$residual > 2e-5 || prof$objective < pmin - 1e-8) {
            stop(sprintf("%s has an invalid profile fit at log OR %.9f, score %.7g, ell gap %.7g", model, beta, prof$residual, prof$objective - pmin))
        }
        max(0, 2 * (prof$objective - pmin))
    }
    target <- qchisq(.95, 1)
    roots <- vapply(c(-1, 1), function(side) {
        choices <- c(.25, .5, 1, 2, 4, 8, 12, 20)
        for (step in choices) {
            test <- mid + side * step
            if (profile_lr(test) >= target) {
                return(uniroot(function(v) profile_lr(v) - target, sort(c(mid, test)),
                               tol = 1e-9)$root)
            }
        }
        stop(sprintf("%s has no bracketed profile endpoint", model))
    }, numeric(1L))
    raw_p <- pchisq(profile_lr(0), 1, lower.tail = FALSE)
    list(model = model, n = nrow(x), cases = sum(asd), plus_cases = sum(asd * y),
         minus_cases = sum(asd * (1 - y)), or = exp(mid), lo = exp(roots[1L]),
         hi = exp(roots[2L]), raw_p = raw_p, rank = qr_rank, full_score = optimum$residual,
         second_score = shifted$residual, profile_score = largest_profile_score,
         optimizer = optimum$optimizer, optimizer_code = optimum$convergence)
}

fits <- list(fit(base_design, "DISCORDANT_UNADJUSTED_FIRTH"),
             fit(design, "DISCORDANT_ADJUSTED_FIRTH"))
comparison <- if (length(args) == 3L) {
    read_tsv(file.path(args[[3L]], "PHENOTYPE_PREFERENCE_R9.tsv"))
} else {
    NULL
}
audit <- lapply(fits, function(f) {
    pp <- c(NA_real_, NA_real_, NA_real_, NA_real_)
    status <- "FROZEN_INDEPENDENT"
    if (!is.null(comparison)) {
        row <- comparison[comparison$analysis == f$model, , drop = FALSE]
        stopifnot(nrow(row) == 1L, row$N == f$n, row$A == f$cases,
                  row$a_plus == f$plus_cases, row$a_minus == f$minus_cases)
        pp <- c(row$OR, row$CI_low, row$CI_high, row$raw_P)
        stopifnot(all(abs(log(c(f$or, f$lo, f$hi) / pp[1:3])) < 2e-6),
                  abs(f$raw_p - pp[4L]) < 2e-7)
        status <- "PASS"
    }
    data.frame(
        analysis = f$model, N = f$n, A = f$cases, a_plus = f$plus_cases, a_minus = f$minus_cases,
        OR_R = f$or, CI_low_R = f$lo, CI_high_R = f$hi, P_R = f$raw_p,
        OR_Python = pp[1L], CI_low_Python = pp[2L], CI_high_Python = pp[3L], P_Python = pp[4L],
        rank = f$rank, full_score = f$full_score, second_score = f$second_score,
        max_profile_score = f$profile_score, optimizer = f$optimizer, stop_code = f$optimizer_code,
        dropped_constant_nuisance = paste(names(droppable)[droppable], collapse = ","),
        R_runtime = R.version.string, method = "base R weighted QR; full-design Jeffreys penalty; profile LR",
        status = status
    )
})
result <- do.call(rbind, audit)
write.table(result, file.path(output_dir, "INDEPENDENT_FIRTH_BASE_R.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
print(result[, c("analysis", "N", "A", "OR_R", "CI_low_R", "CI_high_R", "P_R", "status")],
      row.names = FALSE)
