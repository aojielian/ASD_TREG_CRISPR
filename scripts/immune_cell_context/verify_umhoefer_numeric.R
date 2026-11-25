#!/usr/bin/env Rscript
# Independently check the R10 Umhoefer outputs with base R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
root <- args[[1]]
read <- function(file) {
    read.delim(file.path(root, "03_results", file), na.strings = "NA", check.names = FALSE)
}
u <- read("R10_TREG_TCONV_UNIVERSE.tsv.gz")
primary <- read("R10_TREG_TCONV_PRIMARY.tsv")
continuous <- read("R10_TREG_TCONV_CONTINUOUS.tsv")
adj <- read("R10_TREG_TCONV_ADJUSTED.tsv")
sensitivity <- read("R10_TREG_TCONV_SENSITIVITY.tsv")
distribution <- read("R10_TREG_TCONV_CONTINUOUS_NULL.tsv.gz")
exact_null <- read("R10_TREG_TCONV_ADJUSTED_EXACT_NULL.tsv")
as_flag <- function(x) tolower(as.character(x)) %in% c("true", "1")
d <- subset(u, category %in% c("Treg-only", "Tconv-only"))
scores <- list()
row <- function(t, name) t[t$analysis == name, ]
mark <- function(check, ours, recorded, tolerance = 1e-5) {
    err <- max(abs(ours - recorded))
    status <- if (all(abs(ours - recorded) <= tolerance * (1 + abs(recorded)))) "PASS" else "FAIL"
    scores[[length(scores) + 1L]] <<- data.frame(check, max_abs_difference = err, status)
}
matrix_2x2 <- function(asd, treg) table(factor(asd, levels = c(TRUE, FALSE)),
                                       factor(treg, levels = c(TRUE, FALSE)))
fisher_values <- function(m) {
    N <- sum(m)
    n <- sum(m[1, ])
    s <- sum(m[, 1])
    observed <- m[1, 1]
    support <- seq.int(max(0, n + s - N), min(n, s))
    base <- lchoose(s, support) + lchoose(N - s, n - support)
    prob <- function(log_odds) {
        v <- base + log_odds * support
        v <- exp(v - max(v))
        v / sum(v)
    }
    root <- function(f) uniroot(f, c(-60, 60), tol = 1e-11)$root
    estimate <- if (observed == min(support)) {
        0
    } else if (observed == max(support)) {
        Inf
    } else {
        exp(root(function(t) sum(support * prob(t)) - observed))
    }
    lower <- if (observed == min(support)) {
        0
    } else {
        exp(root(function(t) sum(prob(t)[support >= observed]) - .025))
    }
    upper <- if (observed == max(support)) {
        Inf
    } else {
        exp(root(function(t) sum(prob(t)[support <= observed]) - .025))
    }
    c(estimate, lower, upper, fisher.test(m)$p.value)
}
result_names <- c("OR", "CI_low", "CI_high", "P")
v <- as.numeric(row(primary, "discordant_primary")[result_names])
mark("direct_fisher_OR_CI_P", fisher_values(matrix_2x2(d$asd == 1, d$treg_hit == 1)), v, 1e-8)
for (side in c("treg", "tconv")) {
    v <- as.numeric(row(primary, paste0(side, "_enrichment_shared_U"))[result_names])
    mark(paste0(side, "_single_fisher_OR_CI_P"), fisher_values(matrix_2x2(
        u$asd == 1, u[[paste0(side, "_hit")]] == 1)), v, 1e-8)
}

for (i in seq_len(nrow(sensitivity))) {
    r <- sensitivity[i, ]
    cutoff <- function(x) sort(abs(x), decreasing = TRUE)[ceiling(r$target_fraction * nrow(u))]
    a <- abs(u$treg_signed_foxp3_bin_lfc) >= cutoff(u$treg_signed_foxp3_bin_lfc)
    b <- abs(u$tconv_signed_foxp3_bin_lfc) >= cutoff(u$tconv_signed_foxp3_bin_lfc)
    compare <- xor(a, b)
    mark(paste0(r$analysis, "_group_margins"),
         c(sum(a), sum(b), sum(compare), sum(u$asd[compare])), c(r$treg_n, r$tconv_n, r$n, r$asd_n))
    mark(paste0(r$analysis, "_fisher_OR_CI_P"), fisher_values(matrix_2x2(u$asd[compare] == 1, a[compare])),
         as.numeric(r[result_names]), 1e-8)
}
mark("BH_three_thresholds", p.adjust(sensitivity$P, method = "BH"), sensitivity$P_BH, 1e-12)

cc_d <- d[as_flag(d$complete_case), ]
cc_u <- u[as_flag(u$complete_case), ]
for (kind in c("abs", "signed")) {
    transform <- if (kind == "abs") abs else identity
    trank <- rank(transform(u$treg_signed_foxp3_bin_lfc), ties.method = "average")
    crank <- rank(transform(u$tconv_signed_foxp3_bin_lfc), ties.method = "average")
    delta <- (trank - crank) / (nrow(u) - 1)
    field <- paste0("delta_", kind, "_rank")
    mark(paste0(field, "_recompute"), delta, u[[field]], 1e-12)
    obs <- mean(delta[u$asd == 1]) - mean(delta[u$asd == 0])
    med <- median(delta[u$asd == 1]) - median(delta[u$asd == 0])
    r <- row(continuous, paste0(kind, "_rank_mean"))
    mark(paste0(field, "_mean"), obs, r$effect, 1e-12)
    mark(paste0(field, "_median"), med, row(continuous, paste0(kind, "_rank_median"))$effect, 1e-12)
    check_tail <- function(null, obs, reported) {
        b <- sum(abs(null) >= abs(obs) - 1e-12)
        c(b - reported$mc_extreme_n, (b + 1) / (length(null) + 1) - reported$P)
    }
    mark(paste0(field, "_null_p"), check_tail(distribution[[paste0("unadjusted_", field)]], obs, r), c(0, 0), 1e-12)
    strata <- Filter(function(g) length(unique(g$asd)) == 2L, split(cc_u, cc_u$stratum))
    weights <- vapply(strata, function(g) sum(g$asd), numeric(1))
    differences <- vapply(strata, function(g) mean(g[[field]][g$asd == 1]) - mean(g[[field]][g$asd == 0]), numeric(1))
    r <- row(adj, paste0("continuous_", field, "_conditioned"))
    obs <- weighted.mean(differences, weights)
    mark(paste0(field, "_stratum_weighted_mean"), obs, r$effect, 1e-12)
    mark(paste0(field, "_stratified_null_p"), check_tail(distribution[[paste0("adjusted_", field)]], obs, r), c(0, 0), 1e-12)
}
# Verify the finite conditional null using R's exact hypergeometric distribution
# and an explicit sum over the Cartesian product of feasible local counts.
parts <- lapply(split(cc_d, cc_d$stratum), function(g) {
    successes <- sum(g$treg_hit)
    n <- nrow(g)
    a <- sum(g$asd)
    k <- seq.int(max(0, a - (n - successes)), min(a, successes))
    cbind(k = k, p = dhyper(k, successes, n - successes, a))
})
counts <- expand.grid(lapply(parts, function(x) x[, "k"]))
mass_by_row <- rep(1, nrow(counts))
for (i in seq_along(parts)) {
    local <- parts[[i]]
    mass_by_row <- mass_by_row * local[match(counts[[i]], local[, "k"]), "p"]
}
pmf <- tapply(mass_by_row, rowSums(counts), sum)
r <- row(adj, "discordant_complete_case_conditioned_exact")
local_mean <- sum(as.numeric(names(pmf)) * pmf)
local_var <- sum((as.numeric(names(pmf)) - local_mean)^2 * pmf)
k_obs <- sum(cc_d$asd * cc_d$treg_hit)
p <- sum(pmf[abs(as.numeric(names(pmf)) - local_mean) >= abs(k_obs - local_mean) - 1e-12])
mark("conditioned_discordant_enumerated_exact_p", c(k_obs, local_mean, local_var, p),
     c(r$observed_count, r$expected_count, r$null_variance, r$P), 1e-12)
mark("conditioned_discordant_enumerated_exact_pmf", as.numeric(pmf), exact_null$null_mass, 1e-12)

# Independent penalized likelihood using numeric gradients and exact
# conditional support. These are gene-level profiles, not Wald intervals.
scale_pop <- function(x) (x - mean(x)) / sqrt(mean((x - mean(x))^2))
basic <- cbind(1, cc_d$asd)
rich <- cbind(basic, cc_d$TF, cc_d$EpiFactors, scale_pop(cc_d$expression), scale_pop(cc_d$LOEUF))
firth <- function(x, y) {
    objective <- function(beta) {
        eta <- drop(x %*% beta)
        weights <- plogis(eta) * plogis(-eta)
        det <- determinant(crossprod(x * sqrt(weights)), logarithm = TRUE)
        if (det$sign <= 0 || !is.finite(det$modulus)) return(1e99)
        -sum(y * eta - pmax(0, eta) - log1p(exp(-abs(eta)))) - as.numeric(det$modulus) / 2
    }
    start <- rep(0, ncol(x))
    start[1] <- qlogis((sum(y) + .5) / (length(y) + 1))
    control <- list(iter.max = 400, eval.max = 1500, rel.tol = 1e-12)
    out <- nlminb(start, objective, control = control)
    beta <- out$par
    profile <- function(value) {
        restricted <- function(theta) objective(append(theta, value, after = 1L))
        2 * (nlminb(beta[-2], restricted, control = control)$objective - out$objective)
    }
    criterion <- function(value) profile(value) - qchisq(.95, 1)
    left <- uniroot(criterion, c(beta[2] - 10, beta[2]), tol = 1e-9)$root
    right <- uniroot(criterion, c(beta[2], beta[2] + 10), tol = 1e-9)$root
    c(exp(beta[2]), exp(left), exp(right), pchisq(profile(0), 1, lower.tail = FALSE))
}
for (name in c("ASD_only", "adjusted")) {
    r <- row(adj, paste0("discordant_complete_case_Firth_", name))
    x <- if (name == "ASD_only") basic else rich
    mark(paste0("independent_Firth_", name), firth(x, cc_d$treg_hit), as.numeric(r[result_names]), 5e-5)
}
qa <- do.call(rbind, scores)
write.table(qa, file.path(root, "08_logs", "R10_UMHOEFER_INDEPENDENT_R_QA.tsv"), row.names = FALSE,
            quote = FALSE, sep = "\t")
print(qa, row.names = FALSE)
if (!all(qa$status == "PASS")) quit(save = "no", status = 1)
