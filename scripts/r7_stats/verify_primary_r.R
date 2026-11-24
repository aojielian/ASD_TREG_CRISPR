#!/usr/bin/env Rscript
# Independent base-R implementation of the same full-design Jeffreys objective.
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 2L)
data_dir <- a[[1]]
o <- a[[2]]
rd <- function(s) read.delim(s, check.names = FALSE, na.strings = c("", "NA"),
                             quote = "", comment.char = "", stringsAsFactors = FALSE)
d <- rd(file.path(data_dir, "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz"))
q <- rd(file.path(o, "results/LOEUF_ADJUSTED_MODELS.tsv"))
q <- q[q$model=="M2" & q$threshold==1 & q$units=="CPM", ]
encode <- function(z) tolower(as.character(z)) %in% c("true", "1")
scale_pop <- function(z) (z-mean(z))/sqrt(mean((z-mean(z))^2))
audit <- list()
for (i in seq_len(nrow(q))) {
  row <- q[i, ]
  v <- d[d$screen == row$screen & encode(d$measured) &
         d$primary_status_1 %in% c("DETECTABLE_NONHIT", "DETECTABLE_HIT"), ]
  ok <- complete.cases(v[, c("hgnc_id","class_label","TF_LISTED","EPIFACTORS_GENES_LISTED",
      "primary_mean_log_expression","LOEUF_v2_1_1","GENE_SPAN_BP","num")])
  v <- v[ok, ]
  stopifnot(nrow(v)==row$n_genes, sum(encode(v$sfari_score1)) == row$asd_n)
  x <- cbind(k=1, ASD=encode(v$sfari_score1),TF=v$TF_LISTED, EPI=v$EPIFACTORS_GENES_LISTED)
  if (row$screen == "CHEN") x <- cbind(x,library=encode(v$library_member))
  x <- cbind(x,expression=scale_pop(v$primary_mean_log_expression), LOEUF=scale_pop(v$LOEUF_v2_1_1))
  y <- as.numeric(encode(v$author_hit))
  neg <- function(b) {
    eta <- drop(x %*% b)
    weight <- plogis(eta) * plogis(-eta)
    I <- crossprod(x, weight * x)
    determinant_I <- determinant(I, logarithm = TRUE)
    if (determinant_I$sign != 1) return(1e20)
    -sum(y * eta - pmax(eta, 0) - log1p(exp(-abs(eta)))) -
      .5 * as.numeric(determinant_I$modulus)
  }
  init <- rep(0, ncol(x))
  init[1] <- qlogis((sum(y)+.5)/(length(y)+1))
  full <- nlminb(init, neg, control=list(iter.max=500, eval.max=2000, rel.tol=1e-13))
  score_check <- function(b, objective) sapply(seq_along(b), function(j) {
    e <- rep(0, length(b)); e[j] <- 1e-5
    (objective(b+e)-objective(b-e))/(2e-5)
  })
  stopifnot(max(abs(score_check(full$par, neg))) < 0.001)
  pr <- function(b) {
    objective <- function(z) neg(append(z,b,after=1))
    result <- nlminb(full$par[-2], objective,
                     control=list(iter.max=500,eval.max=2000,rel.tol=1e-13))
    refine <- optim(result$par, objective, method="BFGS",
                    control=list(maxit=500, ndeps=rep(1e-5,length(result$par)), reltol=1e-14))
    z <- if (refine$value < result$objective) refine$par else result$par
    res <- max(abs(score_check(z,objective)))
    if (res >= .001) stop(sprintf("%s profile beta %.6g score %.6g",row$screen,b,res))
    2 * (objective(z) - full$objective)
  }
  ci_lrt <- c(pr(log(row$CI_low)),pr(log(row$CI_high)))
  target <- qchisq(.95,1)
  P <- pchisq(pr(0),1,lower.tail=FALSE)
  independently_fit_OR <- exp(full$par[2])
  stopifnot(abs(log(independently_fit_OR/row$OR)) < .0001, max(abs(ci_lrt-target)) < .0001,
            abs(P-row$P) < .00001)
  audit[[i]] <- data.frame(screen=row$screen,R_base_OR=independently_fit_OR,Python_OR=row$OR,
        R_base_P=P,Python_P=row$P,lower_endpoint_LR=ci_lrt[1],upper_endpoint_LR=ci_lrt[2],
        chi_square_limit=target,base_R_optimizer_code=full$convergence,
        base_R_optimizer_message=full$message,status="PASS")
}
write.table(do.call(rbind,audit), file.path(o,"checks/INDEPENDENT_R_FIRTH_QA.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE)
cat("Independently verified both primary M2 ORs, penalized LR P values and interval endpoints in base R.\n")
