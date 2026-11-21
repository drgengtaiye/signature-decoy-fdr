# Round 15 part C (validation/PLAN_round15_literature_cohorts_controls.md C1; amendment 3): in each cohort of the
# independent set, the first K principal components of its z-scored genes, survival regressed on them jointly (Cox,
# gamma_c), and for each null model scorable there its standardized score regressed on the components (a_mc), its
# predicted log hazard ratio per SD beta_mc = a_mc . gamma_c (null models: higher score, shorter survival), and the
# predicted probability that it validates under (L), p_mc = Phi(0.80 sqrt(D_c) beta_mc - 1.96). K = 10 (primary), 5, 20.
# -> icb/v20_round15/law_predictions.csv.gz (one row per cohort, null model and K)
#   DFDR_LIB=<frozen install> Rscript v20_law.R [cores]
source("/path/to/analysis_root/signature-decoy-fdr/applications/v20_common.R")
args <- commandArgs(TRUE); NC <- if (length(args)) as.integer(args[1]) else 8L
ONLY <- if (length(args) > 1) args[2] else NA      # one cohort, run alone (a forked worker that dies returns nothing)
CO <- read.csv(file.path(G, "signature-decoy-fdr", "validation", "round15_cohorts", "built_cohorts_r15.csv"),
               stringsAsFactors = FALSE)
CO <- CO[CO$built & CO$independent_set, ]
if (!is.na(ONLY)) CO <- CO[CO$cohort == ONLY, ]
load_cohort <- function(id, dir) {
  e <- as.matrix(read.delim(gzfile(file.path(dir, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  cl <- read.delim(file.path(dir, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample), name = id)
}
res <- parallel::mclapply(seq_len(nrow(CO)), function(j) {
  tp <- CO$type[j]; co <- load_cohort(CO$cohort[j], CO$dir[j])
  nm <- read.csv(gzfile(file.path(OUT, sprintf("null_models_%s.csv.gz", tp))), stringsAsFactors = FALSE)
  nm <- nm[nm$n_model >= 2, ]
  pc <- prcomp(t(co$Z), center = TRUE, scale. = FALSE, rank. = 20)$x
  pc <- scale(pc)
  D <- sum(co$event); Y <- Surv(co$time, co$event)
  out <- list()
  for (K in c(10L, 5L, 20L)) {
    P <- pc[, 1:K, drop = FALSE]
    gam <- coef(coxph(Y ~ P))
    for (i in seq_len(nrow(nm))) {
      g <- strsplit(nm$genes[i], ";")[[1]]; w <- as.numeric(strsplit(nm$coefs[i], ";")[[1]])
      k <- g %in% rownames(co$Z)
      if (mean(k) < 0.8 || sum(k) < 2) next
      s <- lin_score(co, g[k], w[k]); s <- (s - mean(s)) / sd(s)
      a <- coef(lm(s ~ P))[-1]
      beta <- sum(a * gam)
      out[[length(out) + 1]] <- data.frame(type = tp, cohort = CO$cohort[j], id = sprintf("%s_null%03d", tp, nm$pool_index[i]),
                                           class = nm$class[i], K = K, events = D, beta = beta,
                                           p_pred = pnorm(0.80 * sqrt(D) * beta - 1.96))
    }
  }
  do.call(rbind, out)
}, mc.cores = NC, mc.preschedule = FALSE)
bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad)) cat("failed cohorts (rerun alone):", paste(CO$cohort[bad], collapse = " "), "\n")
L <- do.call(rbind, res[!bad])
f <- if (is.na(ONLY)) "law_predictions.csv.gz" else sprintf("law_predictions_%s.csv.gz", ONLY)
write.csv(L, gzfile(file.path(OUT, f)), row.names = FALSE)
cat("done", nrow(L), "rows over", length(unique(L$cohort)), "cohorts\n")
