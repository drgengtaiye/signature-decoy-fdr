# Round 18 (validation/PLAN_round18_law_out_of_sample.md): gamma fitted on half of a cohort's patients; the rate of validation predicted and observed in the other half.
# -> icb/v23_law_split/law_split[_<cohort>].csv.gz, one row per cohort, split, K, null model
#   DFDR_LIB=<frozen install> Rscript v23_law_split.R [cores] [cohort]
source("/path/to/analysis_root/signature-decoy-fdr/applications/v20_common.R")
args <- commandArgs(TRUE); NC <- if (length(args)) as.integer(args[1]) else 8L
ONLY <- if (length(args) > 1) args[2] else NA
O <- file.path(G, "icb", "v23_law_split"); dir.create(O, showWarnings = FALSE)
CO <- read.csv(file.path(G, "signature-decoy-fdr", "validation", "round15_cohorts", "built_cohorts_r15.csv"), stringsAsFactors = FALSE)
CO <- CO[CO$built & CO$independent_set, ]; CO$index <- seq_len(nrow(CO)) - 1L
if (!is.na(ONLY)) CO <- CO[CO$cohort == ONLY, ]
load_cohort <- function(id, dir) {
  e <- as.matrix(read.delim(gzfile(file.path(dir, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  cl <- read.delim(file.path(dir, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample), name = id)
}
SPLITS <- 20L
res <- parallel::mclapply(seq_len(nrow(CO)), function(j) {
  tp <- CO$type[j]; co <- load_cohort(CO$cohort[j], CO$dir[j])
  nm <- read.csv(gzfile(file.path(OUT, sprintf("null_models_%s.csv.gz", tp))), stringsAsFactors = FALSE)
  nm <- nm[nm$n_model >= 2, ]
  pc <- scale(prcomp(t(co$Z), center = TRUE, scale. = FALSE, rank. = 20)$x)
  time <- co$time; event <- co$event; n <- length(time)
  # scores and loadings on all patients, without outcome
  S <- list(); LD <- list()
  for (i in seq_len(nrow(nm))) {
    g <- strsplit(nm$genes[i], ";")[[1]]; w <- as.numeric(strsplit(nm$coefs[i], ";")[[1]]); k <- g %in% rownames(co$Z)
    if (mean(k) < 0.8 || sum(k) < 2) next
    s <- lin_score(co, g[k], w[k]); s <- (s - mean(s)) / sd(s)
    id <- sprintf("%s_null%03d", tp, nm$pool_index[i]); S[[id]] <- s
    LD[[id]] <- lapply(c(`10` = 10L, `5` = 5L, `20` = 20L), function(K) coef(lm(s ~ pc[, 1:K, drop = FALSE]))[-1])
  }
  if (!length(S)) return(NULL)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261301L + CO$index[j])
  out <- list(); dropped <- 0L
  for (sp in seq_len(SPLITS)) {
    A <- logical(n)
    for (ev in 0:1) { ix <- which(event == ev); A[sample(ix, floor(length(ix) / 2))] <- TRUE }   # stratified by event
    B <- !A; DB <- sum(event[B])
    Y <- Surv(time[A], event[A])
    for (K in c(10L, 5L, 20L)) {
      if (DB < K + 1 || sum(event[A]) < K + 1) { dropped <- dropped + 1L; next }
      P <- pc[A, 1:K, drop = FALSE]
      fit <- tryCatch(suppressWarnings(coxph(Y ~ P)), error = function(e) NULL)
      gam <- if (is.null(fit)) NULL else coef(fit)
      if (is.null(gam) || anyNA(gam) || any(abs(gam) > 20)) { dropped <- dropped + 1L; next }
      for (id in names(S)) {
        beta <- sum(LD[[id]][[as.character(K)]] * gam)
        kv <- km_valid(S[[id]][B], time[B], event[B])
        out[[length(out) + 1]] <- data.frame(type = tp, cohort = CO$cohort[j], split = sp, K = K, id = id, D_B = DB, beta = beta,
                                             p_pred = pnorm(0.80 * sqrt(DB) * beta - 1.96),
                                             validated = as.numeric(isTRUE(kv[["p"]] < 0.05) && isTRUE(kv[["worse"]] == 1)))
      }
    }
  }
  r <- do.call(rbind, out); attr(r, "dropped") <- dropped; cat(CO$cohort[j], "rows", nrow(r), "dropped K-splits", dropped, "\n"); r
}, mc.cores = NC, mc.preschedule = FALSE)
bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad)) cat("failed cohorts (rerun alone):", paste(CO$cohort[bad], collapse = " "), "\n")
L <- do.call(rbind, res[!bad])
f <- if (is.na(ONLY)) "law_split.csv.gz" else sprintf("law_split_%s.csv.gz", ONLY)
write.csv(L, gzfile(file.path(O, f)), row.names = FALSE)
cat("done", nrow(L), "rows over", length(unique(L$cohort)), "cohorts\n")
