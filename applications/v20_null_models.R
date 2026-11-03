# Round 15 part B3 (validation/PLAN_round15_literature_cohorts_controls.md): the null models of a type. For each of the 15
# theme pools of round 13, 6 co-expressed and 2 random-gene decoy pools grown with decoyfdr over TCGA alone from genes in
# no theme pool and no eligible model of the type (seed 20261204 + type index), each run through the round-13 pipeline
# (univariate Cox P < 0.05, then LASSO-Cox at lambda.min) fitted in TCGA alone; fold ids of decoy pool i (0-based, in the
# order grow_decoys returns them) use seed 20261205 + i. A fitted model is a null model: what the pipeline makes from
# genes of no theme. Only TCGA is read here; the null models are scored in the validation cohorts by v20_null_scores.R.
#   DFDR_LIB=<frozen install> Rscript v20_null_models.R <TYPE> [cores]
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/v20_common.R")
args <- commandArgs(TRUE); tp <- args[1]
NC <- if (length(args) > 1) as.integer(args[2]) else 8L
ti <- match(tp, TYPES) - 1L; stopifnot(!is.na(ti))
tc <- tcga_cohorts("OS", types = tp)[[tp]]; print(tc)
pools <- lapply(theme_pools(), intersect, rownames(tc$Z))
pm <- read.csv(MODELS, stringsAsFactors = FALSE); pm <- pm[pm$type == tp, ]
excl <- unique(c(unlist(theme_pools()), unlist(strsplit(pm$genes_current[pm$genes_current != ""], ";")),
                 unlist(strsplit(pm$genes_printed, ";"))))
lib <- dfdr_library(pools, engine = "mean")
dec <- grow_decoys(list(tc), lib, exclude = excl, n_coex = 6L, n_rand = 2L, n_sets = 1L, seed = 20261204L + ti,
                   cores = NC, max_bytes = 100e9)
print(dec)
Y <- Surv(tc$time, tc$event)
fit <- function(genes, seed) {
  X <- t(tc$Z[genes, , drop = FALSE])
  pu <- vapply(seq_along(genes), function(j) {
    f <- try(coxph(Y ~ X[, j]), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else summary(f)$coefficients[1, 5] }, 1.0)
  keep <- genes[which(pu < 0.05)]
  out <- list(n_pool = length(genes), n_uni = length(keep), n_model = 0L, genes = "", coefs = "")
  if (length(keep) < 2) return(out)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  cv <- try(cv.glmnet(X[, keep, drop = FALSE], Y, family = "cox", alpha = 1, foldid = foldid), silent = TRUE)
  if (inherits(cv, "try-error")) return(out)
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]; b <- b[b != 0]
  out$n_model <- length(b)
  if (length(b) < 2) return(out)
  out$genes <- paste(names(b), collapse = ";"); out$coefs <- paste(signif(b, 8), collapse = ";")
  out
}
info <- dec$info
r <- parallel::mclapply(seq_along(dec$decoys), function(i) fit(dec$decoys[[i]]$genes, 20261205L + i - 1L),
                        mc.cores = NC, mc.preschedule = FALSE)
R <- data.frame(type = tp, pool_index = seq_along(dec$decoys) - 1L, theme_pool = info$target, class = info$class,
                seed = 20261205L + seq_along(dec$decoys) - 1L, do.call(rbind, lapply(r, as.data.frame)))
write.csv(R, gzfile(file.path(OUT, sprintf("null_models_%s.csv.gz", tp))), row.names = FALSE)
cat("done", tp, nrow(R), "decoy pools;", sum(R$n_model >= 2), "null models\n")
