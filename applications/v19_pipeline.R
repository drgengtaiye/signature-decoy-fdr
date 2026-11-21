# Round 13 part B (validation/PLAN_round13_published_risk_models.md, B1-B3): the published pipeline run on each theme
# pool and on 20 co-expressed and 5 random-gene decoy pools of it, in TCGA (training) and the external cohort of the type.
#   DFDR_LIB=<frozen install> Rscript v19_pipeline.R <TYPE> [cores] [targets|all]
# Fold ids: seed 20261117 + unit index, the unit index of theme pool p (0-based) in type t (plan order, 0-based) being
# 100 t + p; decoy j (1-20 co-expressed, 21-25 random) of that pool uses 20261117 + 100000 + 10000 t + 100 p + j.
source("/path/to/analysis_root/signature-decoy-fdr/applications/v19_common.R")
args <- commandArgs(TRUE); tp <- args[1]
NC <- if (length(args) > 1) as.integer(args[2]) else 8L
stage <- if (length(args) > 2) args[3] else "all"
ti <- match(tp, TYPES) - 1L; stopifnot(!is.na(ti))
t0 <- Sys.time(); mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
P <- load_pair(tp); print(P$tcga); print(P$ext)
univ <- intersect(rownames(P$tcga$Z), rownames(P$ext$Z))
pools <- lapply(theme_pools(), intersect, univ)
cat("pool sizes:", paste(names(pools), vapply(pools, length, 1L), sep = "=", collapse = " "), "\n")
Y <- Surv(P$tcga$time, P$tcga$event)

pipe <- function(genes, seed) {
  X <- t(P$tcga$Z[genes, , drop = FALSE])
  pu <- vapply(seq_along(genes), function(j) {
    f <- try(coxph(Y ~ X[, j]), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else summary(f)$coefficients[1, 5] }, 1.0)
  keep <- genes[which(pu < 0.05)]
  out <- list(n_pool = length(genes), n_uni = length(keep), n_model = 0L, genes = "", coefs = "")
  na <- c(tcga_cdir = NA, ext_cdir = NA, tcga_p = NA, tcga_worse = NA, ext_p = NA, ext_worse = NA, validated = 0)
  if (length(keep) < 2) return(c(out, as.list(na)))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  cv <- try(cv.glmnet(X[, keep, drop = FALSE], Y, family = "cox", alpha = 1, foldid = foldid), silent = TRUE)
  if (inherits(cv, "try-error")) return(c(out, as.list(na)))
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]; b <- b[b != 0]
  out$n_model <- length(b)
  if (length(b) < 2) return(c(out, as.list(na)))
  out$genes <- paste(names(b), collapse = ";"); out$coefs <- paste(signif(b, 8), collapse = ";")
  st <- lin_score(P$tcga, names(b), b); se <- lin_score(P$ext, names(b), b)
  c(out, as.list(score_stats(st, se, P, dirn = 1)))
}

run_units <- function(U) {
  r <- parallel::mclapply(seq_len(nrow(U)), function(i) pipe(U$genes[[i]], U$seed[i]), mc.cores = NC,
                          mc.preschedule = FALSE)
  bad <- vapply(r, function(x) is.null(x) || inherits(x, "try-error"), TRUE)
  if (any(bad)) stop("pipeline failed for units ", paste(which(bad), collapse = ","))
  cbind(U[, c("type", "pool", "kind", "class", "j", "seed")], do.call(rbind, lapply(r, as.data.frame)))
}

U <- data.frame(type = tp, pool = names(pools), kind = "target", class = NA_character_, j = 0L,
                seed = 20261117L + 100L * ti + seq_along(pools) - 1L)
U$genes <- unname(pools)
RT <- run_units(U)
write.csv(RT, gzfile(file.path(OUT, sprintf("pipeline_targets_%s.csv.gz", tp))), row.names = FALSE)
cat("targets", mins(), "\n"); print(RT[, c("pool", "n_pool", "n_uni", "n_model", "validated")])

if (stage == "all") {
  pm <- published_models(tp); stopifnot(!is.null(pm))   # the frozen table must exist before decoys are grown
  excl <- unique(unlist(strsplit(pm$genes_current, ";")))
  lib <- dfdr_library(pools, engine = "mean")
  dec <- grow_decoys(list(P$tcga, P$ext), lib, exclude = excl, n_coex = 20L, n_rand = 5L, n_sets = 1L,
                     seed = 20261118L + ti, cores = NC, max_bytes = 100e9)
  print(dec)
  info <- dec$info
  info$j <- ave(seq_len(nrow(info)), info$target, FUN = seq_along)   # coex 1-20 then rand 21-25 within a pool
  write.csv(info, gzfile(file.path(OUT, sprintf("pipeline_decoy_info_%s.csv.gz", tp))), row.names = FALSE)
  pidx <- match(info$target, names(pools)) - 1L
  D <- data.frame(type = tp, pool = info$target, kind = "decoy", class = info$class, j = info$j,
                  seed = 20261117L + 100000L + 10000L * ti + 100L * pidx + info$j)
  D$genes <- lapply(dec$decoys, `[[`, "genes")
  RD <- run_units(D)
  write.csv(RD, gzfile(file.path(OUT, sprintf("pipeline_decoys_%s.csv.gz", tp))), row.names = FALSE)
  cat("decoys", mins(), "\n")
}
