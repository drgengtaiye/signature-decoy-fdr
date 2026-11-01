# Round 11 (validation/PLAN_round11_external_replication.md): the 52 immune signatures in independent cohorts, with
# decoys grown by the package exactly as for TCGA (ten sets of 32 co-expressed and 8 random-gene decoys; 1,000-draw
# exact-size banks per target and cohort).
#   DFDR_LIB=<frozen install> Rscript v16_external.R survival|response [cores]
#   survival: the external survival cohorts of external/cohorts (one endpoint each), decoys from the genes measured
#             in all of them and in no curated signature
#   response: GSE176307 (metastatic urothelial carcinoma, checkpoint blockade; validation/cohorts), alone
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
args <- commandArgs(TRUE); app <- args[1]; NC <- if (length(args) > 1) as.integer(args[2]) else 8L
OUT <- file.path(G, "icb", "v16_external"); dir.create(OUT, showWarnings = FALSE)
K <- 10L; t0 <- Sys.time()
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
L <- immune_library(key = "symbol")
EXT <- file.path(G, "external", "cohorts")

if (app == "survival") {
  man <- read.csv(file.path(EXT, "manifest.csv"), stringsAsFactors = FALSE)
  man <- man[as.logical(man$included), ]
  co <- lapply(man$cohort, function(id) {
    e <- read.delim(gzfile(file.path(EXT, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE)
    M <- as.matrix(e); M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
    cl <- read.delim(file.path(EXT, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
    dfdr_cohort(M[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample),
                name = id)
  })
  # AMENDMENT_1_round11: sets with fewer than three genes measured in some cohort are left out (as for the breast library)
  lb <- L$lib
  nmeas <- sapply(lb$sets, function(g) min(vapply(co, function(c) sum(g %in% rownames(c$Z)), 1L)))
  keep <- names(nmeas)[nmeas >= 3]; cat("left out:", setdiff(names(nmeas), keep), "\n")
  wm <- keep[lb$engine[keep] == "wmean"]
  L$lib <- dfdr_library(lb$sets[keep], lb$engine[keep], lb$weights[wm], lb$direction[keep]); print(L$lib)
  dec <- grow_decoys(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, n_sets = K, seed = 20261101L,
                     cores = NC, max_bytes = 200e9)
  print(dec); saveRDS(dec, file.path(OUT, "external_decoys.rds"))
  write.csv(dec$info, gzfile(file.path(OUT, "external_decoy_info.csv.gz")), row.names = FALSE)
  cat("decoys", mins(), "\n")
  fit <- fit_units(co, L$lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261102L, cores = NC)
  write.csv(fit$units, gzfile(file.path(OUT, "external_surv_units.csv.gz")), row.names = FALSE)
  calls <- call_units(fit$units)
  write.csv(calls, gzfile(file.path(OUT, "external_surv_calls.csv.gz")), row.names = FALSE)
  for (st in c("fold")) {
    f <- decoy_fdr(calls, st, "coex", B = 6000L)
    write.csv(f, file.path(OUT, sprintf("external_surv_%s_fdr.csv", st)), row.names = FALSE); print(f)
  }
  cat("fit", mins(), "\n")
}

if (app == "response") {
  V <- file.path(G, "validation", "cohorts")
  e <- read.delim(gzfile(file.path(V, "GSE176307.expr.tsv.gz")), row.names = 1, check.names = FALSE)
  M <- as.matrix(e); M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  cl <- read.delim(file.path(V, "GSE176307.clin.tsv"), stringsAsFactors = FALSE)
  y <- setNames(cl$y, cl$sample)
  co <- list(dfdr_cohort(M[, cl$sample, drop = FALSE], y, name = "GSE176307"))
  lb <- L$lib
  nmeas <- sapply(lb$sets, function(g) sum(g %in% rownames(co[[1]]$Z)))
  keep <- names(nmeas)[nmeas >= 3]; cat("left out:", setdiff(names(nmeas), keep), "\n")
  wm <- keep[lb$engine[keep] == "wmean"]
  L$lib <- dfdr_library(lb$sets[keep], lb$engine[keep], lb$weights[wm], lb$direction[keep]); print(L$lib)
  dec <- grow_decoys(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, n_sets = K, seed = 20261103L,
                     cores = NC, max_bytes = 200e9)
  print(dec); write.csv(dec$info, gzfile(file.path(OUT, "external_resp_decoy_info.csv.gz")), row.names = FALSE)
  fit <- fit_units(co, L$lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261104L, cores = NC)
  write.csv(fit$units, gzfile(file.path(OUT, "external_resp_units.csv.gz")), row.names = FALSE)
  cat("fit", mins(), "\n")
}
