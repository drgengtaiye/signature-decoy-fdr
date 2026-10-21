#!/usr/bin/env Rscript
# Step 1 of validation/PLAN_closed_form_null.md (SHA-256 edd47e92343ea3bc…, 2026-09-23 21:08 UTC).
# Two cheap statistics per cohort: the spread of single-gene folded AUC, and the mean absolute
# correlation between random gene pairs. One pass over each expression matrix.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
NCORE <- as.integer(Sys.getenv("NCORE", "6"))
NUL <- read.csv(file.path(D, "native_response_matched_v2.csv"))
COH <- sort(unique(NUL$cohort))
auc <- function(y, v) { ok <- is.finite(v); y <- y[ok]; v <- v[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  a <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0); max(a, 1 - a) }

run <- function(nm) {
  x <- try(suppressWarnings(readRDS(file.path(D, paste0(nm, ".rds")))), silent = TRUE)
  if (inherits(x, "try-error")) return(NULL)
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  if (is.na(pick)) pick <- grep("^expr", an, value = TRUE)[1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]; yv <- y[keep]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); s[s == 0] <- 1; Z <- (M - rowMeans(M)) / s
  set.seed(sum(utf8ToInt(nm)))
  g1 <- sample(rownames(Z), min(2000, nrow(Z)))
  a1 <- vapply(g1, function(g) auc(yv, Z[g, ]), 1)
  g2 <- sample(rownames(Z), min(600, nrow(Z)))
  C <- suppressWarnings(cor(t(Z[g2, , drop = FALSE])))
  rho <- mean(abs(C[upper.tri(C)]), na.rm = TRUE)
  data.frame(cohort = nm, n = length(yv), n_resp = sum(yv), n_genes = nrow(Z),
             sigma1 = sd(a1, na.rm = TRUE), auc1_med = median(a1, na.rm = TRUE), rho = rho)
}
R <- do.call(rbind, mclapply(COH, function(nm) tryCatch(run(nm), error = function(e) NULL), mc.cores = NCORE))
write.csv(R, file.path(D, "cf_cohort_stats.csv"), row.names = FALSE)
print(R, digits = 3)
