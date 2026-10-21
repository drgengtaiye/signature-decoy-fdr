#!/usr/bin/env Rscript
# Test B of validation/PLAN_three_tests.md (SHA-256 e2edd3db3cbcb4e4…): the first principal component
# of the 52 signature scores in each held-out cohort, against the signature fixed in advance as the
# best in development. Scores use each publication's own engine, as everywhere else in this paper.
suppressMessages({library(GSVA); library(matrixStats); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
NCORE <- as.integer(Sys.getenv("NCORE", "4"))
meta <- read.csv(file.path(D, "sig_methods_clean.csv"), stringsAsFactors = FALSE)
SPEC <- meta$signature[grepl("Specific", meta$method)]
FAM <- read.csv(file.path(D, "v3_dissociation.csv"), stringsAsFactors = FALSE)
INF <- FAM$signature[FAM$family == "Inflamed/IFN"]
load_sig <- function(nm) {
  f <- file.path(D, "SignatureSets-main/data-raw/IO", paste0(nm, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "latin1")
  w <- if ("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_, nrow(d))
  ok <- !is.na(d$gene_name) & d$gene_name != ""
  list(genes = d$gene_name[ok], weights = w[ok])
}
engine <- function(m) if (grepl("Weighted", m)) "wmean" else if (grepl("ssGSEA", m)) "ssgsea" else "gsva"
auc <- function(y, v) { ok <- is.finite(v); y <- y[ok]; v <- v[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  a <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0); max(a, 1 - a) }

run <- function(nm) {
  M <- as.matrix(read.delim(gzfile(file.path(V, "cohorts", paste0(nm, ".expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  cl <- read.delim(file.path(V, "cohorts", paste0(nm, ".clin.tsv")), stringsAsFactors = FALSE)
  yv <- setNames(cl$y, cl$sample); keep <- intersect(colnames(M), names(yv)[!is.na(yv)])
  M <- M[, keep, drop = FALSE]; yv <- yv[keep]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); s[s == 0] <- 1; Z <- (M - rowMeans(M)) / s
  pool <- rownames(M)
  obs <- list(); gs <- list(); ss <- list()
  for (i in seq_len(nrow(meta))) {
    s_ <- meta$signature[i]; if (s_ %in% SPEC) next
    S <- load_sig(s_); if (is.null(S)) next
    g <- intersect(S$genes, pool); if (length(g) < 3) next
    e <- engine(meta$method[i])
    if (e == "wmean") {
      w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g))
      w[is.na(w)] <- 0; if (sum(abs(w)) == 0) next
      obs[[s_]] <- as.numeric(crossprod(Z[g, , drop = FALSE], w)) / sum(abs(w))
    } else if (e == "ssgsea") ss[[s_]] <- g else gs[[s_]] <- g
  }
  if (length(gs)) { rr <- gsva(gsvaParam(M, gs, minSize = 3), verbose = FALSE); for (n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_, ]) }
  if (length(ss)) { rr <- gsva(ssgseaParam(M, ss, minSize = 3), verbose = FALSE); for (n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_, ]) }
  X <- do.call(cbind, obs); colnames(X) <- names(obs)
  Xs <- scale(X)
  Xs <- Xs[, apply(Xs, 2, function(v) all(is.finite(v))), drop = FALSE]
  p <- prcomp(Xs, center = FALSE, scale. = FALSE)
  pc1 <- as.numeric(p$x[, 1])
  infs <- intersect(colnames(Xs), INF)
  if (length(infs) && cor(pc1, rowMeans(Xs[, infs, drop = FALSE])) < 0) pc1 <- -pc1   # sign fixed without outcome
  a <- vapply(colnames(Xs), function(s_) auc(yv, Xs[, s_]), 1)
  data.frame(cohort = nm, n = length(yv), n_resp = sum(yv), n_sig = ncol(Xs),
             var_pc1 = summary(p)$importance[2, 1], auc_pc1 = auc(yv, pc1),
             signature = colnames(Xs), auc_sig = as.numeric(a), row.names = NULL)
}
REG <- unique(read.csv(file.path(V, "validation_by_signature.csv"))$cohort)
R <- do.call(rbind, mclapply(REG, function(n) tryCatch(run(n), error = function(e) { cat("ERR", n, conditionMessage(e), "\n"); NULL }), mc.cores = NCORE))
write.csv(R, file.path(D, "tb_axis_scores.csv"), row.names = FALSE)
cat("cohorts:", length(unique(R$cohort)), "| rows:", nrow(R), "\n")
print(unique(R[, c("cohort", "n", "n_sig", "var_pc1", "auc_pc1")]), digits = 3)
