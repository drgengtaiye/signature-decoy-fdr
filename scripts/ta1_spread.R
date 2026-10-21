#!/usr/bin/env Rscript
# Test A of validation/PLAN_three_tests.md (SHA-256 e2edd3db3cbcb4e4…, 2026-09-23 21:33 UTC):
# the predictor only - the spread of the inflamed score across each cohort's patients. No outcome.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
NCORE <- as.integer(Sys.getenv("NCORE", "6"))
FAM <- read.csv(file.path(D, "v3_dissociation.csv"), stringsAsFactors = FALSE)
INF <- FAM$signature[FAM$family == "Inflamed/IFN"]
load_sig <- function(nm) {
  f <- file.path(D, "SignatureSets-main/data-raw/IO", paste0(nm, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "latin1")
  e <- if ("gene_id" %in% names(d)) sub("\\..*$", "", d$gene_id) else rep(NA_character_, nrow(d))
  ok <- !is.na(d$gene_name) & d$gene_name != ""
  data.frame(symbol = d$gene_name[ok], ensembl = e[ok], stringsAsFactors = FALSE)
}
G <- unique(do.call(rbind, lapply(INF, load_sig)))
cat("inflamed family:", length(INF), "signatures,", nrow(G), "genes\n")

spread <- function(Z) {
  ens <- mean(grepl("^ENSG", rownames(Z))) > 0.5
  key <- if (ens) G$ensembl else G$symbol
  g <- intersect(key[!is.na(key)], rownames(Z))
  if (length(g) < 10) return(c(NA, NA))
  s <- colMeans(Z[g, , drop = FALSE])
  c(sd(s), length(g))
}
dev1 <- function(nm) {
  x <- try(suppressWarnings(readRDS(file.path(D, paste0(nm, ".rds")))), silent = TRUE)
  if (inherits(x, "try-error")) return(NULL)
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  if (is.na(pick)) pick <- grep("^expr", an, value = TRUE)[1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); s[s == 0] <- 1; Z <- (M - rowMeans(M)) / s
  v <- spread(Z); data.frame(cohort = nm, set = "development", n = ncol(M), spread = v[1], n_genes = v[2])
}
val1 <- function(nm) {
  M <- as.matrix(read.delim(gzfile(file.path(V, "cohorts", paste0(nm, ".expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  cl <- read.delim(file.path(V, "cohorts", paste0(nm, ".clin.tsv")), stringsAsFactors = FALSE)
  keep <- intersect(colnames(M), cl$sample[!is.na(cl$y)]); M <- M[, keep, drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); s[s == 0] <- 1; Z <- (M - rowMeans(M)) / s
  v <- spread(Z); data.frame(cohort = nm, set = "validation", n = ncol(M), spread = v[1], n_genes = v[2])
}
NUL <- read.csv(file.path(D, "native_response_matched_v2.csv"))
DEV <- setdiff(unique(NUL$cohort), "ICB_VanDenEnde")
REG <- unique(read.csv(file.path(V, "validation_by_signature.csv"))$cohort)
R <- rbind(do.call(rbind, mclapply(DEV, function(n) tryCatch(dev1(n), error = function(e) NULL), mc.cores = NCORE)),
           do.call(rbind, mclapply(REG, function(n) tryCatch(val1(n), error = function(e) NULL), mc.cores = NCORE)))
write.csv(R, file.path(D, "ta_cohort_spread.csv"), row.names = FALSE)
print(R, digits = 3)
