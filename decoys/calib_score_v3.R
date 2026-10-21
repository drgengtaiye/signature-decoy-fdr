#!/usr/bin/env Rscript
# Target-decoy calibration, step 2, per cohort: score every decoy exactly as its target's null
# sets were scored in null_matched_v3_coexpr.R (same filters, same engine, the target's own
# weight vector randomly permuted for weighted-mean signatures), then read it against the
# target's own 1,000-draw bank under N0, N1 and N2. The target is re-read from the same bank as
# a check that the stored summaries are reproduced.
#
#   Rscript calib_score_v3.R <ICB_cohort>  -> calib_v3/<cohort>_v3.csv
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment)})
COH <- commandArgs(trailingOnly=TRUE)[1]
D   <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
set.seed(91L + sum(utf8ToInt(COH)))

meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
SPEC <- meta$signature[grepl("Specific", meta$method)]
load_sig <- function(nm, key) {
  f <- file.path(D,"SignatureSets-main/data-raw/IO", paste0(nm,".csv"))
  if(!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors=FALSE, fileEncoding="latin1")
  g <- if(key=="ensembl" && "gene_id" %in% names(d)) sub("\\..*$","",d$gene_id) else d$gene_name
  w <- if("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_,nrow(d))
  ok <- !is.na(g) & g!=""
  list(genes=g[ok], weights=w[ok])
}
auc_s <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_)
  r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
engine <- function(m) if (grepl("Weighted",m)) "wmean" else if (grepl("ssGSEA",m)) "ssgsea" else "gsva"
mean_pair_r <- function(Zg) { if (nrow(Zg)<2) return(NA_real_); C <- tcrossprod(Zg)/(ncol(Zg)-1)
  (sum(C)-sum(diag(C)))/(nrow(Zg)*(nrow(Zg)-1)) }

# ---------------------------------------------------------------- cohort, as in null_matched_v3
x <- readRDS(file.path(D, paste0(COH, ".rds")))
ex <- experiments(x); an <- names(ex)
pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]
if (is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
r <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y) <- rownames(cd)
keep <- intersect(colnames(M), names(y)[!is.na(y)])
stopifnot(length(keep) >= 25)
M <- M[,keep,drop=FALSE]; yv <- y[keep]
stopifnot(sum(yv) >= 5, sum(1-yv) >= 5)
rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
key <- if (mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }
pool <- rownames(M)

wts <- list(); eng_of <- list()
for (i in seq_len(nrow(meta))) {
  s_ <- meta$signature[i]; if (s_ %in% SPEC) next
  S <- load_sig(s_, key); if (is.null(S)) next
  g <- intersect(S$genes, pool); if (length(g) < 3) next
  eng <- engine(meta$method[i]); w <- NULL
  if (eng == "wmean") {
    w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g))
    w[is.na(w)] <- 0; if (sum(abs(w)) == 0) next
  }
  eng_of[[s_]] <- eng; wts[s_] <- list(w)
}

BK  <- readRDS(file.path(D, "null_v2", paste0(COH, "_v2.rds")))
DV  <- readRDS(file.path(D, "calib_decoys_v3.rds"))
stopifnot(DV$key == key)
dec <- Filter(function(d) d$target %in% names(BK$use), DV$decoys)
cat(sprintf("[%s] n=%d decoys=%d (targets present %d)\n", COH, ncol(M), length(dec),
            length(unique(vapply(dec, `[[`, "", "target")))))

# ---------------------------------------------------------------- score the decoys
obs <- setNames(rep(NA_real_, length(dec)), names(dec))
for (eng in c("wmean","gsva","ssgsea")) {
  nm <- names(dec)[vapply(dec, function(d) eng_of[[d$target]], "") == eng]
  if (!length(nm)) next
  if (eng == "wmean") {
    for (n_ in nm) { g <- intersect(dec[[n_]]$genes, pool); ws <- wts[[dec[[n_]]$target]]
      w <- sample(ws, length(g), replace=length(ws)!=length(g))
      obs[n_] <- auc_s(yv, as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))) }
  } else {
    sets <- lapply(dec[nm], function(d) intersect(d$genes, pool))
    par  <- if (eng=="ssgsea") ssgseaParam(M, sets, minSize=3) else gsvaParam(M, sets, minSize=3)
    rr   <- gsva(par, verbose=FALSE)
    obs[rownames(rr)] <- apply(rr, 1, function(z) auc_s(yv, z))
  }
  cat(sprintf("[%s] %s: %d decoys scored\n", COH, eng, length(nm)))
}

# ---------------------------------------------------------------- read against the target's bank
rd <- function(a, b, dirn) {           # a: observed signed AUC; b: bank of signed AUCs
  bd <- if (dirn == 1L) b else 1 - b; bf <- pmax(b, 1 - b)
  c(fold_obs=max(a, 1-a), fold_med=median(bf), fold_p95=unname(quantile(bf,.95)),
    dir_obs=if (dirn==1L) a else 1-a, dir_med=median(bd), dir_p95=unname(quantile(bd,.95)))
}
rows <- list()
for (n_ in names(dec)) {
  d <- dec[[n_]]; g <- intersect(d$genes, pool)
  for (nl in c("N0","N1","N2")) {
    b <- BK$null_auc[[paste(d$target, nl)]]
    rows[[length(rows)+1]] <- data.frame(cohort=COH, n=ncol(M), unit=n_, kind="decoy", target=d$target,
      pool=d$pool, structure=d$structure, null=nl, k=length(g), dirn=d$dirn,
      r_int=mean_pair_r(Z[g,,drop=FALSE]), t(rd(obs[[n_]], b, d$dirn)))
  }
}
for (s_ in names(BK$use)) for (nl in c("N0","N1","N2")) {
  b <- BK$null_auc[[paste(s_, nl)]]
  rows[[length(rows)+1]] <- data.frame(cohort=COH, n=ncol(M), unit=s_, kind="target", target=s_,
    pool="-", structure="signature", null=nl, k=BK$use[[s_]]$k, dirn=BK$dirn[[s_]],
    r_int=mean_pair_r(Z[intersect(load_sig(s_, key)$genes, pool),,drop=FALSE]),
    t(rd(BK$obs[[s_]], b, BK$dirn[[s_]])))
}
dir.create(file.path(D,"calib_v3"), showWarnings=FALSE)
write.csv(do.call(rbind, rows), file.path(D, "calib_v3", paste0(COH, "_v3.csv")), row.names=FALSE)
cat(sprintf("[%s] DONE\n", COH))
