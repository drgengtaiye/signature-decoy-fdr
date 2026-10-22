#!/usr/bin/env Rscript
# Response benchmark under nulls matched on more than gene-set size and scoring engine.
#
# null_matched.R matched a signature's null on size and engine, drawing uniformly from the
# cohort's measured genes (sample(pool, k)). Two properties of a published signature were
# left unmatched. Its genes are not a random sample of the transcriptome in expression
# level, and they are co-expressed, which widens the spread of any score built from them;
# with a folded AUC that wider spread is read at its favourable tail. Three nulls are drawn
# here per signature, each adding one matching dimension:
#   N0  size + engine, uniform draw                        (the published construction)
#   N1  + the signature's expression-decile composition
#   N2  + the signature's internal mean pairwise correlation (seed-and-grow)
# Weighted-mean signatures hand their own weight vector, randomly permuted, to every null
# set, so the null carries weights as the signature does.
# AUC is stored SIGNED, so both the folded statistic of the published analysis and a
# directed test in each signature's curated direction (sig_methods.csv 'association') follow.
#
# Usage: null_matched_v2.R <ICB_cohort_name> [ndraw]
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment)})
args    <- commandArgs(trailingOnly=TRUE)
COH     <- args[1]
NDRAW   <- if (length(args) > 1) as.integer(args[2]) else 500L
D       <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
OUT     <- file.path(D, "null_v2")
POOLCAP <- 8000L     # candidate pool for the correlation matrix
TOL     <- 0.02      # band around the target correlation when growing a set
CHUNK   <- 4000L     # gene sets per GSVA call
dir.create(OUT, showWarnings=FALSE)
BK0 <- readRDS(file.path(OUT, paste0(COH, "_v2.rds")))
set.seed(42L + sum(utf8ToInt(COH)))

meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
meta$association <- trimws(meta$association)
SPEC <- meta$signature[grepl("Specific", meta$method)]
DIRN <- setNames(ifelse(meta$association == "sensitive", 1L, -1L), meta$signature)

load_sig <- function(nm, key) {
  f <- file.path(D,"SignatureSets-main/data-raw/IO", paste0(nm,".csv"))
  if(!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors=FALSE, fileEncoding="latin1")
  g <- if(key=="ensembl" && "gene_id" %in% names(d)) sub("\\..*$","",d$gene_id) else d$gene_name
  w <- if("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_,nrow(d))
  ok <- !is.na(g) & g!=""
  list(genes=g[ok], weights=w[ok])
}
# SIGNED AUC: P(score higher in responders). The published analysis folded this.
auc_s <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_)
  r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
engine <- function(m) if (grepl("Weighted",m)) "wmean" else if (grepl("ssGSEA",m)) "ssgsea" else "gsva"

# ---- grow a gene set with a target expression-decile composition and internal correlation
grow_set <- function(dec_need, rho_target, C, idx_by_dec, S0=integer(0)) {
  need <- dec_need; S <- S0; csum <- numeric(ncol(C))
  if (length(S)) for (s in S) csum <- csum + C[, s]
  k <- sum(dec_need) + length(S0)
  while (length(S) < k) {
    av <- names(need)[need > 0]
    if (!length(av)) break
    d <- if (length(av) == 1L) av else sample(av, 1L)
    cand <- setdiff(idx_by_dec[[d]], S)
    if (!length(cand)) { need[d] <- 0L; next }
    if (!length(S)) {
      pick <- cand[sample.int(length(cand), 1L)]
    } else {
      dd <- abs(csum[cand]/length(S) - rho_target)
      ok <- which(dd <= TOL)
      pick <- if (length(ok)) cand[ok[sample.int(length(ok), 1L)]] else cand[which.min(dd)]
    }
    S <- c(S, pick); csum <- csum + C[, pick]; need[d] <- need[d] - 1L
  }
  S
}
mean_pair_r <- function(C, idx) { if (length(idx) < 2) return(NA_real_)
  sub <- C[idx, idx, drop=FALSE]; (sum(sub) - sum(diag(sub))) / (length(idx)*(length(idx)-1)) }

# ---------------------------------------------------------------- load cohort (v1 filters)
f <- file.path(D, paste0(COH, ".rds"))
x <- readRDS(f)
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
stopifnot(nrow(M) >= 3000)
key <- if (mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }
pool <- rownames(M); nS <- ncol(M)

# expression deciles over the pool
mu  <- rowMeans(M)
brk <- quantile(mu, probs=seq(0,1,by=0.1), na.rm=TRUE); brk[1] <- -Inf; brk[11] <- Inf
dec <- cut(mu, breaks=brk, labels=FALSE, include.lowest=TRUE); names(dec) <- pool

# ---------------------------------------------------------------- signatures present here
use <- list(); obs_sets <- list(); wts <- list()
for (i in seq_len(nrow(meta))) {
  s_ <- meta$signature[i]; if (s_ %in% SPEC) next
  S <- load_sig(s_, key); if (is.null(S)) next
  g <- intersect(S$genes, pool); if (length(g) < 3) next
  eng <- engine(meta$method[i])
  w <- NULL
  if (eng == "wmean") {
    w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g))
    w[is.na(w)] <- 0; if (sum(abs(w)) == 0) next
  }
  use[[s_]] <- list(k=length(g), eng=eng, genes=g); wts[[s_]] <- w
}
cat(sprintf("[%s] n=%d genes=%d sigs=%d ndraw=%d\n", COH, nS, length(pool), length(use), NDRAW))

# ---------------------------------------------------------------- correlation candidate pool
sig_genes <- unique(unlist(lapply(use, `[[`, "genes")))
if (length(pool) > POOLCAP) {
  extra <- POOLCAP - length(sig_genes)
  rest  <- setdiff(pool, sig_genes)
  take  <- if (extra > 0) {
    per <- split(rest, dec[rest])
    unlist(lapply(per, function(v) sample(v, min(length(v), ceiling(extra/10)))), use.names=FALSE)
  } else character(0)
  pool_c <- unique(c(sig_genes, take))
} else pool_c <- pool
Zc <- Z[pool_c,,drop=FALSE]
C  <- tcrossprod(Zc) / (nS - 1)
dimnames(C) <- list(pool_c, pool_c)
dec_c <- dec[pool_c]
idx_by_dec <- split(seq_along(pool_c), dec_c)
cat(sprintf("[%s] correlation pool=%d  (%.0f MB)\n", COH, length(pool_c),
            as.numeric(object.size(C))/1e6))

DRAWC <- file.path(D, "v10_store", paste0(COH, "_draws.rds"))   # v10: the draws are cached, RNG state included
if (file.exists(DRAWC)) { dc <- readRDS(DRAWC); draws <- dc$draws; diag_rows <- dc$diag_rows
  assign(".Random.seed", dc$seed, envir=.GlobalEnv); cat(sprintf("[%s] draws read from cache\n", COH))
} else {
# ---------------------------------------------------------------- draw the three nulls
draws <- list(); diag_rows <- list()
for (s_ in names(use)) {
  u <- use[[s_]]; g <- u$genes; k <- u$k
  gi <- match(g, pool_c)
  rho_t <- mean_pair_r(C, gi[!is.na(gi)])
  dn <- table(factor(dec[g], levels=1:10)); dn <- dn[dn > 0]
  for (nl in c("N0","N1","N2")) {
    sets <- vector("list", NDRAW)
    for (b in seq_len(NDRAW)) {
      sets[[b]] <- if (nl == "N0") sample(pool, k)
        else if (nl == "N1") unlist(lapply(names(dn), function(d)
               sample(pool_c[idx_by_dec[[d]]], min(dn[[d]], length(idx_by_dec[[d]])))), use.names=FALSE)
        else pool_c[grow_set(dn, rho_t, C, idx_by_dec)]
    }
    names(sets) <- sprintf("%s|%s|r%04d", s_, nl, seq_len(NDRAW))
    draws[[paste(s_, nl)]] <- list(sig=s_, null=nl, eng=u$eng, sets=sets)
    ach <- vapply(sets[seq_len(min(60, NDRAW))], function(gg) {
             ii <- match(gg, pool_c); mean_pair_r(C, ii[!is.na(ii)]) }, numeric(1))
    diag_rows[[paste(s_, nl)]] <- data.frame(cohort=COH, signature=s_, null=nl, k=k,
      rho_target=rho_t, rho_achieved=mean(ach, na.rm=TRUE),
      dec_L1=mean(vapply(sets[seq_len(min(60,NDRAW))], function(gg)
        sum(abs(table(factor(dec[gg],levels=1:10)) - table(factor(dec[g],levels=1:10))))/2, numeric(1))))
  }
}
cat(sprintf("[%s] sets drawn: %d\n", COH, sum(vapply(draws, function(d) length(d$sets), integer(1)))))
dir.create(dirname(DRAWC), showWarnings=FALSE)
saveRDS(list(draws=draws, diag_rows=diag_rows, seed=.Random.seed), DRAWC)
}

# ---------------------------------------------------------------- v10: keep the scores, not only the AUCs
# Everything above is null_matched_v3_coexpr.R line for line (same seed, same draw order), so the N0
# sets are the sets of the published bank. Below, the N0 sets are scored and their per-patient scores
# kept; weighted-mean N1/N2 sets are scored only to consume the RNG in the published order. The decoys
# are then scored by calib_score_v3.R's code under its own seed. Ranks across patients are stored,
# because every AUC under any relabelling is a product of those ranks with the labels.
#   Rscript v10_store.R <ICB_cohort>   -> v10_store/<cohort>_{tgt,dcy,n0}.f32 + _meta.rds
ST <- file.path(D, "v10_store"); dir.create(ST, showWarnings=FALSE)
score_mat <- function(sets, eng, weight_src=NULL) {           # sets x patients, published scoring
  if (eng == "wmean") {
    t(vapply(seq_along(sets), function(i) {
      g <- sets[[i]]
      w <- if (is.null(weight_src)) rep(1, length(g)) else sample(weight_src, length(g), replace=length(weight_src)!=length(g))
      as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))
    }, numeric(nS)))
  } else {
    out <- matrix(NA_real_, length(sets), nS, dimnames=list(names(sets), NULL))
    for (i in seq(1, length(sets), by=CHUNK)) {
      ss <- sets[i:min(i+CHUNK-1, length(sets))]
      par <- if (eng=="ssgsea") ssgseaParam(M, ss, minSize=3) else gsvaParam(M, ss, minSize=3)
      rr <- gsva(par, verbose=FALSE); out[rownames(rr), ] <- rr
    }
    out
  }
}
n0_scores <- list(); chk <- c()
for (eng in unique(vapply(use, `[[`, character(1), "eng"))) {
  ks <- names(draws)[vapply(draws, function(d) d$eng, character(1)) == eng]
  for (kk in ks) {
    d <- draws[[kk]]
    if (d$null != "N0" && eng != "wmean") next              # GSVA/ssGSEA draw no random numbers
    S <- score_mat(d$sets, eng, weight_src=wts[[d$sig]])
    if (d$null == "N0") { n0_scores[[d$sig]] <- S
      chk <- c(chk, max(abs(apply(S, 1, function(z) auc_s(yv, z)) - BK0$null_auc[[kk]]))) }
  }
  cat(sprintf("[%s] %s N0 scored\n", COH, eng))
}
cat(sprintf("[%s] N0 bank reproduced: max |AUC difference| = %.3g\n", COH, max(chk)))
stopifnot(max(chk) < 1e-10)

# ---------------------------------------------------------------- the targets themselves
tgt <- matrix(NA_real_, length(use), nS, dimnames=list(names(use), NULL))
gs <- list(); ss <- list()
for (s_ in names(use)) { u <- use[[s_]]
  if (u$eng == "wmean") tgt[s_, ] <- as.numeric(crossprod(Z[u$genes,,drop=FALSE], wts[[s_]]))/sum(abs(wts[[s_]]))
  else if (u$eng == "ssgsea") ss[[s_]] <- u$genes else gs[[s_]] <- u$genes }
if (length(gs)) { rr <- gsva(gsvaParam(M, gs, minSize=3), verbose=FALSE); tgt[rownames(rr), ] <- rr }
if (length(ss)) { rr <- gsva(ssgseaParam(M, ss, minSize=3), verbose=FALSE); tgt[rownames(rr), ] <- rr }
d_t <- max(abs(apply(tgt, 1, function(z) auc_s(yv, z)) - BK0$obs[rownames(tgt)]))
cat(sprintf("[%s] targets reproduced: max |AUC difference| = %.3g\n", COH, d_t)); stopifnot(d_t < 1e-10)

# ---------------------------------------------------------------- the decoys, as calib_score_v3.R scores them
set.seed(91L + sum(utf8ToInt(COH)))
DV  <- readRDS(file.path(D, "calib_decoys_v3.rds")); stopifnot(DV$key == key)
dec <- Filter(function(d) d$target %in% names(use), DV$decoys)
eng_of <- lapply(use, `[[`, "eng")
dcy <- matrix(NA_real_, length(dec), nS, dimnames=list(names(dec), NULL))
for (eng in c("wmean","gsva","ssgsea")) {
  nm <- names(dec)[vapply(dec, function(d) eng_of[[d$target]], "") == eng]
  if (!length(nm)) next
  if (eng == "wmean") {
    for (n_ in nm) { g <- intersect(dec[[n_]]$genes, pool); ws <- wts[[dec[[n_]]$target]]
      w <- sample(ws, length(g), replace=length(ws)!=length(g))
      dcy[n_, ] <- as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w)) }
  } else {
    sets <- lapply(dec[nm], function(d) intersect(d$genes, pool))
    par  <- if (eng=="ssgsea") ssgseaParam(M, sets, minSize=3) else gsvaParam(M, sets, minSize=3)
    rr   <- gsva(par, verbose=FALSE); dcy[rownames(rr), ] <- rr
  }
}
CV <- read.csv(file.path(D, "calib_v3", paste0(COH, "_v3.csv")))
CV <- CV[CV$kind == "decoy" & CV$null == "N0", ]
a_ref <- ifelse(CV$dirn == 1, CV$dir_obs, 1 - CV$dir_obs); names(a_ref) <- CV$unit
d_d <- max(abs(apply(dcy, 1, function(z) auc_s(yv, z)) - a_ref[rownames(dcy)]))
cat(sprintf("[%s] decoys reproduced: max |AUC difference| = %.3g\n", COH, d_d)); stopifnot(d_d < 1e-10)

# ---------------------------------------------------------------- store ranks (float32) and metadata
rk <- function(S) t(apply(S, 1, rank))
wr <- function(S, f) { con <- file(f, "wb"); writeBin(as.vector(t(S)), con, size=4); close(con) }  # row-major
N0 <- do.call(rbind, n0_scores[names(use)])                      # (52 x 1000) x patients, target-major
wr(rk(tgt), file.path(ST, paste0(COH, "_tgt.f32")))
wr(rk(dcy), file.path(ST, paste0(COH, "_dcy.f32")))
wr(rk(N0),  file.path(ST, paste0(COH, "_n0.f32")))
gk <- intersect(pool, unique(c(unlist(lapply(use, `[[`, "genes")), unlist(lapply(DV$decoys, `[[`, "genes")))))
wr(Z[gk, , drop=FALSE], file.path(ST, paste0(COH, "_z.f32")))
saveRDS(list(cohort=COH, n=nS, y=as.integer(yv), targets=rownames(tgt), decoys=rownames(dcy),
             dirn_t=as.integer(DIRN[rownames(tgt)]), dirn_d=vapply(dec, `[[`, 0L, "dirn")[rownames(dcy)],
             target_of=vapply(dec, `[[`, "", "target")[rownames(dcy)],
             pool_d=vapply(dec, `[[`, "", "pool")[rownames(dcy)],
             structure_d=vapply(dec, `[[`, "", "structure")[rownames(dcy)],
             ndraw=NDRAW, genes=gk, tgt_genes=lapply(use, `[[`, "genes"),
             dcy_genes=lapply(dec, function(d) intersect(d$genes, pool))),
        file.path(ST, paste0(COH, "_meta.rds")))
cat(sprintf("[%s] STORED\n", COH))
