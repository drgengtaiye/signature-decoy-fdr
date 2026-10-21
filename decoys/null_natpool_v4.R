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
# v4nat: N1 and N2 only, drawn from a pool in which signature genes appear at their NATURAL share.
# null_matched_v3_coexpr.R built its matched-draw pool as every signature gene plus a stratified
# sample of the rest (8,000), so signature genes - which on average run higher in responders - made
# up 18.8% of the pool against 4.3% of the transcriptome, and the N1/N2 nulls carried part of the
# signal under test (bank signed median +0.024-0.027 vs N0 -0.003). Here the pool is a
# decile-stratified sample of ALL measured genes, a signature's internal correlation is computed
# from its own genes directly rather than through the pool, and GSVA/ssGSEA sets are scored in
# chunks across signatures.
#
# Usage: null_natpool_v4.R <ICB_cohort_name> [ndraw=300]   -> null_v4nat/<cohort>_nat.rds
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment)})
args    <- commandArgs(trailingOnly=TRUE)
COH     <- args[1]
NDRAW   <- if (length(args) > 1) as.integer(args[2]) else 300L
D       <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
OUT     <- file.path(D, "null_v4nat")
POOLCAP <- 8000L     # candidate pool for the correlation matrix
TOL     <- 0.02      # band around the target correlation when growing a set
CHUNK   <- 4000L     # gene sets per GSVA call
dir.create(OUT, showWarnings=FALSE)
set.seed(4242L + sum(utf8ToInt(COH)))

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
pool_c <- if (length(pool) > POOLCAP) unlist(lapply(split(pool, dec[pool]),
            function(v) sample(v, min(length(v), ceiling(POOLCAP/10)))), use.names=FALSE) else pool
cat(sprintf("[%s] natural pool: signature-gene share %.3f (transcriptome %.3f)\n", COH,
            mean(pool_c %in% sig_genes), mean(pool %in% sig_genes)))
Zc <- Z[pool_c,,drop=FALSE]
C  <- tcrossprod(Zc) / (nS - 1)
dimnames(C) <- list(pool_c, pool_c)
dec_c <- dec[pool_c]
idx_by_dec <- split(seq_along(pool_c), dec_c)
cat(sprintf("[%s] correlation pool=%d  (%.0f MB)\n", COH, length(pool_c),
            as.numeric(object.size(C))/1e6))

# ---------------------------------------------------------------- draw the three nulls
draws <- list(); diag_rows <- list()
for (s_ in names(use)) {
  u <- use[[s_]]; g <- u$genes; k <- u$k
  Zg <- Z[g,,drop=FALSE]; rho_t <- mean_pair_r(tcrossprod(Zg)/(nS-1), seq_along(g))
  dn <- table(factor(dec[g], levels=1:10)); dn <- dn[dn > 0]
  for (nl in c("N1","N2")) {
    sets <- vector("list", NDRAW)
    for (b in seq_len(NDRAW)) {
      sets[[b]] <- if (nl == "N1") unlist(lapply(names(dn), function(d)
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

# ---------------------------------------------------------------- score everything
score_sets <- function(sets, eng, weight_src=NULL) {
  if (eng == "wmean") {
    vapply(seq_along(sets), function(i) {
      g <- sets[[i]]
      w <- if (is.null(weight_src)) rep(1, length(g)) else sample(weight_src, length(g), replace=length(weight_src)!=length(g))
      auc_s(yv, as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w)))
    }, numeric(1))
  } else {
    out <- numeric(0)
    for (i in seq(1, length(sets), by=CHUNK)) {
      ss <- sets[i:min(i+CHUNK-1, length(sets))]
      par <- if (eng=="ssgsea") ssgseaParam(M, ss, minSize=3) else gsvaParam(M, ss, minSize=3)
      rr <- gsva(par, verbose=FALSE)
      v <- apply(rr, 1, function(z) auc_s(yv, z))
      out <- c(out, v[names(ss)])
    }
    out
  }
}
null_auc <- list()
for (eng in unique(vapply(use, `[[`, character(1), "eng"))) {
  ks <- names(draws)[vapply(draws, function(d) d$eng, character(1)) == eng]
  if (!length(ks)) next
  if (eng == "wmean") { for (kk in ks) null_auc[[kk]] <- score_sets(draws[[kk]]$sets, eng, weight_src=wts[[draws[[kk]]$sig]]) }
  else { all_s <- unlist(lapply(draws[ks], `[[`, "sets"), recursive=FALSE, use.names=FALSE)
         names(all_s) <- unlist(lapply(draws[ks], function(d) names(d$sets)))
         v <- score_sets(all_s, eng)
         for (kk in ks) null_auc[[kk]] <- v[names(draws[[kk]]$sets)] }
  cat(sprintf("[%s] %s nulls scored\n", COH, eng))
}
# observed
obs <- list()
gsva_sets <- list(); ssg_sets <- list()
for (s_ in names(use)) {
  u <- use[[s_]]
  if (u$eng == "wmean") {
    w <- wts[[s_]]
    obs[[s_]] <- auc_s(yv, as.numeric(crossprod(Z[u$genes,,drop=FALSE], w))/sum(abs(w)))
  } else if (u$eng == "ssgsea") ssg_sets[[s_]] <- u$genes else gsva_sets[[s_]] <- u$genes
}
if (length(gsva_sets)) { rr <- gsva(gsvaParam(M, gsva_sets, minSize=3), verbose=FALSE)
  for (n_ in rownames(rr)) obs[[n_]] <- auc_s(yv, as.numeric(rr[n_,])) }
if (length(ssg_sets))  { rr <- gsva(ssgseaParam(M, ssg_sets, minSize=3), verbose=FALSE)
  for (n_ in rownames(rr)) obs[[n_]] <- auc_s(yv, as.numeric(rr[n_,])) }

saveRDS(list(cohort=COH, n=nS, y=yv, ndraw=NDRAW, use=lapply(use, function(u) u[c("k","eng")]),
             dirn=DIRN[names(use)], obs=unlist(obs), null_auc=null_auc,
             diag=do.call(rbind, diag_rows)),
        file.path(OUT, paste0(COH, "_nat.rds")))
cat(sprintf("[%s] DONE -> %s\n", COH, file.path(OUT, paste0(COH,"_nat.rds"))))
