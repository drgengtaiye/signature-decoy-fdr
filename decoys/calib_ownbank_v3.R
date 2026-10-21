#!/usr/bin/env Rscript
# Target-decoy calibration v3, check: does reading a decoy against its TARGET's null bank give the
# same answer as building the decoy its OWN bank by the published procedure? For N0 the bank has
# no structure and the two must agree; under N1/N2 the bank is not centred at 0.5 (co-expressed
# null sets sit on real biological axes), so the shared-bank design has to be checked, not assumed.
# The first PER_T pool-A coex decoys of every target get their own N0/N1/N2 banks here, drawn and
# scored exactly as in null_matched_v3_coexpr.R; the decoy's observed AUC is taken from
# calib_v3/<cohort>_v3.csv so that only the reference changes.
#
#   Rscript calib_ownbank_v3.R <ICB_cohort> [nbank=300] [per_target=4] -> calib_v3/<cohort>_own.csv
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment)})
a <- commandArgs(trailingOnly=TRUE); COH <- a[1]
NB    <- if (length(a) > 1) as.integer(a[2]) else 300L
PER_T <- if (length(a) > 2) as.integer(a[3]) else 4L
PMODE <- if (length(a) > 3) a[4] else "enriched"   # "natural": signature genes at their natural share, N1/N2 only
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
POOLCAP <- 8000L; TOL <- 0.02; CHUNK <- 4000L
NW <- as.integer(Sys.getenv("GSVA_WORKERS", "1"))   # GSVA is single-threaded unless given workers
set.seed(313L + sum(utf8ToInt(COH)))

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
grow_set <- function(dec_need, rho_target, C, idx_by_dec) {
  need <- dec_need; S <- integer(0); csum <- numeric(ncol(C)); k <- sum(dec_need)
  while (length(S) < k) {
    av <- names(need)[need > 0]; if (!length(av)) break
    d <- if (length(av)==1L) av else sample(av,1L)
    cand <- setdiff(idx_by_dec[[d]], S)
    if (!length(cand)) { need[d] <- 0L; next }
    pick <- if (!length(S)) cand[sample.int(length(cand),1L)] else {
      dd <- abs(csum[cand]/length(S) - rho_target); ok <- which(dd <= TOL)
      if (length(ok)) cand[ok[sample.int(length(ok),1L)]] else cand[which.min(dd)] }
    S <- c(S, pick); csum <- csum + C[, pick]; need[d] <- need[d] - 1L
  }
  S
}
mean_pair_r <- function(C, idx) { if (length(idx)<2) return(NA_real_)
  sub <- C[idx,idx,drop=FALSE]; (sum(sub)-sum(diag(sub)))/(length(idx)*(length(idx)-1)) }

x <- readRDS(file.path(D, paste0(COH, ".rds")))
ex <- experiments(x); an <- names(ex)
pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]
if (is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
r <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y) <- rownames(cd)
keep <- intersect(colnames(M), names(y)[!is.na(y)])
M <- M[,keep,drop=FALSE]; yv <- y[keep]
rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
key <- if (mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }
pool <- rownames(M); nS <- ncol(M)
mu  <- rowMeans(M)
brk <- quantile(mu, probs=seq(0,1,by=0.1), na.rm=TRUE); brk[1] <- -Inf; brk[11] <- Inf
dec <- cut(mu, breaks=brk, labels=FALSE, include.lowest=TRUE); names(dec) <- pool

wts <- list(); eng_of <- list(); sig_genes <- character(0)
for (i in seq_len(nrow(meta))) {
  s_ <- meta$signature[i]; if (s_ %in% SPEC) next
  S <- load_sig(s_, key); if (is.null(S)) next
  g <- intersect(S$genes, pool); if (length(g) < 3) next
  eng <- engine(meta$method[i]); w <- NULL
  if (eng == "wmean") {
    w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g))
    w[is.na(w)] <- 0; if (sum(abs(w)) == 0) next
  }
  eng_of[[s_]] <- eng; wts[s_] <- list(w); sig_genes <- union(sig_genes, g)
}

V  <- read.csv(file.path(D, "calib_v3", paste0(COH, "_v3.csv")), stringsAsFactors=FALSE)
V  <- V[V$kind=="decoy" & V$null=="N0" & V$pool=="A" & V$structure=="coex",]
V$b <- as.integer(sub(".*_", "", V$unit))
V  <- V[V$b <= PER_T,]
DV <- readRDS(file.path(D, "calib_decoys_v3.rds"))$decoys
sub_ <- DV[V$unit]
cat(sprintf("[%s] n=%d decoys checked=%d nbank=%d\n", COH, nS, length(sub_), NB))

dg <- unique(unlist(lapply(sub_, `[[`, "genes")))
extra <- POOLCAP - length(union(dg, sig_genes)); rest <- setdiff(pool, union(dg, sig_genes))
take <- if (extra > 0) unlist(lapply(split(rest, dec[rest]),
          function(v) sample(v, min(length(v), ceiling(extra/10)))), use.names=FALSE) else character(0)
pool_c <- if (PMODE == "natural") unlist(lapply(split(pool, dec[pool]),
            function(v) sample(v, min(length(v), ceiling(POOLCAP/10)))), use.names=FALSE) else unique(c(sig_genes, dg, take))
cat(sprintf("[%s] pool mode %s: signature-gene share %.3f\n", COH, PMODE, mean(pool_c %in% sig_genes)))
C <- tcrossprod(Z[pool_c,,drop=FALSE])/(nS-1); dimnames(C) <- list(pool_c, pool_c)
idx_by_dec <- split(seq_along(pool_c), dec[pool_c])

score_sets <- function(sets, eng, weight_src=NULL) {
  if (eng == "wmean") {
    vapply(seq_along(sets), function(i) { g <- sets[[i]]
      w <- if (is.null(weight_src)) rep(1, length(g)) else sample(weight_src, length(g), replace=length(weight_src)!=length(g))
      auc_s(yv, as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))) }, numeric(1))
  } else {
    out <- numeric(0)
    for (i in seq(1, length(sets), by=CHUNK)) {
      ss <- sets[i:min(i+CHUNK-1, length(sets))]
      par <- if (eng=="ssgsea") ssgseaParam(M, ss, minSize=3) else gsvaParam(M, ss, minSize=3)
      bp <- if (NW > 1) BiocParallel::MulticoreParam(NW) else BiocParallel::SerialParam()
      rr <- gsva(par, verbose=FALSE, BPPARAM=bp); v <- apply(rr, 1, function(z) auc_s(yv, z))
      out <- c(out, v[names(ss)])
    }
    out
  }
}
# draw every bank first, then score each engine in chunks, as null_matched_v3_coexpr.R does
# (one GSVA call per bank re-estimates the gene-wise kernels every time; batching does not change
# any set's score, since GSVA scores each set independently of the others in the call)
meta_b <- list(); banks <- list()
for (i in seq_along(sub_)) {
  d <- sub_[[i]]; g <- intersect(d$genes, pool); k <- length(g)
  rho_t <- { Zg <- Z[g,,drop=FALSE]; Cg <- tcrossprod(Zg)/(nS-1); (sum(Cg)-sum(diag(Cg)))/(k*(k-1)) }
  dn <- table(factor(dec[g], levels=1:10)); dn <- dn[dn > 0]
  for (nl in (if (PMODE == "natural") c("N1","N2") else c("N0","N1","N2"))) {
    sets <- vector("list", NB)
    for (b in seq_len(NB)) sets[[b]] <- if (nl=="N0") sample(pool, k)
      else if (nl=="N1") unlist(lapply(names(dn), function(dd)
             sample(pool_c[idx_by_dec[[dd]]], min(dn[[dd]], length(idx_by_dec[[dd]])))), use.names=FALSE)
      else pool_c[grow_set(dn, rho_t, C, idx_by_dec)]
    key_ <- paste(i, nl); names(sets) <- sprintf("%s|%s|s%04d", i, nl, seq_len(NB))
    banks[[key_]] <- sets; meta_b[[key_]] <- list(i=i, nl=nl, eng=eng_of[[d$target]], tgt=d$target)
  }
}
cat(sprintf("[%s] banks drawn: %d\n", COH, length(banks)))
bankauc <- list()
for (eng in c("wmean","gsva","ssgsea")) {
  ks <- names(banks)[vapply(meta_b, `[[`, "", "eng") == eng]; if (!length(ks)) next
  if (eng == "wmean") { for (kk in ks) bankauc[[kk]] <- score_sets(banks[[kk]], eng, weight_src=wts[[meta_b[[kk]]$tgt]]) }
  else { v <- score_sets(unlist(banks[ks], recursive=FALSE, use.names=FALSE) |> setNames(unlist(lapply(banks[ks], names))), eng)
         for (kk in ks) bankauc[[kk]] <- v[names(banks[[kk]])] }
  cat(sprintf("[%s] %s banks scored (%d)\n", COH, eng, length(ks)))
}
rows <- list()
for (kk in names(banks)) {
  i <- meta_b[[kk]]$i; nl <- meta_b[[kk]]$nl; d <- sub_[[i]]; bk <- bankauc[[kk]]
  a_obs <- with(V[i,], if (dirn == 1L) dir_obs else 1 - dir_obs)
  bd <- if (d$dirn == 1L) bk else 1 - bk; bf <- pmax(bk, 1 - bk)
  rows[[length(rows)+1]] <- data.frame(cohort=COH, unit=d$name, kind="decoy_own", target=d$target,
    null=if (PMODE == "natural") paste0(nl, "nat") else nl, dirn=d$dirn, fold_obs=max(a_obs,1-a_obs), fold_med=median(bf), fold_p95=unname(quantile(bf,.95)),
    dir_obs=if (d$dirn==1L) a_obs else 1-a_obs, dir_med=median(bd), dir_p95=unname(quantile(bd,.95)))
}
write.csv(do.call(rbind, rows), file.path(D, "calib_v3", paste0(COH, if (PMODE == "natural") "_ownnat.csv" else "_own.csv")), row.names=FALSE)
cat(sprintf("[%s] DONE\n", COH))
