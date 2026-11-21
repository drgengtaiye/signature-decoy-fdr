#!/usr/bin/env Rscript
# Analysis 1 of validation/PLAN_calibration_round5.md: a signed null for the weighted-mean signatures whose
# curated weights are not all equal.
#
# null_matched_v2.R, null_matched_prog.R and null_matched_icbsurv.R score a weighted-mean signature as
# crossprod(Z, w)/sum(|w|) but its random sets as colMeans(Z), the same score only when every weight is 1.
# Four signatures carry weights of +1 and -1 (EMT_Thompson, IMPRES_Auslander, IRG_Yang, MPS_Guijarro). For
# each of them, in each cohort (and TCGA type), 200 sets are drawn at the signature's own rung K from the same
# pool; round(K m / k) genes chosen at random get weight -1 and the rest +1 (m of the k measured genes are
# negative), and the set is scored by the signature's own formula. Only the null columns of those rows change;
# the observed statistic is recomputed and must equal the stored one.
#
#   Rscript v10_signed_null.R response|prognosis|icbsurv   -> v10_signed_<setting>.csv
#   (response also writes scores/<cohort>_signed.csv.gz, the signed sets' per-patient scores)
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment); library(survival); library(parallel)})
MODE <- commandArgs(trailingOnly=TRUE)[1]
D <- "/path/to/analysis_root/icb"
P <- "/path/to/analysis_root/pancan"
NRAND <- 200
LADDER <- c(4, 7, 12, 20, 32, 50, 80, 130, 260, 440)
SIGNED <- c("EMT_Thompson", "IMPRES_Auslander", "IRG_Yang", "MPS_Guijarro")
meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature)
load_sig <- function(nm, key="symbol") {
  d <- read.csv(file.path(D,"SignatureSets-main/data-raw/IO", paste0(nm,".csv")), stringsAsFactors=FALSE,
                fileEncoding="latin1")
  g <- if(key=="ensembl" && "gene_id" %in% names(d)) sub("\\..*$","",d$gene_id) else d$gene_name
  w <- suppressWarnings(as.numeric(d$weight)); ok <- !is.na(g) & g!=""
  list(genes=g[ok], weights=w[ok])
}
auc <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_)
  r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  a<-(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0); max(a,1-a) }
cidx_prog <- function(v,y){ np <- -as.numeric(v); c <- as.numeric(concordance(y~np)$concordance); max(c,1-c) }
cidx_icb <- function(v, y) { ok <- is.finite(v); if (length(unique(y[ok, 1])) < 2) return(NA_real_)
  c_ <- as.numeric(concordance(y[ok] ~ I(-as.numeric(v[ok])))$concordance); max(c_, 1 - c_) }

# signature score and its signed null scores (sets x patients) in one cohort
signed <- function(Z, pool, s_, key, cohort) {
  S <- load_sig(s_, key); g <- intersect(S$genes, pool); if (length(g) < 3) return(NULL)
  w <- S$weights[match(g, S$genes)]; w[is.na(w)] <- 0
  obs <- as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))
  k <- length(g); m <- sum(w < 0); K <- LADDER[which.min(abs(LADDER - k))]
  set.seed(20261021L + sum(utf8ToInt(cohort)) + sum(utf8ToInt(s_)))
  nul <- t(vapply(seq_len(NRAND), function(i) {
    gg <- sample(pool, K); ww <- rep(1, K); ww[sample.int(K, round(K*m/k))] <- -1
    as.numeric(crossprod(Z[gg,,drop=FALSE], ww))/K }, numeric(ncol(Z))))
  list(obs=obs, null=nul, k=k, m=m, K=K)
}
row_of <- function(v) c(rand_med=median(v,na.rm=TRUE), rand_p95=unname(quantile(v,.95,na.rm=TRUE)),
                        rand_mean=mean(v,na.rm=TRUE), rand_sd=sd(v,na.rm=TRUE))

load_icb <- function(f) {
  x <- suppressWarnings(readRDS(f)); ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]
  if (is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
  list(M=as.matrix(assay(ex[[pick]])), cd=as.data.frame(colData(x)))
}
finish <- function(M) {
  rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
  M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
  key <- if(mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
  Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }
  list(M=M, Z=Z, key=key, pool=rownames(M))
}

rows <- list()
if (MODE == "response") {
  ref <- read.csv(file.path(D, "native_response_matched_v2.csv"))
  for (nm in unique(ref$cohort)) {
    L <- load_icb(file.path(D, paste0(nm, ".rds"))); cd <- L$cd
    r <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y) <- rownames(cd)
    keep <- intersect(colnames(L$M), names(y)[!is.na(y)]); yv <- y[keep]
    F <- finish(L$M[, keep, drop=FALSE]); sc <- list()
    for (s_ in SIGNED) {
      q <- signed(F$Z, F$pool, s_, F$key, nm); if (is.null(q)) next
      rr <- ref[ref$cohort==nm & ref$signature==s_, ]; stopifnot(nrow(rr)==1)
      a <- auc(yv, q$obs); stopifnot(abs(a - rr$AUC) < 1e-10, q$k == rr$n_genes)
      v <- apply(q$null, 1, function(z) auc(yv, z))
      sc[[s_]] <- setNames(as.data.frame(t(q$null)), sprintf("wmean_signed|%s|r%04d", s_, seq_len(NRAND)))
      rows[[length(rows)+1]] <- data.frame(cohort=nm, signature=s_, n_genes=q$k, n_negative=q$m,
        null_size=q$K, AUC=a, t(row_of(v)), old_med=rr$rand_med, old_p95=rr$rand_p95, old_sd=rr$rand_sd)
    }
    # per-patient scores of the signed sets, for the label-permutation scripts (v10_signed_ranks.py)
    if (length(sc)) write.csv(cbind(sample=colnames(F$Z), do.call(cbind, unname(sc))),
                              gzfile(file.path(D, "scores", paste0(nm, "_signed.csv.gz"))), row.names=FALSE)
    cat("response", nm, "\n")
  }
} else if (MODE == "icbsurv") {
  ref <- read.csv(file.path(D, "native_icbsurv_matched.csv"))
  for (nm in unique(ref$cohort)) {
    L <- load_icb(file.path(D, paste0(nm, ".rds"))); cd <- L$cd; M <- L$M
    eps <- list()
    for (ep in c("OS","PFS")) {
      tc <- paste0("survival_time_", tolower(ep)); ec <- paste0("event_occurred_", tolower(ep))
      if (!all(c(tc, ec) %in% colnames(cd))) next
      tm <- suppressWarnings(as.numeric(cd[[tc]])); ev <- suppressWarnings(as.numeric(cd[[ec]]))
      keep <- intersect(colnames(M), rownames(cd)[is.finite(tm) & tm > 0 & is.finite(ev)])
      if (length(keep) >= 25 && sum(ev[match(keep, rownames(cd))], na.rm=TRUE) >= 10)
        eps[[ep]] <- list(keep=keep, y=Surv(tm[match(keep, rownames(cd))], ev[match(keep, rownames(cd))]))
    }
    keep_all <- sort(unique(unlist(lapply(eps, `[[`, "keep"))))
    F <- finish(M[, keep_all, drop=FALSE])
    for (s_ in SIGNED) {
      q <- signed(F$Z, F$pool, s_, F$key, nm); if (is.null(q)) next
      names(q$obs) <- colnames(F$Z); colnames(q$null) <- colnames(F$Z)
      for (ep in names(eps)) {
        rr <- ref[ref$cohort==nm & ref$signature==s_ & ref$endpoint==ep, ]; if (nrow(rr) != 1) next
        kp <- eps[[ep]]$keep; yy <- eps[[ep]]$y
        cc <- cidx_icb(q$obs[kp], yy); stopifnot(abs(cc - rr$C) < 1e-10)
        v <- apply(q$null[, kp, drop=FALSE], 1, function(z) cidx_icb(z, yy))
        rows[[length(rows)+1]] <- data.frame(cohort=nm, endpoint=ep, signature=s_, n_genes=q$k,
          n_negative=q$m, null_size=q$K, C=cc, t(row_of(v)), old_med=rr$rand_med, old_p95=rr$rand_p95,
          old_sd=rr$rand_sd)
      }
    }
    cat("icbsurv", nm, "\n")
  }
} else if (MODE == "prognosis") {
  ref <- read.csv(file.path(D, "native_prognosis_matched.csv"))
  surv <- read.delim(file.path(P,"Survival_SupplementalTable_S1_20171025_xena_sp"), check.names=FALSE,
                     stringsAsFactors=FALSE)
  ENDPOINTS <- list(OS=c("OS","OS.time"), DSS=c("DSS","DSS.time"), PFI=c("PFI","PFI.time"))
  one <- function(ca) {
    e <- as.matrix(read.delim(gzfile(file.path(P,"expr",paste0(ca,".gz"))), row.names=1, check.names=FALSE))
    colnames(e) <- gsub("\\.","-",colnames(e)); e <- e[,substr(colnames(e),14,15)=="01",drop=FALSE]
    s <- surv[match(colnames(e),surv$sample),]
    ok <- !is.na(s$OS) & !is.na(s$OS.time) & s$OS.time>0
    e <- e[,ok,drop=FALSE]; s <- s[ok,]
    e <- e[rowSums(is.finite(e))==ncol(e) & rowSds(e)>0,,drop=FALSE]
    Z <- { sd <- rowSds(e); sd[sd==0] <- 1; (e-rowMeans(e))/sd }
    out <- list()
    for (s_ in SIGNED) {
      q <- signed(Z, rownames(e), s_, "symbol", ca); if (is.null(q)) next
      for (ep in names(ENDPOINTS)) {
        rr <- ref[ref$cancer==ca & ref$signature==s_ & ref$endpoint==ep, ]; if (nrow(rr) != 1) next
        ev <- suppressWarnings(as.numeric(s[[ENDPOINTS[[ep]][1]]])); tm <- suppressWarnings(as.numeric(s[[ENDPOINTS[[ep]][2]]]))
        keep <- which(!is.na(ev) & !is.na(tm) & tm > 0); yy <- Surv(tm[keep], ev[keep])
        cc <- cidx_prog(q$obs[keep], yy); stopifnot(abs(cc - rr$C) < 1e-10)
        v <- apply(q$null[, keep, drop=FALSE], 1, function(z) cidx_prog(z, yy))
        out[[length(out)+1]] <- data.frame(cancer=ca, endpoint=ep, signature=s_, n_genes=q$k,
          n_negative=q$m, null_size=q$K, C=cc, t(row_of(v)), old_med=rr$rand_med, old_p95=rr$rand_p95,
          old_sd=rr$rand_sd)
      }
    }
    cat("prognosis", ca, "\n"); do.call(rbind, out)
  }
  rows <- mclapply(unique(ref$cancer), one, mc.cores=20, mc.preschedule=FALSE)
  bad <- vapply(rows, function(r) inherits(r, "try-error") || is.null(r), TRUE)
  if (any(bad)) print(rows[bad])
}
R <- do.call(rbind, rows)
write.csv(R, file.path(D, paste0("v10_signed_", MODE, ".csv")), row.names=FALSE)
cat(sprintf("%s: %d rows written\n", MODE, nrow(R)))
