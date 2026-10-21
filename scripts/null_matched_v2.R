#!/usr/bin/env Rscript
# Response benchmark against a null matched on BOTH gene-set size and scoring engine.
#
# Identical in method to null_matched.R, with two changes:
#  1. signature metadata come from sig_methods_clean.csv. read.csv() on sig_methods.csv
#     silently loses three rows (WNT-beta-catenin_Takeuchi, CTL_Jiang and
#     Tcell-inflamed_Ayers) because the free-text description column contains quoted
#     line breaks; the clean file holds only the signature and method columns.
#  2. cohorts run in parallel, each with its own fixed seed, so a result does not depend
#     on the order in which cohorts or signatures are processed.
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
NRAND  <- 200
LADDER <- c(4, 7, 12, 20, 32, 50, 80, 130, 260, 440)
NCORE  <- as.integer(Sys.getenv("NCORE", "7"))
OUTDIR <- file.path(D, "null_matched_v2_parts"); dir.create(OUTDIR, showWarnings=FALSE)

meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(meta) == 59)
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
auc <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_)
  r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  a<-(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0); max(a,1-a) }
engine <- function(m) if (grepl("Weighted",m)) "wmean" else if (grepl("ssGSEA",m)) "ssgsea" else "gsva"

run_cohort <- function(f) {
  nm <- sub("\\.rds$","",basename(f))
  part <- file.path(OUTDIR, paste0(nm, ".csv")); if (file.exists(part)) return(nm)
  set.seed(sum(utf8ToInt(nm)))
  x <- try(suppressWarnings(readRDS(f)), silent=TRUE); if(inherits(x,"try-error")) return(NULL)
  ex <- try(experiments(x),silent=TRUE); if(inherits(ex,"try-error")) return(NULL)
  an <- names(ex); pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]
  if(is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
  if(is.na(pick)) return(NULL)
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  if(!("response" %in% colnames(cd))) return(NULL)
  r <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y)<-rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); if(length(keep)<25) return(NULL)
  M <- M[,keep,drop=FALSE]; yv <- y[keep]
  if(sum(yv)<5 || sum(1-yv)<5) return(NULL)
  rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
  M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
  if(nrow(M)<3000) return(NULL)
  key <- if(mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
  Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }
  pool <- rownames(M)

  obs <- list(); use <- list(); gsva_sets <- list(); ssg_sets <- list()
  for (i in seq_len(nrow(meta))) {
    s_ <- meta$signature[i]; mth <- meta$method[i]
    if (s_ %in% SPEC) next
    S <- load_sig(s_, key); if(is.null(S)) next
    g <- intersect(S$genes, pool); if(length(g) < 3) next
    if (engine(mth)=="wmean") {
      w <- S$weights[match(g, S$genes)]; if(all(is.na(w))) w <- rep(1,length(g))
      w[is.na(w)] <- 0; if(sum(abs(w))==0) next
      obs[[s_]] <- as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))
    } else if (engine(mth)=="ssgsea") ssg_sets[[s_]] <- g else gsva_sets[[s_]] <- g
    use[[s_]] <- list(k=length(g), eng=engine(mth))
  }
  if (length(gsva_sets)) { rr <- gsva(gsvaParam(M, gsva_sets, minSize=3), verbose=FALSE)
                           for(n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_,]) }
  if (length(ssg_sets))  { rr <- gsva(ssgseaParam(M, ssg_sets, minSize=3), verbose=FALSE)
                           for(n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_,]) }

  nulls <- list()
  for (eng in sort(unique(sapply(use[names(obs)], `[[`, "eng")))) {
    for (k in LADDER) {
      if (k > length(pool)) next
      sets <- lapply(seq_len(NRAND), function(i) sample(pool, k))
      names(sets) <- sprintf("r%04d", seq_along(sets))
      if (eng == "wmean") {
        v <- sapply(sets, function(g) auc(yv, colMeans(Z[g,,drop=FALSE])))
      } else {
        prm <- if (eng=="ssgsea") ssgseaParam(M, sets, minSize=3) else gsvaParam(M, sets, minSize=3)
        rr <- gsva(prm, verbose=FALSE); v <- apply(rr, 1, function(z) auc(yv, z))
      }
      nulls[[paste(eng,k)]] <- as.numeric(v)
    }
  }
  rows <- list()
  for (s_ in names(obs)) {
    u <- use[[s_]]; k <- LADDER[which.min(abs(LADDER - u$k))]
    v <- nulls[[paste(u$eng, k)]]; if (is.null(v)) next
    rows[[s_]] <- data.frame(cohort=nm, n=length(yv), signature=s_, n_genes=u$k,
        engine=u$eng, null_size=k, AUC=auc(yv, obs[[s_]]),
        rand_med=median(v,na.rm=TRUE), rand_p95=quantile(v,.95,na.rm=TRUE,names=FALSE),
        # mean and SD of the same null draw, so an effect can be expressed in units of
        # its own null spread rather than as a raw difference in concordance
        rand_mean=mean(v,na.rm=TRUE), rand_sd=sd(v,na.rm=TRUE))
  }
  out <- do.call(rbind, rows); write.csv(out, part, row.names=FALSE)
  cat(sprintf("RESP %-18s n=%3d sigs=%2d\n", nm, length(yv), nrow(out))); nm
}

files <- list.files(D, pattern="^ICB_.*\\.rds$", full.names=TRUE)
invisible(mclapply(files, function(f) tryCatch(run_cohort(f), error=function(e) {
  cat("FAILED", f, conditionMessage(e), "\n"); NULL }), mc.cores=NCORE, mc.preschedule=FALSE))
R <- do.call(rbind, lapply(list.files(OUTDIR, pattern="\\.csv$", full.names=TRUE), read.csv))
write.csv(R, file.path(D,"native_response_matched_v2.csv"), row.names=FALSE)
cat(sprintf("\ndone: %d cohorts, %d signatures, %d rows\n", length(unique(R$cohort)),
            length(unique(R$signature)), nrow(R)))
