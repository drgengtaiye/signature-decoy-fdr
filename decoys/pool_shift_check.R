#!/usr/bin/env Rscript
# Is the N1/N2 null shift caused by the draw pool? null_matched_v3_coexpr.R draws its matched
# null sets from pool_c = every signature gene + a decile-stratified sample of the rest, capped
# at 8,000, so signature (largely immune) genes are over-represented in the very pool the null
# is drawn from. Here the N1 draw (decile-matched, the step N2 builds on) is repeated from three
# pools - (i) that pool, (ii) the same size drawn with signature genes at their natural share,
# (iii) the same size with no signature gene - and the signed AUC of each null set is read with
# an unweighted mean of z-scores (cheap; the shift is a property of the pool, not the engine).
#   Rscript pool_shift_check.R -> pool_shift_check.csv
suppressMessages({library(matrixStats); library(MultiAssayExperiment); library(SummarizedExperiment); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
COH <- c("ICB_Mariathasan","ICB_Braun","ICB_Gide","ICB_Liu","ICB_Kim","ICB_Jung",
         "ICB_Hugo","ICB_Riaz","ICB_Van_Allen","ICB_Miao1","ICB_Puch","ICB_Limagne1")
NB <- 300L; POOLCAP <- 8000L
meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
SPEC <- meta$signature[grepl("Specific", meta$method)]
rg <- function(nm, key) { f <- file.path(D,"SignatureSets-main/data-raw/IO",paste0(nm,".csv")); if(!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors=FALSE, fileEncoding="latin1")
  g <- if(key=="ensembl" && "gene_id" %in% names(d)) sub("\\..*$","",d$gene_id) else d$gene_name; unique(g[!is.na(g)&g!=""]) }
auc_s <- function(y,v){ r<-rank(v); n1<-sum(y==1); n0<-sum(y==0); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
one <- function(co) {
  set.seed(5L + sum(utf8ToInt(co)))
  x <- readRDS(file.path(D, paste0(co, ".rds"))); ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]; if (is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[,keep,drop=FALSE]; yv <- y[keep]
  rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
  M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
  key <- if (mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
  Z <- { s<-rowSds(M); s[s==0]<-1; (M-rowMeans(M))/s }; pool <- rownames(M)
  mu <- rowMeans(M); brk <- quantile(mu, probs=seq(0,1,by=.1)); brk[1] <- -Inf; brk[11] <- Inf
  dec <- cut(mu, breaks=brk, labels=FALSE, include.lowest=TRUE); names(dec) <- pool
  sigs <- setdiff(meta$signature, SPEC); G <- lapply(setNames(sigs,sigs), function(s) intersect(rg(s,key), pool))
  G <- G[vapply(G, length, 1L) >= 3]; sg <- unique(unlist(G))
  strat <- function(v, n) unlist(lapply(split(v, dec[v]), function(u) sample(u, min(length(u), ceiling(n/10)))), use.names=FALSE)
  P <- list(enriched = unique(c(sg, strat(setdiff(pool, sg), POOLCAP - length(sg)))),
            natural  = strat(pool, POOLCAP),
            sig_free = strat(setdiff(pool, sg), POOLCAP))
  g1 <- vapply(pool, function(g) auc_s(yv, Z[g,]), 1.0)
  out <- list()
  for (pn in names(P)) { pc <- P[[pn]]; ib <- split(pc, dec[pc])
    for (s in names(G)) { dn <- table(factor(dec[G[[s]]], levels=1:10)); dn <- dn[dn>0]
      a <- vapply(seq_len(NB), function(b) { gg <- unlist(lapply(names(dn), function(d) sample(ib[[d]], min(dn[[d]], length(ib[[d]])))))
             auc_s(yv, colMeans(Z[gg,,drop=FALSE])) }, 1.0)
      out[[length(out)+1]] <- data.frame(cohort=co, pool=pn, signature=s, sig_share=mean(pc %in% sg), shift=median(a)-0.5) } }
  out[[length(out)+1]] <- data.frame(cohort=co, pool="single_gene_sig", signature=NA, sig_share=NA, shift=mean(g1[sg])-0.5)
  out[[length(out)+1]] <- data.frame(cohort=co, pool="single_gene_rest", signature=NA, sig_share=NA, shift=mean(g1[setdiff(pool,sg)])-0.5)
  do.call(rbind, out)
}
R <- do.call(rbind, mclapply(COH, one, mc.cores=12))
write.csv(R, file.path(D, "pool_shift_check.csv"), row.names=FALSE)
cat("signature-gene share of the draw pool and median shift of the N1 bank (signed AUC - 0.5):\n")
for (pn in c("enriched","natural","sig_free")) { r <- R[R$pool==pn,]
  pc <- tapply(r$shift, r$cohort, mean)
  cat(sprintf("  %-9s share %.3f  mean shift %+.4f  (>0 in %d of 12 cohorts; tests >0 %.0f%%)\n", pn, mean(r$sig_share),
              mean(r$shift), sum(pc > 0), 100*mean(r$shift > 0))) }
r <- R[R$pool=="single_gene_sig",]; q <- R[R$pool=="single_gene_rest",]
cat(sprintf("single genes: signature genes mean AUC-0.5 %+.4f (>0 in %d/12 cohorts); all other genes %+.4f\n",
            mean(r$shift), sum(r$shift>0), mean(q$shift)))
