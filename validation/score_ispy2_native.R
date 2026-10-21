# Native scores of the 52 signatures for the I-SPY2 contrasts, with the development scorer (icb/native_score.R, score_all).
#!/usr/bin/env Rscript
# Re-score every signature with the METHOD ITS OWN PUBLICATION USED
#   GSVA (35), ssGSEA (6), weighted mean (12)   -> 53 of 59
#   6 signatures need bespoke algorithms from the original papers
#   (COX-IS, IPRES, IPS, IPSOV, PassON, PredictIO) and are excluded, noted in the report.
#
# Runs both sides: ICB response (13 cohorts) and TCGA prognosis (20 cancer types),
# with the same matched-random-gene-set null.
suppressMessages({library(GSVA); library(matrixStats); library(survival)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
P <- "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
NRAND <- 200

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

score_all <- function(M, key) {
  # M: genes x samples, log-scale. Returns samples x signatures.
  Z <- { s <- rowSds(M,na.rm=TRUE); s[s==0|is.na(s)] <- 1; (M-rowMeans(M,na.rm=TRUE))/s }
  out <- list(); gsva_sets <- list(); ssgsea_sets <- list()
  for (i in seq_len(nrow(meta))) {
    nm <- meta$signature[i]; mth <- meta$method[i]
    if (nm %in% SPEC) next
    S <- load_sig(nm, key); if(is.null(S)) next
    g <- intersect(S$genes, rownames(M)); if(length(g) < 3) next
    if (grepl("Weighted", mth)) {
      idx <- match(g, S$genes); w <- S$weights[idx]
      if (all(is.na(w))) w <- rep(1, length(g))
      w[is.na(w)] <- 0
      if (sum(abs(w)) == 0) next
      out[[nm]] <- as.numeric(crossprod(Z[g,,drop=FALSE], w)) / sum(abs(w))
    } else if (grepl("ssGSEA", mth)) ssgsea_sets[[nm]] <- g
    else gsva_sets[[nm]] <- g                                  # GSVA (default)
  }
  if (length(gsva_sets))  { r <- gsva(gsvaParam(M, gsva_sets, minSize=3), verbose=FALSE)
                            for(n in rownames(r)) out[[n]] <- as.numeric(r[n,]) }
  if (length(ssgsea_sets)){ r <- gsva(ssgseaParam(M, ssgsea_sets, minSize=3), verbose=FALSE)
                            for(n in rownames(r)) out[[n]] <- as.numeric(r[n,]) }
  as.data.frame(out, check.names=FALSE)
}

auc <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_)
  r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  a<-(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0); max(a,1-a) }
cidx <- function(v,y){ np<--as.numeric(v); c<-as.numeric(concordance(y~np)$concordance); max(c,1-c) }


V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
for (nm in c("ISPY2_pembro", "ISPY2_durva")) {
  M <- as.matrix(read.delim(gzfile(file.path(V, "cohorts", paste0(nm, ".expr.tsv.gz"))), row.names=1, check.names=FALSE))
  M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
  Sc <- score_all(M, "symbol"); rownames(Sc) <- colnames(M)
  write.csv(Sc, file.path(V, "cohorts", paste0(nm, ".native_scores.csv")))
  cat(nm, nrow(M), "genes,", ncol(Sc), "signatures scored\n")
}
