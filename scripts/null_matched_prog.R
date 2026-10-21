#!/usr/bin/env Rscript
# Prognosis benchmark with a null matched on BOTH gene-set size and scoring engine.
#
# native_score.R compared every signature in a cancer type with one shared null of
# 48-gene sets scored as mean z. That is the construction the response side abandoned
# (null_matched.R): signatures have a median of 17 genes and are scored by GSVA, ssGSEA
# or a weighted mean, and the spread of a null depends on set size. This script applies
# the identical correction to survival in the TCGA types, so the two sides of the
# dissociation are measured against the same kind of null.
#
# Signatures and random sets are scored once per cancer type, on the same primary-tumour
# samples used by native_score.R; overall survival, disease-specific survival and the
# progression-free interval are then each evaluated on the samples that have that endpoint.
suppressMessages({library(GSVA); library(matrixStats); library(survival); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
P <- "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
NRAND  <- 200
LADDER <- c(4, 7, 12, 20, 32, 50, 80, 130, 260, 440)
NCORE  <- as.integer(Sys.getenv("NCORE", "12"))
ENDPOINTS <- list(OS=c("OS","OS.time"), DSS=c("DSS","DSS.time"), PFI=c("PFI","PFI.time"))
OUTDIR <- file.path(D, "null_matched_prog_parts"); dir.create(OUTDIR, showWarnings=FALSE)

meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
SPEC <- meta$signature[grepl("Specific", meta$method)]
load_sig <- function(nm) {
  f <- file.path(D,"SignatureSets-main/data-raw/IO", paste0(nm,".csv"))
  if(!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors=FALSE, fileEncoding="latin1")
  g <- d$gene_name
  w <- if("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_,nrow(d))
  ok <- !is.na(g) & g!=""
  list(genes=g[ok], weights=w[ok])
}
engine <- function(m) if (grepl("Weighted",m)) "wmean" else if (grepl("ssGSEA",m)) "ssgsea" else "gsva"
cidx <- function(v,y){ np <- -as.numeric(v); c <- as.numeric(concordance(y~np)$concordance); max(c,1-c) }

surv <- read.delim(file.path(P,"Survival_SupplementalTable_S1_20171025_xena_sp"),
                   check.names=FALSE, stringsAsFactors=FALSE)
cancers <- sub("\\.gz$","",list.files(file.path(P,"expr"),pattern="\\.gz$"))

run_one <- function(ca) {
  set.seed(42)
  part <- file.path(OUTDIR, paste0(ca, ".csv"))
  if (file.exists(part)) return(ca)
  e <- try(as.matrix(read.delim(gzfile(file.path(P,"expr",paste0(ca,".gz"))),
                                row.names=1, check.names=FALSE)), silent=TRUE)
  if (inherits(e,"try-error")) return(NULL)
  colnames(e) <- gsub("\\.","-",colnames(e)); e <- e[,substr(colnames(e),14,15)=="01",drop=FALSE]
  s <- surv[match(colnames(e),surv$sample),]
  ok <- !is.na(s$OS) & !is.na(s$OS.time) & s$OS.time>0          # as native_score.R
  e <- e[,ok,drop=FALSE]; s <- s[ok,]
  e <- e[rowSums(is.finite(e))==ncol(e) & rowSds(e)>0,,drop=FALSE]
  if (ncol(e)<80 || sum(s$OS)<25) return(NULL)
  Z <- { sd <- rowSds(e); sd[sd==0] <- 1; (e-rowMeans(e))/sd }
  pool <- rownames(e)

  # observed scores, each signature by its own publication's method
  obs <- list(); use <- list(); gsva_sets <- list(); ssg_sets <- list()
  for (i in seq_len(nrow(meta))) {
    s_ <- meta$signature[i]; mth <- meta$method[i]
    if (s_ %in% SPEC) next
    S <- load_sig(s_); if (is.null(S)) next
    g <- intersect(S$genes, pool); if (length(g) < 3) next
    if (engine(mth)=="wmean") {
      w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1,length(g))
      w[is.na(w)] <- 0; if (sum(abs(w))==0) next
      obs[[s_]] <- as.numeric(crossprod(Z[g,,drop=FALSE], w))/sum(abs(w))
    } else if (engine(mth)=="ssgsea") ssg_sets[[s_]] <- g else gsva_sets[[s_]] <- g
    use[[s_]] <- list(k=length(g), eng=engine(mth))
  }
  if (length(gsva_sets)) { rr <- gsva(gsvaParam(e, gsva_sets, minSize=3), verbose=FALSE)
                           for (n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_,]) }
  if (length(ssg_sets))  { rr <- gsva(ssgseaParam(e, ssg_sets, minSize=3), verbose=FALSE)
                           for (n_ in rownames(rr)) obs[[n_]] <- as.numeric(rr[n_,]) }

  # null scores: a ladder of sizes, per engine, scored exactly as the signatures are
  nullscore <- list()                       # engine -> sets x samples
  for (eng in unique(sapply(use[names(obs)], `[[`, "eng"))) {
    sets <- list()
    for (k in LADDER) for (i in seq_len(NRAND)) sets[[sprintf("k%03d_r%04d", k, i)]] <- sample(pool, k)
    if (eng == "wmean") {
      nullscore[[eng]] <- t(sapply(sets, function(g) colMeans(Z[g,,drop=FALSE])))
    } else {
      prm <- if (eng=="ssgsea") ssgseaParam(e, sets, minSize=3) else gsvaParam(e, sets, minSize=3)
      nullscore[[eng]] <- gsva(prm, verbose=FALSE)
    }
  }

  rows <- list()
  for (ep in names(ENDPOINTS)) {
    ev <- suppressWarnings(as.numeric(s[[ENDPOINTS[[ep]][1]]]))
    tm <- suppressWarnings(as.numeric(s[[ENDPOINTS[[ep]][2]]]))
    keep <- which(!is.na(ev) & !is.na(tm) & tm > 0)
    if (length(keep) < 80 || sum(ev[keep]) < 25) next
    y <- Surv(tm[keep], ev[keep])
    nulls <- list()
    for (eng in names(nullscore)) {
      M <- nullscore[[eng]]
      v <- apply(M[, keep, drop=FALSE], 1, function(z) cidx(z, y))
      for (k in LADDER) nulls[[paste(eng,k)]] <- as.numeric(v[grepl(sprintf("^k%03d_", k), names(v))])
    }
    for (s_ in names(obs)) {
      u <- use[[s_]]; k <- LADDER[which.min(abs(LADDER - u$k))]
      v <- nulls[[paste(u$eng, k)]]; if (is.null(v) || !length(v)) next
      rows[[paste(ep, s_)]] <- data.frame(cancer=ca, endpoint=ep, n=length(keep),
          events=sum(ev[keep]), signature=s_, n_genes=u$k, engine=u$eng, null_size=k,
          C=cidx(obs[[s_]][keep], y), rand_med=median(v, na.rm=TRUE),
          rand_p95=quantile(v, .95, na.rm=TRUE, names=FALSE),
          # mean and SD of the same null draw, so the prognostic and response sides can
          # be compared in units of their own null spread
          rand_mean=mean(v, na.rm=TRUE), rand_sd=sd(v, na.rm=TRUE))
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, part, row.names=FALSE)
  cat(sprintf("PROG %-6s n=%4d rows=%3d\n", ca, ncol(e), nrow(out)))
  ca
}

done <- mclapply(cancers, function(ca) tryCatch(run_one(ca), error=function(err) {
  cat("FAILED", ca, conditionMessage(err), "\n"); NULL }), mc.cores=NCORE, mc.preschedule=FALSE)
parts <- list.files(OUTDIR, pattern="\\.csv$", full.names=TRUE)
R <- do.call(rbind, lapply(parts, read.csv))
write.csv(R, file.path(D,"native_prognosis_matched.csv"), row.names=FALSE)
cat(sprintf("\ndone: %d cancer types, %d rows\n", length(parts), nrow(R)))
