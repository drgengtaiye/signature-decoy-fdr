#!/usr/bin/env Rscript
# Target-decoy calibration, step 1: decoys matched to each published signature.
#
# calib_decoys_v2.R grew 120 decoys to the median internal correlation of the library, by size
# bin, in the mean of the 12 per-cohort correlation matrices. They reached r 0.20-0.48 per
# cohort against the signatures' 0.31-0.60, so any false-positive rate read from them was a
# lower bound. Here every one of the 52 signatures (the targets) receives its own decoys,
# matched to it on size, engine, expression-decile composition and its internal correlation
# in EACH of the 12 cohorts separately. A decoy is then judged against its target's own null
# bank (null_v2/*_v2.rds), so target and decoy differ only in which genes carry the structure.
#
# Pools   A  genes in none of the 52 signatures
#         B  A minus every gene of the 90-signature TME library (no known immune programme)
# Per target: 32 coex decoys from A, 8 uniform random-gene decoys from A, 16 coex decoys from B.
#
#   Rscript calib_decoys_v3.R   -> calib_decoys_v3.rds, calib_decoys_v3_info.csv
suppressMessages({library(matrixStats); library(MultiAssayExperiment); library(SummarizedExperiment)
                  library(parallel)})
D   <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
COH <- c("ICB_Mariathasan","ICB_Braun","ICB_Gide","ICB_Liu","ICB_Kim","ICB_Jung",
         "ICB_Hugo","ICB_Riaz","ICB_Van_Allen","ICB_Miao1","ICB_Puch","ICB_Limagne1")
NA_COEX <- 32L; NA_RAND <- 8L; NB_COEX <- 16L
TOL     <- 0.02          # near-best band: any candidate within TOL of the best fit may be taken
NCORE   <- 24L
RNGkind("L'Ecuyer-CMRG"); set.seed(20261019L)

meta <- read.csv(file.path(D,"sig_methods_clean.csv"), stringsAsFactors=FALSE)
meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
SPEC <- meta$signature[grepl("Specific", meta$method)]
read_genes <- function(f, key) {
  d <- read.csv(f, stringsAsFactors=FALSE, fileEncoding="latin1")
  g <- if (key=="ensembl" && "gene_id" %in% names(d)) sub("\\..*$","",d$gene_id) else d$gene_name
  unique(g[!is.na(g) & g!=""])
}
prep <- function(co) {                                   # the null_matched_v3_coexpr.R filters
  x  <- readRDS(file.path(D, paste0(co, ".rds")))
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm","expr")[c("expr_gene_tpm","expr") %in% an][1]
  if (is.na(pick)) pick <- grep("^expr", an, value=TRUE)[1]
  M  <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r  <- as.character(cd$response); y <- ifelse(r=="R",1,ifelse(r=="NR",0,NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)])
  M <- M[,keep,drop=FALSE]
  rownames(M) <- sub("\\..*$","",rownames(M)); M <- M[!duplicated(rownames(M)),,drop=FALSE]
  M <- M[rowSums(is.finite(M))==ncol(M) & rowSds(M)>0,,drop=FALSE]
  mu  <- rowMeans(M)
  brk <- quantile(mu, probs=seq(0,1,by=.1), na.rm=TRUE); brk[1] <- -Inf; brk[11] <- Inf
  dec <- cut(mu, breaks=brk, labels=FALSE, include.lowest=TRUE); names(dec) <- rownames(M)
  s <- rowSds(M); s[s==0] <- 1
  list(Z=(M - rowMeans(M))/s, dec=dec, n=ncol(M))
}
mean_pair_r <- function(C, idx) { if (length(idx)<2) return(NA_real_)
  sub <- C[idx,idx,drop=FALSE]; (sum(sub)-sum(diag(sub)))/(length(idx)*(length(idx)-1)) }

# ---------------------------------------------------------------- cohorts, universe, pools
P <- list(); univ <- NULL
for (co in COH) { P[[co]] <- prep(co); g <- rownames(P[[co]]$Z)
                  univ <- if (is.null(univ)) g else intersect(univ, g) }
key <- if (mean(grepl("^ENSG", univ)) > 0.5) "ensembl" else "symbol"
T0  <- read.csv(file.path(D,"null_v3_tidy.csv"), stringsAsFactors=FALSE)
T0  <- T0[T0$null=="N0" & T0$cohort %in% COH,]
TGT <- sort(unique(T0$signature))
sgenes <- lapply(setNames(TGT,TGT), function(s) read_genes(file.path(D,"SignatureSets-main/data-raw/IO",paste0(s,".csv")), key))
allsig <- unique(unlist(lapply(setdiff(meta$signature, SPEC), function(s) {
  f <- file.path(D,"SignatureSets-main/data-raw/IO",paste0(s,".csv")); if (file.exists(f)) read_genes(f,key) })))
tme <- unique(unlist(lapply(list.files(file.path(D,"SignatureSets-main/data-raw/TME"), "\\.csv$", full.names=TRUE),
                            read_genes, key=key)))
poolA <- setdiff(univ, allsig); poolB <- setdiff(poolA, tme)
cat(sprintf("universe %d | signature genes %d | TME genes %d | pool A %d | pool B %d | targets %d\n",
            length(univ), length(intersect(allsig,univ)), length(intersect(tme,univ)),
            length(poolA), length(poolB), length(TGT)))

# reference decile: the median of each gene's decile across the 12 cohorts
DM <- sapply(COH, function(co) P[[co]]$dec[univ]); rownames(DM) <- univ
refdec <- setNames(as.integer(round(rowMedians(DM))), univ)

# per-cohort correlation among pool-A genes
Cl <- lapply(COH, function(co) { Zc <- P[[co]]$Z[poolA,,drop=FALSE]; tcrossprod(Zc)/(ncol(Zc)-1) })
names(Cl) <- COH
cat(sprintf("correlation matrices: 12 x %d genes (%.1f GB)\n", length(poolA),
            12*as.numeric(object.size(Cl[[1]]))/1e9))

# ---------------------------------------------------------------- the targets
tinfo <- lapply(TGT, function(s) {
  rows <- T0[T0$signature==s,]
  rt <- sapply(COH, function(co) { g <- intersect(sgenes[[s]], rownames(P[[co]]$Z))
          if (length(g)<2) return(NA_real_); Zg <- P[[co]]$Z[g,,drop=FALSE]
          mean_pair_r(tcrossprod(Zg)/(ncol(Zg)-1), seq_along(g)) })
  gu <- intersect(sgenes[[s]], univ)
  k  <- as.integer(round(median(rows$k)))
  dprof <- refdec[gu]
  if (length(dprof) < k) dprof <- c(dprof, sample(dprof, k-length(dprof), replace=TRUE))
  if (length(dprof) > k) dprof <- sample(dprof, k)
  list(sig=s, k=k, engine=rows$engine[1], rt=rt, dn=table(factor(dprof, levels=1:10)))
}); names(tinfo) <- TGT

# ---------------------------------------------------------------- multi-cohort seed-and-grow
grow_multi <- function(dn, rt, idx_by_dec) {
  use_c <- which(is.finite(rt)); rt <- rt[use_c]
  need <- dn[dn > 0]; S <- integer(0); k <- sum(need)
  cs <- matrix(0, length(use_c), length(poolA))
  while (length(S) < k) {
    av <- names(need)[need > 0]
    if (!length(av)) {                                    # a decile ran dry: borrow a neighbour
      av <- as.character(1:10); need[av[!av %in% names(need)]] <- 0L }
    d <- if (length(av)==1L) av else sample(av, 1L, prob=pmax(need[av],1e-9))
    cand <- setdiff(idx_by_dec[[d]], S)
    if (!length(cand)) { need[d] <- 0L; next }
    pick <- if (!length(S)) cand[sample.int(length(cand),1L)] else {
      L <- sqrt(colMeans((cs[, cand, drop=FALSE]/length(S) - rt)^2))
      ok <- which(L <= min(L) + TOL); cand[ok[sample.int(length(ok),1L)]] }
    S <- c(S, pick)
    for (j in seq_along(use_c)) cs[j,] <- cs[j,] + Cl[[use_c[j]]][, pick]
    need[d] <- max(0L, need[d] - 1L)
  }
  S
}
ibA <- split(seq_along(poolA), refdec[poolA])
inB <- poolA %in% poolB
ibB <- lapply(ibA, function(v) v[inB[v]])

one_target <- function(s) {
  t <- tinfo[[s]]; out <- list()
  mk <- function(idx, pool, struct, b) {
    g <- poolA[idx]
    r_ach <- sapply(COH, function(co) mean_pair_r(Cl[[co]], idx))
    list(name=sprintf("dcy_%s_%s_%s_%02d", pool, struct, s, b), target=s, pool=pool, structure=struct,
         genes=g, k=length(g), engine=t$engine, dirn=sample(c(1L,-1L),1L),
         r_target=t$rt, r_achieved=r_ach)
  }
  for (b in seq_len(NA_COEX)) out[[length(out)+1]] <- mk(grow_multi(t$dn, t$rt, ibA), "A", "coex", b)
  for (b in seq_len(NA_RAND)) out[[length(out)+1]] <- mk(sample(seq_along(poolA), t$k), "A", "rand", b)
  for (b in seq_len(NB_COEX)) out[[length(out)+1]] <- mk(grow_multi(t$dn, t$rt, ibB), "B", "coex", b)
  out
}
t0 <- Sys.time()
DEC <- unlist(mclapply(TGT, one_target, mc.cores=NCORE, mc.set.seed=TRUE), recursive=FALSE)
names(DEC) <- vapply(DEC, `[[`, character(1), "name")
cat(sprintf("%d decoys grown in %.1f min\n", length(DEC), as.numeric(difftime(Sys.time(),t0,units="mins"))))

# ---------------------------------------------------------------- matching report
info <- do.call(rbind, lapply(DEC, function(d) data.frame(decoy=d$name, target=d$target, pool=d$pool,
  structure=d$structure, k=d$k, engine=d$engine, dirn=d$dirn,
  r_target_med=median(d$r_target, na.rm=TRUE), r_achieved_med=median(d$r_achieved, na.rm=TRUE),
  r_profile_cor=suppressWarnings(cor(d$r_target, d$r_achieved, use="complete.obs")),
  r_absdev_med=median(abs(d$r_achieved - d$r_target), na.rm=TRUE))))
jac <- do.call(rbind, lapply(split(DEC, paste(info$target, info$pool, info$structure)), function(L) {
  if (length(L) < 2) return(NULL)
  jj <- combn(length(L), 2, function(ij) { a <- L[[ij[1]]]$genes; b <- L[[ij[2]]]$genes
          length(intersect(a,b))/length(union(a,b)) })
  data.frame(target=L[[1]]$target, pool=L[[1]]$pool, structure=L[[1]]$structure,
             jac_med=median(jj), jac_max=max(jj)) }))
cat("\nper-cohort internal r, median over decoys (target = the matched signature):\n")
for (ps in list(c("A","coex"), c("B","coex"), c("A","rand"))) {
  ii <- info$pool==ps[1] & info$structure==ps[2]
  cat(sprintf("  %s %-4s  target %.3f  achieved %.3f  |dev| %.3f  profile cor %.2f  within-target Jaccard %.3f (max %.2f)\n",
      ps[1], ps[2], median(info$r_target_med[ii]), median(info$r_achieved_med[ii]),
      median(info$r_absdev_med[ii]), median(info$r_profile_cor[ii], na.rm=TRUE),
      median(jac$jac_med[jac$pool==ps[1] & jac$structure==ps[2]]),
      max(jac$jac_max[jac$pool==ps[1] & jac$structure==ps[2]])))
}
saveRDS(list(decoys=DEC, info=info, jaccard=jac, cohorts=COH, key=key, poolA=poolA, poolB=poolB,
             targets=tinfo), file.path(D, "calib_decoys_v3.rds"))
write.csv(info, file.path(D, "calib_decoys_v3_info.csv"), row.names=FALSE)
cat("saved -> calib_decoys_v3.rds\n")
