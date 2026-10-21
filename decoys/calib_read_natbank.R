#!/usr/bin/env Rscript
# Read the v3 decoys (observed AUC from calib_v3/<cohort>_v3.csv, unchanged) and the 52 targets
# against the natural-pool N1/N2 banks of null_natpool_v4.R, in the format of calib_score_v3.R.
#   Rscript calib_read_natbank.R <ICB_cohort> -> calib_v3/<cohort>_nat.csv
COH <- commandArgs(trailingOnly=TRUE)[1]
D   <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
BK  <- readRDS(file.path(D, "null_v4nat", paste0(COH, "_nat.rds")))
V   <- read.csv(file.path(D, "calib_v3", paste0(COH, "_v3.csv")), stringsAsFactors=FALSE)
V   <- V[V$null == "N0",]
rd <- function(a, b, dirn) { bd <- if (dirn == 1L) b else 1 - b; bf <- pmax(b, 1 - b)
  c(fold_obs=max(a,1-a), fold_med=median(bf), fold_p95=unname(quantile(bf,.95)),
    dir_obs=if (dirn==1L) a else 1-a, dir_med=median(bd), dir_p95=unname(quantile(bd,.95))) }
rows <- list()
for (i in seq_len(nrow(V))) {
  v <- V[i,]; a <- if (v$dirn == 1L) v$dir_obs else 1 - v$dir_obs
  if (v$kind == "target") a <- BK$obs[[v$unit]]
  for (nl in c("N1","N2")) { b <- BK$null_auc[[paste(v$target, nl)]]; if (is.null(b)) next
    dn <- if (v$kind == "target") BK$dirn[[v$unit]] else v$dirn
    rows[[length(rows)+1]] <- data.frame(v[c("cohort","n","unit","kind","target","pool","structure")],
      null=paste0(nl, "nat"), k=v$k, dirn=dn, r_int=v$r_int, t(rd(a, b, dn))) }
}
write.csv(do.call(rbind, rows), file.path(D, "calib_v3", paste0(COH, "_nat.csv")), row.names=FALSE)
cat(sprintf("[%s] DONE natural-pool reading\n", COH))
