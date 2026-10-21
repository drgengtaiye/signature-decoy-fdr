#!/usr/bin/env Rscript
# Flatten the three nulls of null_matched_v3_coexpr.R into one tidy table.
# For every cohort x signature x null, both statistics are formed:
#   folded    max(a, 1-a)        -- what the published analysis uses
#   directed  a or 1-a           -- by the curation's 'association' column
# The null is transformed the same way, so each statistic has its own null.
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
fs <- list.files(file.path(D,"null_v2"), pattern="_v2\\.rds$", full.names=TRUE)
rows <- list()
for (f in fs) {
  r <- readRDS(f)
  for (s in names(r$obs)) {
    a <- r$obs[[s]]; if (!is.finite(a)) next
    ds <- if (is.null(r$dirn[[s]]) || is.na(r$dirn[[s]])) 1L else r$dirn[[s]]
    for (nl in c("N0","N1","N2")) {
      v <- r$null_auc[[paste(s, nl)]]; if (is.null(v) || all(!is.finite(v))) next
      fo <- max(a, 1-a); fn <- pmax(v, 1-v)
      do_ <- if (ds > 0) a else 1-a
      dn  <- if (ds > 0) v else 1-v
      rows[[length(rows)+1]] <- data.frame(
        cohort=r$cohort, n=r$n, signature=s, null=nl,
        k=r$use[[s]]$k, engine=r$use[[s]]$eng, dirn=ds, ndraw=length(v),
        fold_obs=fo, fold_med=median(fn,na.rm=TRUE),
        fold_p95=as.numeric(quantile(fn,.95,na.rm=TRUE)), fold_sd=sd(fn,na.rm=TRUE),
        dir_obs=do_, dir_med=median(dn,na.rm=TRUE),
        dir_p95=as.numeric(quantile(dn,.95,na.rm=TRUE)), dir_sd=sd(dn,na.rm=TRUE))
    }
  }
  cat(sprintf("  %s\n", basename(f)))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(D,"null_v3_tidy.csv"), row.names=FALSE)
cat(sprintf("tidy: %d rows, %d cohorts, %d signatures\n", nrow(out),
            length(unique(out$cohort)), length(unique(out$signature))))
