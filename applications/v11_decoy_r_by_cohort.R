#!/usr/bin/env Rscript
# Per signature and cohort: the signature's internal correlation and the median over its co-expressed decoys of
# theirs, from the decoys saved by grow_decoys().  Rscript v11_decoy_r_by_cohort.R tcga|breast|response
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg"); app <- commandArgs(trailingOnly = TRUE)[1]
dec <- readRDS(file.path(OUT, sprintf("%s_decoys.rds", app)))
co <- dec$cohorts
rows <- do.call(rbind, lapply(dec$decoys[vapply(dec$decoys, `[[`, "", "class") == "coex"], function(d)
  data.frame(target = d$target, cohort = co, r_target = d$r_target, r_decoy = d$r_achieved)))
R <- aggregate(cbind(r_target, r_decoy) ~ target + cohort, rows, median)
names(R)[names(R) == "r_decoy"] <- "r_decoy_med"
write.csv(R, file.path(OUT, sprintf("%s_decoy_r_by_type.csv", app)), row.names = FALSE)
cat(sprintf("%s: %d signature x cohort pairs; median target r %.3f, decoy r %.3f; decoys within 0.05: %.0f%%\n", app,
            nrow(R), median(R$r_target), median(R$r_decoy_med), 100 * mean(abs(rows$r_decoy - rows$r_target) <= 0.05, na.rm = TRUE)))
