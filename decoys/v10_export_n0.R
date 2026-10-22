#!/usr/bin/env Rscript
# Write the N0 banks of null_matched_v3_coexpr.R (signed AUC, 1,000 draws per test) as one table for Python.
#   Rscript v10_export_n0.R  -> icb/v10_n0_bank.csv.gz  (cohort, signature, a0001..a1000)
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
out <- list()
for (f in list.files(file.path(D, "null_v2"), pattern="_v2\\.rds$", full.names=TRUE)) {
  r <- readRDS(f)
  for (nm in grep(" N0$", names(r$null_auc), value=TRUE)) {
    v <- r$null_auc[[nm]]; stopifnot(length(v) == 1000)
    out[[length(out)+1]] <- c(list(cohort=r$cohort, signature=sub(" N0$", "", nm)),
                              setNames(as.list(v), sprintf("a%04d", 1:1000)))
  }
}
R <- do.call(rbind.data.frame, out)
write.csv(R, gzfile(file.path(D, "v10_n0_bank.csv.gz")), row.names=FALSE)
cat(nrow(R), "tests\n")
