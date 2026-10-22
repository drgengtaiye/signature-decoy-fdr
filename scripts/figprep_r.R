# Data for Figures 1d and 4b: internal mean pairwise correlation of every signature and of its set-1 decoys in each
# response and TCGA cohort, from the decoy objects of the package runs (applications/v11_*.R).
#   Rscript figprep_r.R   -> icb/v13_fig/r_by_cohort_<app>.csv
for (app in c("response", "tcga")) {
  d <- readRDS(sprintf("/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v11_pkg/%s_decoys.rds", app))
  rows <- lapply(d$decoys, function(x) data.frame(decoy = x$name, target = x$target, class = x$class,
                                                   cohort = d$cohorts, r_target = x$r_target, r_decoy = x$r_achieved))
  out <- do.call(rbind, rows)
  write.csv(out, sprintf("/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v13_fig/r_by_cohort_%s.csv", app), row.names = FALSE)
  cat(app, nrow(out), "\n")
}
