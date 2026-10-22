#!/usr/bin/env Rscript
# Application B of validation/PLAN_round6_applications.md: the breast-cancer prognostic signatures of Venet et al.
# (2011) in the breast cohorts, scored by their first principal component (primary) and by the mean z (sensitivity);
# folded C-index. Seeds derive from 20261022: decoys 20261022, null banks 20261023, global null 20261024.
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg"); NC <- 13L
co <- breast_cohorts(); for (c in co) print(c)
sets <- breast_sets()
nmeas <- sapply(sets, function(g) min(vapply(co, function(c) sum(g %in% rownames(c$Z)), 1L)))
cat("sets with fewer than 3 genes in some cohort (excluded):", paste(names(nmeas)[nmeas < 3], collapse = ", "), "\n")
keep <- names(nmeas)[nmeas >= 3]
lib <- dfdr_library(sets[keep], engine = "pc1"); print(lib)
t0 <- Sys.time()
dec <- grow_decoys(co, lib, exclude = unique(unlist(sets)), n_coex = 32, n_rand = 8, seed = 20261022L, cores = NC)
print(dec); saveRDS(dec, file.path(OUT, "breast_decoys.rds"))
write.csv(dec$info, file.path(OUT, "breast_decoy_info.csv"), row.names = FALSE)
cat(sprintf("decoys in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
close_ <- function(calls) {
  i <- dec$info; ok <- abs(i$r_achieved_med - i$r_target_med) <= 0.05
  calls$kind == "target" | calls$unit %in% i$decoy[ok]
}
fdr_table <- function(calls) rbind(decoy_fdr(calls, "fold", "coex"), decoy_fdr(calls, "fold", "rand"),
  cbind(decoy_fdr(calls, "fold", "coex", decoy_filter = close_(calls))[, ], class = "coex, closest")[, -2])
fit <- fit_units(co, lib, dec, n_draw = 1000L, keep_draws = 200L, seed = 20261023L, cores = NC)
saveRDS(fit, file.path(OUT, "breast_pc1_fit.rds"))
write.csv(fit$units, gzfile(file.path(OUT, "breast_pc1_units.csv.gz")), row.names = FALSE)
calls <- call_units(fit$units, stats = "fold")
write.csv(calls, file.path(OUT, "breast_pc1_calls.csv"), row.names = FALSE)
F <- fdr_table(calls); print(F); write.csv(F, file.path(OUT, "breast_pc1_fdr.csv"), row.names = FALSE)
cat(sprintf("fit in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cal <- calibrate_null(fit, B = 200L, stats = "fold", seed = 20261024L, cores = NC)
print(cal$summary)
for (nm in names(cal)) write.csv(cal[[nm]], file.path(OUT, sprintf("breast_pc1_globalnull_%s.csv", nm)), row.names = FALSE)
rm(fit); gc()
libm <- dfdr_library(sets[keep], engine = "mean")
fm <- fit_units(co, libm, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
write.csv(fm$units, gzfile(file.path(OUT, "breast_mean_units.csv.gz")), row.names = FALSE)
cm <- call_units(fm$units, stats = "fold")
write.csv(cm, file.path(OUT, "breast_mean_calls.csv"), row.names = FALSE)
Fm <- fdr_table(cm); print(Fm); write.csv(Fm, file.path(OUT, "breast_mean_fdr.csv"), row.names = FALSE)
cat(sprintf("all done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
