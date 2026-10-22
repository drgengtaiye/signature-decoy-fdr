#!/usr/bin/env Rscript
# Application A of validation/PLAN_round6_applications.md: the 52 immune signatures and their decoys for
# overall survival in 20 TCGA cancer types (primary), then disease-specific survival and the
# progression-free interval with the same decoys (secondary). Seeds derive from 20261022: decoys 20261022,
# null banks 20261023, global null 20261024.
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg"); NC <- 20L
L <- immune_library(key = "symbol")
CL <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
clusters <- setNames(CL$cluster_fixed, CL$signature)
fdr_table <- function(calls) do.call(rbind, c(
  lapply(c("fold", "dir"), function(st) do.call(rbind, lapply(c("coex", "rand"), function(cl) decoy_fdr(calls, st, cl)))),
  lapply(c("fold", "dir"), function(st) cbind(decoy_fdr(calls, st, "coex", clusters = clusters), boot = "clusters"))))

co <- tcga_cohorts("OS")
cat(length(co), "cancer types\n"); for (c in co) print(c)
t0 <- Sys.time()
dec <- grow_decoys(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, seed = 20261022L, cores = NC,
                   max_bytes = 200e9)
print(dec); saveRDS(dec, file.path(OUT, "tcga_decoys.rds"))
write.csv(dec$info, file.path(OUT, "tcga_decoy_info.csv"), row.names = FALSE)
cat(sprintf("decoys in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
fit <- fit_units(co, L$lib, dec, n_draw = 1000L, keep_draws = 200L, seed = 20261023L, cores = NC)
saveRDS(fit, file.path(OUT, "tcga_OS_fit.rds"))
write.csv(fit$units, gzfile(file.path(OUT, "tcga_OS_units.csv.gz")), row.names = FALSE)
calls <- call_units(fit$units)
write.csv(calls, file.path(OUT, "tcga_OS_calls.csv"), row.names = FALSE)
F <- fdr_table(calls); print(F)
write.csv(F, file.path(OUT, "tcga_OS_fdr.csv"), row.names = FALSE)
cat(sprintf("fit in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cal <- calibrate_null(fit, B = 200L, seed = 20261024L, cores = NC)
print(cal$summary)
for (nm in names(cal)) write.csv(cal[[nm]], file.path(OUT, sprintf("tcga_OS_globalnull_%s.csv", nm)), row.names = FALSE)
cat(sprintf("calibration in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rm(fit); gc()
for (ep in c("DSS", "PFI")) {
  ce <- tcga_cohorts(ep)
  fe <- fit_units(ce, L$lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
  write.csv(fe$units, gzfile(file.path(OUT, sprintf("tcga_%s_units.csv.gz", ep))), row.names = FALSE)
  ca <- call_units(fe$units)
  write.csv(ca, file.path(OUT, sprintf("tcga_%s_calls.csv", ep)), row.names = FALSE)
  Fe <- fdr_table(ca); cat(ep, "\n"); print(Fe)
  write.csv(Fe, file.path(OUT, sprintf("tcga_%s_fdr.csv", ep)), row.names = FALSE)
}
cat(sprintf("all done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
