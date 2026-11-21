#!/usr/bin/env Rscript
# Application A, resumed from the saved fit: v11_tcga_prognosis.R stopped after saving the fit, the per-unit table and
# the calls, at the FDR table (its two bootstrap rows did not bind); this script builds that table, runs the global
# null with decoyfdr as installed in Rlib2, and fits disease-specific survival and the progression-free interval with
# the same decoys. Seeds as in v11_tcga_prognosis.R.
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg"); NC <- 20L
L <- immune_library(key = "symbol")
CL <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
clusters <- setNames(CL$cluster_fixed, CL$signature)
fdr_table <- function(calls) do.call(rbind, c(
  lapply(c("fold", "dir"), function(st) do.call(rbind, lapply(c("coex", "rand"), function(cl)
    cbind(decoy_fdr(calls, st, cl), boot = "targets")))),
  lapply(c("fold", "dir"), function(st) cbind(decoy_fdr(calls, st, "coex", clusters = clusters), boot = "clusters"))))
dec <- readRDS(file.path(OUT, "tcga_decoys.rds"))
calls <- read.csv(file.path(OUT, "tcga_OS_calls.csv"), stringsAsFactors = FALSE)
F <- fdr_table(calls); print(F)
write.csv(F, file.path(OUT, "tcga_OS_fdr.csv"), row.names = FALSE)
t0 <- Sys.time()
fit <- readRDS(file.path(OUT, "tcga_OS_fit.rds"))
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
