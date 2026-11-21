#!/usr/bin/env Rscript
# V3 of validation/PLAN_round6_applications.md: the response benchmark rerun end to end with decoyfdr,
# decoys and null banks drawn afresh (seed 20261022), beside the published calls and FDRs.
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg")
co <- response_cohorts()
for (c in co) print(c)
L <- immune_library(key = "ensembl")
print(L$lib)
t0 <- Sys.time()
res <- run_decoy_fdr(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, n_draw = 1000, keep_draws = 0,
                     seed = 20261022L, cores = 12L, B = 6000L)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(res)
saveRDS(res$decoys, file.path(OUT, "response_decoys.rds"))
write.csv(res$fit$units, gzfile(file.path(OUT, "response_units.csv.gz")), row.names = FALSE)
write.csv(res$calls, file.path(OUT, "response_calls.csv"), row.names = FALSE)
write.csv(res$fdr, file.path(OUT, "response_fdr.csv"), row.names = FALSE)
write.csv(res$decoys$info, file.path(OUT, "response_decoy_info.csv"), row.names = FALSE)
