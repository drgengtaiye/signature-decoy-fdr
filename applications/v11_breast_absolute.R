#!/usr/bin/env Rscript
# Post hoc sensitivity analysis of Application B (not in the plan; prompted by v11_breast_pc1_strength.R): decoys
# matched on the mean absolute pairwise correlation of their genes, the coherence a first-principal-component score
# depends on, instead of the mean signed correlation. Null banks are drawn with the seed of the primary analysis and
# are therefore the same; only the decoys differ. Decoy seed 20261025.
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
detach("package:decoyfdr", unload = TRUE)          # this analysis needs grow_decoys(coherence =), installed in Rlib4
suppressMessages(library(decoyfdr, lib.loc = "/path/to/analysis_root/icb/v11_pkg/Rlib4"))
stopifnot("coherence" %in% names(formals(grow_decoys)))
OUT <- file.path(D, "v11_pkg"); NC <- 13L
co <- breast_cohorts(); sets <- breast_sets()
nmeas <- sapply(sets, function(g) min(vapply(co, function(c) sum(g %in% rownames(c$Z)), 1L)))
keep <- names(nmeas)[nmeas >= 3]
lib <- dfdr_library(sets[keep], engine = "pc1")
dec <- grow_decoys(co, lib, exclude = unique(unlist(sets)), n_coex = 32, n_rand = 8, seed = 20261025L, cores = NC,
                   coherence = "absolute")
print(dec); saveRDS(dec, file.path(OUT, "breast_abs_decoys.rds"))
write.csv(dec$info, file.path(OUT, "breast_abs_decoy_info.csv"), row.names = FALSE)
fit <- fit_units(co, lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
write.csv(fit$units, gzfile(file.path(OUT, "breast_abs_units.csv.gz")), row.names = FALSE)
calls <- call_units(fit$units, stats = "fold")
write.csv(calls, file.path(OUT, "breast_abs_calls.csv"), row.names = FALSE)
i <- dec$info; ok <- abs(i$r_achieved_med - i$r_target_med) <= 0.05
F <- rbind(decoy_fdr(calls, "fold", "coex"), decoy_fdr(calls, "fold", "rand"),
           cbind(decoy_fdr(calls, "fold", "coex", decoy_filter = calls$kind == "target" | calls$unit %in% i$decoy[ok]),
                 class = "coex, closest")[, -2])
print(F); write.csv(F, file.path(OUT, "breast_abs_fdr.csv"), row.names = FALSE)
