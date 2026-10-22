#!/usr/bin/env Rscript
# Post hoc diagnostic for Application B (not in the plan): the decoys are matched on the mean pairwise correlation of
# their genes, the coherence that matters for a mean score; a first-principal-component score depends instead on the
# share of variance its first component carries, which a set mixing up- and down-regulated genes can have high with
# a mean correlation near zero. Per signature and cohort: that share for the signature and the median over its
# co-expressed and random-gene decoys.
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg")
co <- breast_cohorts(); sets <- breast_sets()
dec <- readRDS(file.path(OUT, "breast_decoys.rds"))
targets <- unique(dec$info$target)
share <- function(Z, g) { g <- g[g %in% rownames(Z)]; if (length(g) < 3) return(NA_real_)
  X <- Z[g, , drop = FALSE]; G <- if (nrow(X) <= ncol(X)) tcrossprod(X) else crossprod(X)
  ev <- if (nrow(G) > 200) RSpectra::eigs_sym(G, 1L, which = "LA")$values[1] else eigen(G, symmetric = TRUE, only.values = TRUE)$values[1]
  ev / sum(diag(G)) }
rows <- parallel::mclapply(co, function(c) do.call(rbind, lapply(targets, function(s) {
  ds <- dec$decoys[dec$info$target == s]
  cx <- vapply(ds[vapply(ds, `[[`, "", "class") == "coex"], function(d) share(c$Z, d$genes), 1.0)
  rx <- vapply(ds[vapply(ds, `[[`, "", "class") == "rand"], function(d) share(c$Z, d$genes), 1.0)
  data.frame(cohort = c$name, target = s, k = sum(sets[[s]] %in% rownames(c$Z)), pc1_target = share(c$Z, sets[[s]]),
             pc1_coex_med = median(cx, na.rm = TRUE), pc1_rand_med = median(rx, na.rm = TRUE))
})), mc.cores = 13L)
R <- do.call(rbind, rows)
write.csv(R, file.path(OUT, "breast_pc1_strength.csv"), row.names = FALSE)
cat(sprintf("share of variance on the first component, median over signature x cohort: signatures %.3f, co-expressed decoys %.3f, random-gene decoys %.3f; signature above its decoys' median in %.0f%%\n",
            median(R$pc1_target, na.rm = TRUE), median(R$pc1_coex_med, na.rm = TRUE), median(R$pc1_rand_med, na.rm = TRUE),
            100 * mean(R$pc1_target > R$pc1_coex_med, na.rm = TRUE)))
