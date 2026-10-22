#!/usr/bin/env Rscript
# Analysis 8 of validation/PLAN_round7_estimand_and_bias.md: the decoy pool at library scale, on the 20 TCGA types.
#   (a) pool size for the immune library and for MSigDB Hallmark, Hallmark + C2 canonical pathways, and Hallmark + C2 CP
#       + C5 GO biological process; (b) 4 co-expressed decoys per Hallmark set from each MSigDB pool (seed 740000) and
#       how close they come to the set's internal correlation; (c) the 50 Hallmark sets, mean z, on overall survival with
#       ten decoy sets (seeds 741000 decoys, 741001 null banks, 741002 global null). Seeds are arbitrary labels.
#   DFDR_LIB=<frozen 0.2.0 install> Rscript v12_hallmark.R [cores]
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
suppressMessages(library(msigdbr))
stopifnot(packageVersion("decoyfdr", lib.loc = RUNLIB) >= "0.2.0")
NC <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 20L }
OUT <- file.path(D, "v12_pkg"); t0 <- Sys.time()
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
co <- tcga_cohorts("OS")
univ <- Reduce(intersect, lapply(co, function(c) rownames(c$Z)))
ms <- function(...) { x <- msigdbr(species = "Homo sapiens", ...); split(x$gene_symbol, x$gs_name) }
H <- lapply(ms(collection = "H"), unique)
CP <- unlist(lapply(c("CP", "CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                    function(s) ms(collection = "C2", subcollection = s)), recursive = FALSE)
BP <- ms(collection = "C5", subcollection = "GO:BP")
L <- immune_library(key = "symbol")
libs <- list(immune = unique(c(unlist(L$lib$sets), L$exclude)), hallmark = unique(unlist(H)),
             hallmark_cp = unique(c(unlist(H), unlist(CP))), hallmark_cp_gobp = unique(c(unlist(H), unlist(CP), unlist(BP))))
nsets <- c(immune = length(L$lib$sets), hallmark = length(H), hallmark_cp = length(H) + length(CP),
           hallmark_cp_gobp = length(H) + length(CP) + length(BP))
A <- data.frame(library = names(libs), sets = nsets[names(libs)], library_genes_measured = vapply(libs, function(g) sum(univ %in% g), 1L),
                universe = length(univ), pool = vapply(libs, function(g) sum(!univ %in% g), 1L))
cat(sprintf("msigdbr %s\n", packageVersion("msigdbr"))); print(A)
write.csv(cbind(A, msigdbr = as.character(packageVersion("msigdbr"))), file.path(OUT, "hallmark_pool_sizes.csv"), row.names = FALSE)

# (b) decoys for each Hallmark set from the three MSigDB pools
Hs <- lapply(H, function(g) intersect(g, univ)); Hs <- Hs[vapply(Hs, length, 1L) >= 3]
libH <- dfdr_library(Hs, engine = "mean")
B <- list()
for (nm in c("hallmark", "hallmark_cp", "hallmark_cp_gobp")) {
  pool_n <- sum(!univ %in% libs[[nm]])
  if (pool_n < 3 * max(lengths(Hs))) { cat(nm, ": pool of", pool_n, "genes, too small to grow decoys of the largest sets\n") }
  ok <- if (pool_n >= 3 * max(lengths(Hs))) names(Hs) else names(Hs)[lengths(Hs) * 3 <= pool_n]
  if (!length(ok)) next
  d <- grow_decoys(co, dfdr_library(Hs[ok], engine = "mean"), exclude = setdiff(libs[[nm]], unlist(Hs)), n_coex = 4, n_rand = 0,
                   seed = 740000L, cores = NC, max_bytes = 300e9)
  i <- d$info; i$pool <- nm; i$pool_size <- length(d$pool); B[[nm]] <- i
  cat(nm, ": pool", length(d$pool), "genes;", mins(), "\n")
}
Bi <- do.call(rbind, B)
write.csv(Bi, gzfile(file.path(OUT, "hallmark_pool_decoys.csv.gz")), row.names = FALSE)
Bs <- do.call(rbind, lapply(split(Bi, Bi$pool), function(x) {
  per <- tapply(abs(x$r_achieved_med - x$r_target_med) <= 0.05, x$target, all)
  data.frame(pool = x$pool[1], pool_size = x$pool_size[1], targets = length(per),
             median_ratio = median(x$r_achieved_med / x$r_target_med), share_all_within_0.05 = mean(per))
}))
print(Bs); write.csv(Bs, file.path(OUT, "hallmark_pool_summary.csv"), row.names = FALSE)

# (c) the Hallmark library on overall survival
dec <- grow_decoys(co, libH, n_coex = 32, n_rand = 8, n_sets = 10, seed = 741000L, cores = NC, max_bytes = 300e9)
print(dec); write.csv(dec$info, gzfile(file.path(OUT, "hallmark_decoy_info.csv.gz")), row.names = FALSE)
cat("decoys", mins(), "\n")
fit <- fit_units(co, libH, dec, n_draw = 1000L, keep_draws = 200L, seed = 741001L, cores = NC)
write.csv(fit$units, gzfile(file.path(OUT, "hallmark_OS_units.csv.gz")), row.names = FALSE)
R <- simplify2array(lapply(fit$scores, function(s) cor(t(s$targets[names(Hs), ]), method = "spearman")))
hc <- stats::hclust(stats::as.dist(1 - abs(apply(R, c(1, 2), median))), method = "average")
cl <- stats::cutree(hc, h = 0.5); write.csv(data.frame(signature = names(cl), cluster = cl), file.path(OUT, "hallmark_clusters.csv"), row.names = FALSE)
cat(length(unique(cl)), "Hallmark clusters\n")
calls <- call_units(fit$units, stats = "fold")
write.csv(calls, gzfile(file.path(OUT, "hallmark_OS_calls.csv.gz")), row.names = FALSE)
F <- rbind(cbind(decoy_fdr(calls, "fold", "coex"), boot = "targets"), cbind(decoy_fdr(calls, "fold", "coex", clusters = cl), boot = "clusters"),
           cbind(decoy_fdr(calls, "fold", "rand", clusters = cl), boot = "clusters"))
print(F); write.csv(F, file.path(OUT, "hallmark_OS_fdr.csv"), row.names = FALSE)
cal <- calibrate_null(fit, B = 200L, stats = "fold", seed = 741002L, cores = NC)
print(cal$summary)
for (nm in names(cal)) write.csv(cal[[nm]], file.path(OUT, sprintf("hallmark_OS_globalnull_%s.csv", nm)), row.names = FALSE)
cat("all done in", mins(), "\n")
