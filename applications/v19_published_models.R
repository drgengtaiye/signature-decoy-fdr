# Round 13 part A (validation/PLAN_round13_published_risk_models.md, A6-A8): each frozen published model of a type and
# its 32 co-expressed and 8 random-gene decoys scored in TCGA and in the external cohort of the type.
#   DFDR_LIB=<frozen install> Rscript v19_published_models.R <TYPE> [cores]
# A decoy carries its model's coefficients in random order (Mersenne-Twister, seed 20261115 + type index, drawn decoy
# by decoy in the order of grow_decoys) and the direction of its own association with survival in TCGA (plan A7).
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/v19_common.R")
args <- commandArgs(TRUE); tp <- args[1]
NC <- if (length(args) > 1) as.integer(args[2]) else 8L
ti <- match(tp, TYPES) - 1L; stopifnot(!is.na(ti))
pm <- published_models(tp); stopifnot(!is.null(pm))
if (!nrow(pm)) { cat("no scored model for", tp, "\n"); quit(save = "no") }
P <- load_pair(tp); print(P$tcga); print(P$ext)
univ <- intersect(rownames(P$tcga$Z), rownames(P$ext$Z))
sets <- setNames(strsplit(pm$genes_current, ";"), pm$model_id)
W <- setNames(lapply(strsplit(pm$coefs_current, ";"), as.numeric), pm$model_id)
stopifnot(all(lengths(sets) == lengths(W)), all(vapply(sets, function(g) all(g %in% univ), TRUE)))

rows <- list()
for (m in names(sets)) {
  st <- lin_score(P$tcga, sets[[m]], W[[m]]); se <- lin_score(P$ext, sets[[m]], W[[m]])
  rows[[m]] <- data.frame(type = tp, model = m, kind = "target", class = NA, decoy = NA, k = length(sets[[m]]),
                          dirn = 1, t(score_stats(st, se, P, dirn = 1)))
}
big <- names(sets)[lengths(sets) >= 3]                   # grow_decoys needs three genes
lib <- dfdr_library(sets[big], engine = "wmean", weights = W[big])
dec <- grow_decoys(list(P$tcga, P$ext), lib, exclude = unique(unlist(theme_pools())), n_coex = 32L, n_rand = 8L,
                   n_sets = 1L, seed = 20261115L + ti, cores = NC, max_bytes = 100e9)
print(dec)
write.csv(dec$info, gzfile(file.path(OUT, sprintf("models_decoy_info_%s.csv.gz", tp))), row.names = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261115L + ti)
wd <- lapply(dec$decoys, function(d) { w <- W[[d$target]]; stopifnot(length(w) == d$k); w[sample.int(length(w))] })
dr <- parallel::mclapply(seq_along(dec$decoys), function(i) {
  d <- dec$decoys[[i]]
  st <- lin_score(P$tcga, d$genes, wd[[i]]); se <- lin_score(P$ext, d$genes, wd[[i]])
  dirn <- if (cindex_rows(st, P$tcga$time, P$tcga$event) <= 0.5) 1 else -1   # C <= 1/2: higher score, shorter survival
  data.frame(type = tp, model = d$target, kind = "decoy", class = d$class, decoy = d$name, k = d$k, dirn = dirn,
             t(score_stats(st, se, P, dirn = dirn)), weights = paste(signif(wd[[i]], 8), collapse = ";"))
}, mc.cores = NC, mc.preschedule = TRUE)
R <- rbind(do.call(rbind, rows) |> transform(weights = NA), do.call(rbind, dr))
write.csv(R, gzfile(file.path(OUT, sprintf("models_%s.csv.gz", tp))), row.names = FALSE)
cat("done", tp, nrow(R), "rows\n")
