# Acceptance test of check_models() (validation/TOOL_check_models_acceptance.md): round 15 part A and the published
# models and decoys of part B5, rerun through the packaged function from the same inputs and seeds.
#   DFDR_LIB=<install of decoyfdr 0.3.0> Rscript v24_check_models_repro.R <TYPE> [cores]
# -> icb/v24_tool/repro_{units,scores,decoys}_<TYPE>.csv.gz (compared with round 15 by decoys/v24_check_models_repro.py)
source("/path/to/analysis_root/signature-decoy-fdr/applications/v20_common.R")
stopifnot(packageVersion("decoyfdr", lib.loc = RUNLIB) >= "0.3.0")
args <- commandArgs(TRUE); tp <- args[1]
NC <- if (length(args) > 1) as.integer(args[2]) else 8L
ti <- match(tp, TYPES) - 1L; stopifnot(!is.na(ti))
OUT24 <- file.path(G, "icb", "v24_tool")
pm <- published_models(tp); stopifnot(nrow(pm) > 0)
num <- function(s) as.numeric(strsplit(s, ";")[[1]])
models <- setNames(lapply(seq_len(nrow(pm)), function(i) setNames(num(pm$coefs_current[i]), strsplit(pm$genes_current[i], ";")[[1]])),
                   pm$model_id)
PR <- load_pair(tp); print(PR$tcga); print(PR$ext)
CO <- read.csv(file.path(G, "signature-decoy-fdr", "validation", "round15_cohorts", "built_cohorts_r15.csv"), stringsAsFactors = FALSE)
CO <- CO[CO$type == tp & CO$built, ]
load_cohort <- function(id, dir) {                                   # as v20_cohort_scores.R
  e <- as.matrix(read.delim(gzfile(file.path(dir, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  cl <- read.delim(file.path(dir, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample), name = id)
}
stopifnot(PR$ext$name %in% CO$cohort)
val <- lapply(seq_len(nrow(CO)), function(j) if (CO$cohort[j] == PR$ext$name) PR$ext else load_cohort(CO$cohort[j], CO$dir[j]))
t0 <- Sys.time()
chk <- check_models(models, train = PR$tcga, validate = val, match = list(PR$tcga, PR$ext),
                    n_printed = setNames(pm$n_genes, pm$model_id), exclude = unique(unlist(theme_pools())),
                    n_coex = 32L, n_rand = 8L, seed = 20261202L + ti, cores = NC, max_bytes = 100e9)
print(chk$models[chk$models$status != "kept" | chk$models$genes_used != chk$models$genes_given, ])
cat(tp, ": check_models took", format(Sys.time() - t0), "\n")
for (x in c("units", "scores", "decoys"))
  write.csv(chk[[x]], gzfile(file.path(OUT24, sprintf("repro_%s_%s.csv.gz", x, tp))), row.names = FALSE)
cat("done", tp, "\n")
