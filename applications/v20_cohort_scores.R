# Round 15 parts B4-B6 (validation/PLAN_round15_literature_cohorts_controls.md, amendment 3): in every built eligible
# cohort of a type, the round-13 rule (L) and the directed C-index of
#   null models (B3: fitted in TCGA alone; higher score, shorter survival),
#   published models of the type (the frozen table; their published direction),
#   their co-expressed and random-gene decoys of A3 (their coefficients and directions as drawn there).
# A model or decoy is scorable in a cohort when at least 80% of its genes are measured there (round 13 A5), and is
# scored on those genes. -> icb/v20_round15/cohort_scores_<TYPE>.csv.gz
#   DFDR_LIB=<frozen install> Rscript v20_cohort_scores.R <TYPE> [cores]
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/v20_common.R")
args <- commandArgs(TRUE); tp <- args[1]
NC <- if (length(args) > 1) as.integer(args[2]) else 8L
CO <- read.csv(file.path(G, "signature-decoy-fdr", "validation", "round15_cohorts", "built_cohorts_r15.csv"),
               stringsAsFactors = FALSE)
CO <- CO[CO$type == tp & CO$built, ]
if (!nrow(CO)) { cat("no built cohort for", tp, "\n"); quit(save = "no") }
load_cohort <- function(id, dir) {
  e <- as.matrix(read.delim(gzfile(file.path(dir, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  cl <- read.delim(file.path(dir, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample), name = id)
}
num <- function(s) as.numeric(strsplit(s, ";")[[1]])
units <- list()
nm <- read.csv(gzfile(file.path(OUT, sprintf("null_models_%s.csv.gz", tp))), stringsAsFactors = FALSE)
nm <- nm[nm$n_model >= 2, ]
for (i in seq_len(nrow(nm)))
  units[[length(units) + 1]] <- list(kind = "null", id = sprintf("%s_null%03d", tp, nm$pool_index[i]), model = NA,
                                     class = nm$class[i], genes = strsplit(nm$genes[i], ";")[[1]], w = num(nm$coefs[i]),
                                     dirn = 1)
pm <- published_models(tp)
for (i in seq_len(nrow(pm)))
  units[[length(units) + 1]] <- list(kind = "target", id = pm$model_id[i], model = pm$model_id[i], class = NA,
                                     genes = strsplit(pm$genes_current[i], ";")[[1]], w = num(pm$coefs_current[i]),
                                     dirn = 1)
R <- read.csv(gzfile(file.path(OUT, sprintf("models_%s.csv.gz", tp))), stringsAsFactors = FALSE)
DI <- read.csv(gzfile(file.path(OUT, sprintf("models_decoy_info_%s.csv.gz", tp))), stringsAsFactors = FALSE)
R <- R[R$kind == "decoy", ]; R$genes <- DI$genes[match(R$decoy, DI$decoy)]; stopifnot(!anyNA(R$genes))
for (i in seq_len(nrow(R)))
  units[[length(units) + 1]] <- list(kind = "decoy", id = R$decoy[i], model = R$model[i], class = R$class[i],
                                     genes = strsplit(R$genes[i], ";")[[1]], w = num(R$weights[i]), dirn = R$dirn[i])
cat(tp, ":", length(units), "units over", nrow(CO), "cohorts\n")
out <- list()
for (j in seq_len(nrow(CO))) {
  co <- load_cohort(CO$cohort[j], CO$dir[j]); print(co)
  r <- parallel::mclapply(units, function(u) {
    k <- u$genes %in% rownames(co$Z)
    base <- data.frame(type = tp, cohort = CO$cohort[j], kind = u$kind, id = u$id, model = u$model, class = u$class,
                       k = length(u$genes), measured = sum(k))
    if (mean(k) < 0.8 || sum(k) < 2) return(cbind(base, scorable = FALSE, cdir = NA, p = NA, worse = NA, validated = NA))
    s <- u$dirn * lin_score(co, u$genes[k], u$w[k])
    ce <- cindex_rows(s, co$time, co$event); kv <- km_valid(s, co$time, co$event)
    cbind(base, scorable = TRUE, cdir = 1 - ce, p = kv[["p"]], worse = kv[["worse"]],
          validated = as.numeric(isTRUE(kv[["p"]] < 0.05) && isTRUE(kv[["worse"]] == 1)))
  }, mc.cores = NC, mc.preschedule = TRUE)
  out[[j]] <- do.call(rbind, r)
  cat(CO$cohort[j], "done\n")
}
S <- do.call(rbind, out)
write.csv(S, gzfile(file.path(OUT, sprintf("cohort_scores_%s.csv.gz", tp))), row.names = FALSE)
cat("done", tp, nrow(S), "rows\n")
