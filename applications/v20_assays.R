# Round 15 part D (validation/PLAN_round15_literature_cohorts_controls.md D1-D2, amendment 4): each clinical assay and its 32 co-expressed and 8 random-gene
# decoys, scored and judged as in round 13 A6-A9 in TCGA and every built eligible cohort of its type. Decoys are grown over TCGA and the round-11 cohort of the
# type from genes in no assay, no eligible model of the type and no theme pool (seed 20261209 + assay index). A unit is scorable in a cohort when at least 80%
# of its genes are measured there, and is scored on those genes.
# -> icb/v20_round15/assay_scores.csv.gz, assay_decoy_info.csv.gz
#   DFDR_LIB=<frozen install> Rscript v20_assays.R [cores]
source("/path/to/analysis_root/signature-decoy-fdr/applications/v20_common.R")
# (the TCGA pair is held in PR: common.R keeps the TCGA data path in P, which a second load_pair() call reads)
args <- commandArgs(TRUE); NC <- if (length(args)) as.integer(args[1]) else 8L
AD <- file.path(G, "signature-decoy-fdr", "validation", "round15_assays")
A <- read.csv(file.path(AD, "assays_r15_frozen.csv"), stringsAsFactors = FALSE)
ALIAS <- read.csv(file.path(AD, "assays_alias_r15.csv"), stringsAsFactors = FALSE)
ER <- read.csv(file.path(AD, "er_status_r15.csv"), stringsAsFactors = FALSE)
stopifnot(all(is.finite(A$weight)), !anyDuplicated(paste(A$assay, A$gene_symbol)))
CO <- read.csv(file.path(G, "signature-decoy-fdr", "validation", "round15_cohorts", "built_cohorts_r15.csv"), stringsAsFactors = FALSE)
load_cohort <- function(id, dir, keep = NULL) {
  e <- as.matrix(read.delim(gzfile(file.path(dir, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  cl <- read.delim(file.path(dir, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  if (!is.null(keep)) { cl <- cl[cl$sample %in% keep, ]; e <- e[, cl$sample, drop = FALSE] }
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample), name = id)
}
# a gene that a cohort file carries under an older symbol counts as measured (amendment 4): its row takes the current symbol
apply_alias <- function(co, genes) {
  old <- rownames(co$Z); rn <- old
  for (g in genes) if (!g %in% rn) { al <- ALIAS$ALIAS[ALIAS$gene_symbol == g]; hit <- al[al %in% rn]; if (length(hit)) rn[match(hit[1], rn)] <- g }
  if (identical(rn, old)) return(co)
  for (nm in names(co)) { x <- co[[nm]]
    if (is.matrix(x) && identical(rownames(x), old)) rownames(co[[nm]]) <- rn else if (is.atomic(x) && identical(names(x), old)) names(co[[nm]]) <- rn }
  co
}
assays <- unique(A$assay); all_assay_genes <- unique(A$gene_symbol); out <- list(); info <- list()
for (ai in seq_along(assays)) {
  a <- assays[ai]; x <- A[A$assay == a, ]; tp <- x$cancer[1]; stopifnot(tp %in% c("BRCA", "COAD"))
  PR <- load_pair(tp); PR$tcga <- apply_alias(PR$tcga, x$gene_symbol); PR$ext <- apply_alias(PR$ext, x$gene_symbol); print(PR$tcga); print(PR$ext)
  pm <- published_models(tp)
  excl <- unique(c(unlist(theme_pools()), unlist(strsplit(pm$genes_current, ";")), all_assay_genes))
  univ <- intersect(rownames(PR$tcga$Z), rownames(PR$ext$Z))
  k <- x$gene_symbol %in% univ
  cat(a, ":", nrow(x), "genes,", sum(k), "measured in both cohorts used to grow decoys; missing:", x$gene_symbol[!k], "\n")
  g <- x$gene_symbol[k]; w <- x$weight[k]
  lib <- dfdr_library(setNames(list(g), a), engine = "wmean", weights = setNames(list(w), a))
  dec <- grow_decoys(list(PR$tcga, PR$ext), lib, exclude = excl, n_coex = 32L, n_rand = 8L, n_sets = 1L, seed = 20261209L + ai - 1L, cores = NC, max_bytes = 100e9)
  print(dec); info[[a]] <- transform(dec$info, assay = a)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261209L + ai - 1L)
  wd <- lapply(dec$decoys, function(d) { stopifnot(length(w) == d$k); w[sample.int(length(w))] })
  dn <- vapply(seq_along(dec$decoys), function(i) {           # the direction of its own association with survival in TCGA (plan A7)
    if (cindex_rows(lin_score(PR$tcga, dec$decoys[[i]]$genes, wd[[i]]), PR$tcga$time, PR$tcga$event) <= 0.5) 1 else -1 }, 1)
  units <- c(list(list(kind = "target", name = a, class = NA, genes = x$gene_symbol, w = x$weight, dirn = 1)),
             lapply(seq_along(dec$decoys), function(i) list(kind = "decoy", name = dec$decoys[[i]]$name, class = dec$decoys[[i]]$class, genes = dec$decoys[[i]]$genes, w = wd[[i]], dirn = dn[i])))
  cc <- CO[CO$type == tp & CO$built, ]
  cohorts <- c(list(TCGA = PR$tcga), setNames(lapply(seq_len(nrow(cc)), function(j) load_cohort(cc$cohort[j], cc$dir[j])), cc$cohort))
  # ER-positive patients, secondary (amendment 4): annotated cohorts with at least 25 events among them
  for (j in seq_len(nrow(cc))) {
    ep <- ER$sample[ER$cohort == cc$cohort[j] & ER$er == 1]
    if (length(ep)) { cx <- load_cohort(cc$cohort[j], cc$dir[j], keep = ep); if (sum(cx$event) >= 25) cohorts[[paste0(cc$cohort[j], "_ERpos")]] <- cx }
  }
  for (nm in names(cohorts)) {
    co <- apply_alias(cohorts[[nm]], x$gene_symbol)
    r <- parallel::mclapply(units, function(u) {
      m <- u$genes %in% rownames(co$Z)
      base <- data.frame(assay = a, type = tp, cohort = nm, kind = u$kind, unit = u$name, class = u$class, k = length(u$genes), measured = sum(m))
      if (mean(m) < 0.8 || sum(m) < 2) return(cbind(base, scorable = FALSE, cdir = NA, p = NA, worse = NA, validated = NA))
      s <- u$dirn * lin_score(co, u$genes[m], u$w[m])
      ce <- cindex_rows(s, co$time, co$event); kv <- km_valid(s, co$time, co$event)
      cbind(base, scorable = TRUE, cdir = 1 - ce, p = kv[["p"]], worse = kv[["worse"]], validated = as.numeric(isTRUE(kv[["p"]] < 0.05) && isTRUE(kv[["worse"]] == 1)))
    }, mc.cores = NC, mc.preschedule = TRUE)
    out[[length(out) + 1]] <- do.call(rbind, r)
    cat(a, nm, "done\n")
  }
}
write.csv(do.call(rbind, out), gzfile(file.path(OUT, "assay_scores.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, info), gzfile(file.path(OUT, "assay_decoy_info.csv.gz")), row.names = FALSE)
cat("done\n")
