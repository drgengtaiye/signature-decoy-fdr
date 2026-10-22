#!/usr/bin/env Rscript
# Analysis 4 of validation/PLAN_round7_estimand_and_bias.md: every application rerun with decoyfdr 0.2.0 and ten
# independent decoy sets, with the seeds of round 6 (decoys 20261022, null banks 20261023, global null 20261024;
# absolute-coherence decoys 20261025). Set 1 must reproduce round 6 (icb/v11_pkg) exactly; the check is printed.
# Seeds are arbitrary labels.
#   DFDR_LIB=<frozen 0.2.0 install> Rscript v12_rerun.R response|tcga|breast|breast_abs [cores]
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
stopifnot(packageVersion("decoyfdr", lib.loc = RUNLIB) >= "0.2.0")
args <- commandArgs(trailingOnly = TRUE); app <- args[1]; NC <- if (length(args) > 1) as.integer(args[2]) else 20L
OUT <- file.path(D, "v12_pkg"); OLD <- file.path(D, "v11_pkg")
K <- 10L
CL <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
imm_clusters <- setNames(CL$cluster_fixed, CL$signature)
t0 <- Sys.time()
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# FDR table: per statistic, co-expressed decoys with the target and the cluster bootstrap, random-gene decoys, and the
# closest-matched co-expressed decoys (internal r within 0.05 of the target's), each pooled over the ten sets
fdr_table <- function(calls, info, clusters, stats) {
  close_ids <- info$decoy[abs(info$r_achieved_med - info$r_target_med) <= 0.05]
  filt <- calls$kind == "target" | calls$unit %in% close_ids
  do.call(rbind, lapply(stats, function(st) rbind(
    cbind(decoy_fdr(calls, st, "coex"), decoys = "coex", boot = "targets"),
    cbind(decoy_fdr(calls, st, "coex", clusters = clusters), decoys = "coex", boot = "clusters"),
    cbind(decoy_fdr(calls, st, "rand"), decoys = "rand", boot = "targets"),
    cbind(decoy_fdr(calls, st, "rand", clusters = clusters), decoys = "rand", boot = "clusters"),
    cbind(decoy_fdr(calls, st, "coex", decoy_filter = filt), decoys = "coex closest", boot = "targets"),
    cbind(decoy_fdr(calls, st, "coex", decoy_filter = filt, clusters = clusters), decoys = "coex closest", boot = "clusters"))))
}
# the estimate of each set alone (point estimates), and pi0 at the sensitivity thresholds
by_set <- function(calls, stats) do.call(rbind, lapply(stats, function(st) do.call(rbind, lapply(seq_len(K), function(s)
  cbind(set = s, decoy_fdr(calls, st, "coex", B = 2L, decoy_filter = calls$kind == "target" | calls$set %in% s)[,
        c("stat", "called", "decoy_rate", "exp_false", "fdr", "pi0", "fdr_pi0")])))))
pi0_sens <- function(calls, stats, clusters) do.call(rbind, lapply(stats, function(st) do.call(rbind, lapply(c(0.25, 0.5, 0.75),
  function(l) cbind(lambda = l, decoy_fdr(calls, st, "coex", clusters = clusters, lambda = l)[,
        c("stat", "called", "fdr", "pi0", "fdr_pi0", "fdr_pi0_lo", "fdr_pi0_hi")])))))

# set 1 against round 6: identical decoys, identical calls, identical estimate
check_set1 <- function(calls, old_calls_file, old_info_file, info, stats) {
  oc <- read.csv(file.path(OLD, old_calls_file), stringsAsFactors = FALSE)
  oi <- read.csv(file.path(OLD, old_info_file), stringsAsFactors = FALSE)
  i1 <- info[info$set == 1, ]
  g_ok <- identical(i1$decoy, oi$decoy) && identical(i1$genes, oi$genes) && identical(as.integer(i1$dirn), as.integer(oi$dirn))
  n1 <- calls[calls$kind == "target" | calls$set %in% 1L, ]
  key <- function(x) paste(x$unit, x$stat)
  m <- match(key(oc), key(n1))
  c_ok <- !anyNA(m) && identical(as.logical(oc$call), as.logical(n1$call[m])) &&
    isTRUE(all.equal(oc$p, n1$p[m], tolerance = 1e-12))
  f_new <- sapply(stats, function(st) decoy_fdr(calls, st, "coex", B = 2L, decoy_filter = calls$kind == "target" | calls$set %in% 1L)$fdr)
  f_old <- sapply(stats, function(st) decoy_fdr(oc, st, "coex", B = 2L)$fdr)
  cat(sprintf("SET-1 CHECK against round 6: decoys identical %s; calls and P identical %s; FDR %s vs %s\n", g_ok, c_ok,
              paste(sprintf("%.4f", f_new), collapse = "/"), paste(sprintf("%.4f", f_old), collapse = "/")))
  invisible(g_ok && c_ok && isTRUE(all.equal(unname(f_new), unname(f_old), tolerance = 1e-12)))
}
save_all <- function(prefix, calls, info, clusters, stats) {
  write.csv(calls, gzfile(file.path(OUT, sprintf("%s_calls.csv.gz", prefix))), row.names = FALSE)
  F <- fdr_table(calls, info, clusters, stats); print(F[, c("stat", "decoys", "boot", "called", "fdr", "fdr_lo", "fdr_hi",
                                                            "fdr_set_min", "fdr_set_max", "pi0", "fdr_pi0", "fdr_pi0_lo", "fdr_pi0_hi")])
  write.csv(F, file.path(OUT, sprintf("%s_fdr.csv", prefix)), row.names = FALSE)
  write.csv(by_set(calls, stats), file.path(OUT, sprintf("%s_fdr_by_set.csv", prefix)), row.names = FALSE)
  write.csv(pi0_sens(calls, stats, clusters), file.path(OUT, sprintf("%s_pi0_lambda.csv", prefix)), row.names = FALSE)
  for (st in stats) write.csv(target_vs_decoys(calls, st, "coex", same_direction = st == "dir"),
                              file.path(OUT, sprintf("%s_%s_target_vs_decoys.csv", prefix, st)), row.names = FALSE)
}

if (app == "response") {
  co <- response_cohorts(); L <- immune_library(key = "ensembl")
  res <- run_decoy_fdr(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, n_sets = K, n_draw = 1000, keep_draws = 0,
                       seed = 20261022L, cores = NC, B = 6000L)
  cat("fit", mins(), "\n"); print(res$decoys)
  write.csv(res$decoys$info, gzfile(file.path(OUT, "response_decoy_info.csv.gz")), row.names = FALSE)
  write.csv(res$fit$units, gzfile(file.path(OUT, "response_units.csv.gz")), row.names = FALSE)
  check_set1(res$calls, "response_calls.csv", "response_decoy_info.csv", res$decoys$info, c("fold", "dir"))
  save_all("response", res$calls, res$decoys$info, imm_clusters, c("fold", "dir"))
}

if (app == "tcga") {
  L <- immune_library(key = "symbol")
  co <- tcga_cohorts("OS")
  dec <- grow_decoys(co, L$lib, exclude = L$exclude, n_coex = 32, n_rand = 8, n_sets = K, seed = 20261022L, cores = NC,
                     max_bytes = 200e9)
  print(dec); saveRDS(dec, file.path(OUT, "tcga_decoys.rds"))
  write.csv(dec$info, gzfile(file.path(OUT, "tcga_decoy_info.csv.gz")), row.names = FALSE)
  cat("decoys", mins(), "\n")
  fit <- fit_units(co, L$lib, dec, n_draw = 1000L, keep_draws = 200L, seed = 20261023L, cores = NC)
  saveRDS(fit, file.path(OUT, "tcga_OS_fit.rds"))
  write.csv(fit$units, gzfile(file.path(OUT, "tcga_OS_units.csv.gz")), row.names = FALSE)
  calls <- call_units(fit$units)
  check_set1(calls, "tcga_OS_calls.csv", "tcga_decoy_info.csv", dec$info, c("fold", "dir"))
  save_all("tcga_OS", calls, dec$info, imm_clusters, c("fold", "dir"))
  cat("OS fit", mins(), "\n")
  cal <- calibrate_null(fit, B = 200L, seed = 20261024L, cores = NC)
  print(cal$summary)
  for (nm in names(cal)) write.csv(cal[[nm]], file.path(OUT, sprintf("tcga_OS_globalnull_%s.csv", nm)), row.names = FALSE)
  cat("calibration", mins(), "\n")
  rm(fit, cal); gc()
  for (ep in c("DSS", "PFI")) {
    ce <- tcga_cohorts(ep)
    fe <- fit_units(ce, L$lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
    write.csv(fe$units, gzfile(file.path(OUT, sprintf("tcga_%s_units.csv.gz", ep))), row.names = FALSE)
    ca <- call_units(fe$units)
    check_set1(ca, sprintf("tcga_%s_calls.csv", ep), "tcga_decoy_info.csv", dec$info, c("fold", "dir"))
    cat(ep, "\n"); save_all(sprintf("tcga_%s", ep), ca, dec$info, imm_clusters, c("fold", "dir"))
    rm(fe); gc()
  }
}

if (app == "tcga_ep") {             # one secondary endpoint with the saved decoys (the same fit as in app tcga)
  ep <- args[3]; stopifnot(ep %in% c("DSS", "PFI"))
  L <- immune_library(key = "symbol")
  dec <- readRDS(file.path(OUT, "tcga_decoys.rds"))
  ce <- tcga_cohorts(ep)
  fe <- fit_units(ce, L$lib, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
  write.csv(fe$units, gzfile(file.path(OUT, sprintf("tcga_%s_units.csv.gz", ep))), row.names = FALSE)
  ca <- call_units(fe$units)
  check_set1(ca, sprintf("tcga_%s_calls.csv", ep), "tcga_decoy_info.csv", dec$info, c("fold", "dir"))
  cat(ep, "\n"); save_all(sprintf("tcga_%s", ep), ca, dec$info, imm_clusters, c("fold", "dir"))
  app <- paste(app, ep)
}

if (app %in% c("breast", "breast_abs")) {
  co <- breast_cohorts(); sets <- breast_sets()
  nmeas <- sapply(sets, function(g) min(vapply(co, function(c) sum(g %in% rownames(c$Z)), 1L)))
  keep <- names(nmeas)[nmeas >= 3]
  lib <- dfdr_library(sets[keep], engine = "pc1"); print(lib)
  absolute <- app == "breast_abs"
  dec <- grow_decoys(co, lib, exclude = unique(unlist(sets)), n_coex = 32, n_rand = 8, n_sets = K,
                     seed = if (absolute) 20261025L else 20261022L, cores = NC, coherence = if (absolute) "absolute" else "signed")
  print(dec); saveRDS(dec, file.path(OUT, sprintf("%s_decoys.rds", app)))
  write.csv(dec$info, gzfile(file.path(OUT, sprintf("%s_decoy_info.csv.gz", app))), row.names = FALSE)
  cat("decoys", mins(), "\n")
  fit <- fit_units(co, lib, dec, n_draw = 1000L, keep_draws = if (absolute) 0L else 200L, seed = 20261023L, cores = NC)
  write.csv(fit$units, gzfile(file.path(OUT, sprintf("%s_pc1_units.csv.gz", app))), row.names = FALSE)
  # breast clusters (plan, analysis 6): average linkage on 1 - |median over cohorts of the Spearman r between the
  # targets' first-principal-component scores|, cut at 0.5; written once, by the primary run
  cf <- file.path(OUT, "breast_clusters.csv")
  if (!absolute) {
    R <- simplify2array(lapply(fit$scores, function(s) { r <- matrix(NA, length(keep), length(keep), dimnames = list(keep, keep))
      x <- s$targets[intersect(keep, rownames(s$targets)), , drop = FALSE]; r[rownames(x), rownames(x)] <- cor(t(x), method = "spearman"); r }))
    Rm <- apply(R, c(1, 2), median, na.rm = TRUE)
    hc <- stats::hclust(stats::as.dist(1 - abs(Rm)), method = "average")
    cl <- stats::cutree(hc, h = 0.5)
    write.csv(data.frame(signature = names(cl), cluster = cl), cf, row.names = FALSE)
    cat(length(unique(cl)), "breast clusters\n")
    saveRDS(fit, file.path(OUT, "breast_pc1_fit.rds"))
  }
  CLb <- read.csv(cf, stringsAsFactors = FALSE); bclusters <- setNames(CLb$cluster, CLb$signature)
  calls <- call_units(fit$units, stats = "fold")
  check_set1(calls, if (absolute) "breast_abs_calls.csv" else "breast_pc1_calls.csv",
             if (absolute) "breast_abs_decoy_info.csv" else "breast_decoy_info.csv", dec$info, "fold")
  save_all(sprintf("%s_pc1", app), calls, dec$info, bclusters, "fold")
  cat("pc1 fit", mins(), "\n")
  if (!absolute) {
    cal <- calibrate_null(fit, B = 200L, stats = "fold", seed = 20261024L, cores = NC)
    print(cal$summary)
    for (nm in names(cal)) write.csv(cal[[nm]], file.path(OUT, sprintf("breast_pc1_globalnull_%s.csv", nm)), row.names = FALSE)
    rm(fit, cal); gc()
    libm <- dfdr_library(sets[keep], engine = "mean")
    fm <- fit_units(co, libm, dec, n_draw = 1000L, keep_draws = 0L, seed = 20261023L, cores = NC)
    write.csv(fm$units, gzfile(file.path(OUT, "breast_mean_units.csv.gz")), row.names = FALSE)
    cm <- call_units(fm$units, stats = "fold")
    check_set1(cm, "breast_mean_calls.csv", "breast_decoy_info.csv", dec$info, "fold")
    save_all("breast_mean", cm, dec$info, bclusters, "fold")
  }
}
cat(app, "all done in", mins(), "\n")
