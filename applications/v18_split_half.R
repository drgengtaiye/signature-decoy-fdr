#!/usr/bin/env Rscript
# Round 12 (validation/PLAN_round12_reproducibility_by_type.md), sections A and B: every cohort split at random into
# halves stratified by event status, 500 times; in each half the folded C of every target and co-expressed decoy from
# per-patient scores, and rho (the share of the target's co-expressed decoys below it, ties one half).
#   DFDR_LIB=<frozen install> Rscript v18_split_half.R tcga|external [cores]
#   tcga:     scores stored by the round-7 fit (icb/v12_pkg/tcga_OS_fit.rds); splits seed 20261105
#   external: scores recomputed from the round-11 decoys by the steps of fit_units, the null-bank draws made but not
#             scored, so the weighted-mean decoys get the weights of round 11 (reused from external_scores.rds when
#             present); splits seed 20261106
# AMENDMENT_1_round12: folded C values within 1e-12 of each other count as tied (folding can leave mathematically equal
# values one unit apart in the last place), and splits run in blocks of 25 per cohort, in parallel.
# AMENDMENT_3_round12 (post hoc sensitivity): apps tcga_boot and external_boot, 200 nonparametric bootstrap resamples
# of each cohort's patients (seeds 20261107 and 20261108), rho on each resample.
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
args <- commandArgs(TRUE); app <- args[1]; NC <- if (length(args) > 1) as.integer(args[2]) else 8L
OUT <- file.path(G, "icb", "v18_split"); dir.create(OUT, showWarnings = FALSE)
NS <- 500L; NB <- 200L; BLOCK <- 25L; TIE <- 1e-12; t0 <- Sys.time()
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# halves A of every split, drawn serially in cohort order from one stream
draw_splits <- function(cohorts, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  lapply(cohorts, function(co) lapply(seq_len(NS), function(b) {
    a <- integer(0)
    for (e in c(1L, 0L)) { i <- which(co$event == e); i <- i[sample.int(length(i))]; a <- c(a, i[seq_len(length(i) %/% 2L)]) }
    sort(a)
  }))
}

# bootstrap resamples of every cohort, drawn serially in cohort order from one stream
draw_boot <- function(cohorts, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  lapply(cohorts, function(co) lapply(seq_len(NB), function(b) sample.int(co$n, co$n, replace = TRUE)))
}

# rho of every target from folded C of the targets (rows of St) and their co-expressed decoys (rows of Sd)
rho_of <- function(St, Sd, dtarget, time, event) {
  ft <- fold(cindex_rows(St, time, event)); fd <- fold(cindex_rows(Sd, time, event))
  names(ft) <- rownames(St)
  vapply(rownames(St), function(s) { x <- fd[dtarget == s] - ft[s]; mean((x < -TIE) + 0.5 * (abs(x) <= TIE)) }, 1)
}

split_cohort <- function(co, St, Sd, dtarget, A, idx) {
  rows <- list()
  if (0L %in% idx) {
    full <- rho_of(St, Sd, dtarget, co$time, co$event)
    rows[[1]] <- data.frame(cohort = co$name, split = 0L, half = "full", target = names(full), rho = unname(full))
  }
  for (b in setdiff(idx, 0L)) {
    a <- A[[b]]; bb <- setdiff(seq_len(co$n), a)
    ra <- rho_of(St[, a, drop = FALSE], Sd[, a, drop = FALSE], dtarget, co$time[a], co$event[a])
    rb <- rho_of(St[, bb, drop = FALSE], Sd[, bb, drop = FALSE], dtarget, co$time[bb], co$event[bb])
    rows[[length(rows) + 1L]] <- data.frame(cohort = co$name, split = b, half = rep(c("A", "B"), each = length(ra)),
                                            target = c(names(ra), names(rb)), rho = c(unname(ra), unname(rb)))
  }
  do.call(rbind, rows)
}

# folded C of every target and decoy on the whole cohort, against the stored units
check_full <- function(co, St, Sd, units) {
  f <- c(fold(cindex_rows(St, co$time, co$event)), fold(cindex_rows(Sd, co$time, co$event)))
  names(f) <- c(rownames(St), rownames(Sd))
  u <- units[units$cohort == co$name, ]
  m <- match(names(f), u$unit)
  data.frame(cohort = co$name, rows = length(f), unmatched = sum(is.na(m)),
             max_abs_diff = max(abs(f - u$fold_obs[m]), na.rm = TRUE))
}

run_splits <- function(cohorts, scores, dec, units, seed, prefix) {
  A <- draw_splits(cohorts, seed)
  DEC <- dec$decoys
  # tasks: (cohort, block of splits); block 0 also holds the whole cohort and the check; largest cohorts first
  blocks <- split(0:NS, c(0L, (seq_len(NS) - 1L) %/% BLOCK))
  ord <- order(-vapply(cohorts, `[[`, 1L, "n"))
  tasks <- do.call(rbind, lapply(ord, function(i) data.frame(i = i, b = seq_along(blocks))))
  one <- function(k) {
    i <- tasks$i[k]; idx <- blocks[[tasks$b[k]]]
    co <- cohorts[[i]]; sc <- scores[[i]]
    cx <- rownames(sc$decoys)[vapply(DEC[rownames(sc$decoys)], `[[`, "", "class") == "coex"]
    Sd <- sc$decoys[cx, , drop = FALSE]; dtarget <- vapply(DEC[cx], `[[`, "", "target")
    St <- sc$targets[rownames(sc$targets) %in% unique(dtarget), , drop = FALSE]
    list(check = if (0L %in% idx) check_full(co, sc$targets, sc$decoys, units) else NULL,
         rho = split_cohort(co, St, Sd, dtarget, A[[i]], idx))
  }
  res <- parallel::mclapply(seq_len(nrow(tasks)), one, mc.cores = NC, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
  if (any(bad)) stop("failed: ", paste(unique(vapply(cohorts[tasks$i[bad]], `[[`, "", "name")), collapse = ", "),
                     if (inherits(res[[which(bad)[1]]], "try-error")) as.character(res[[which(bad)[1]]]) else "")
  ck <- do.call(rbind, lapply(res, `[[`, "check")); ck <- ck[order(ck$cohort), ]; print(ck)
  write.csv(ck, file.path(OUT, sprintf("%s_check_full.csv", prefix)), row.names = FALSE)
  R <- do.call(rbind, lapply(res, `[[`, "rho")); R <- R[order(R$cohort, R$split, R$half), ]
  write.csv(R, gzfile(file.path(OUT, sprintf("rho_split_%s.csv.gz", prefix))), row.names = FALSE)
  sz <- data.frame(cohort = vapply(cohorts, `[[`, "", "name"), n = vapply(cohorts, `[[`, 1L, "n"),
                   events = vapply(cohorts, function(c) sum(c$event), 1),
                   half_a = vapply(A, function(a) length(a[[1]]), 1L))
  write.csv(sz, file.path(OUT, sprintf("%s_sizes.csv", prefix)), row.names = FALSE)
  cat(prefix, "splits", mins(), "\n")
}

run_boot <- function(cohorts, scores, dec, seed, prefix) {
  Bx <- draw_boot(cohorts, seed)
  DEC <- dec$decoys
  blocks <- split(seq_len(NB), (seq_len(NB) - 1L) %/% BLOCK)
  ord <- order(-vapply(cohorts, `[[`, 1L, "n"))
  tasks <- do.call(rbind, lapply(ord, function(i) data.frame(i = i, b = seq_along(blocks))))
  one <- function(k) {
    i <- tasks$i[k]; co <- cohorts[[i]]; sc <- scores[[i]]
    cx <- rownames(sc$decoys)[vapply(DEC[rownames(sc$decoys)], `[[`, "", "class") == "coex"]
    Sd <- sc$decoys[cx, , drop = FALSE]; dtarget <- vapply(DEC[cx], `[[`, "", "target")
    St <- sc$targets[rownames(sc$targets) %in% unique(dtarget), , drop = FALSE]
    do.call(rbind, lapply(blocks[[tasks$b[k]]], function(b) {
      x <- Bx[[i]][[b]]
      r <- rho_of(St[, x, drop = FALSE], Sd[, x, drop = FALSE], dtarget, co$time[x], co$event[x])
      data.frame(cohort = co$name, boot = b, target = names(r), rho = unname(r))
    }))
  }
  res <- parallel::mclapply(seq_len(nrow(tasks)), one, mc.cores = NC, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
  if (any(bad)) stop("failed: ", paste(unique(vapply(cohorts[tasks$i[bad]], `[[`, "", "name")), collapse = ", "))
  R <- do.call(rbind, res); R <- R[order(R$cohort, R$boot), ]
  write.csv(R, gzfile(file.path(OUT, sprintf("rho_boot_%s.csv.gz", prefix))), row.names = FALSE)
  cat(prefix, "bootstrap", mins(), "\n")
}

if (app == "tcga_boot") {
  fit <- readRDS(file.path(D, "v12_pkg", "tcga_OS_fit.rds"))
  run_boot(fit$cohorts, fit$scores, fit$decoys, 20261107L, "tcga")
}

if (app == "tcga") {
  fit <- readRDS(file.path(D, "v12_pkg", "tcga_OS_fit.rds"))
  cat("loaded", mins(), "\n")
  units <- read.csv(gzfile(file.path(D, "v12_pkg", "tcga_OS_units.csv.gz")), stringsAsFactors = FALSE)
  run_splits(fit$cohorts, fit$scores, fit$decoys, units, 20261105L, "tcga")
}

if (app %in% c("external", "external_boot")) {
  # cohorts and library exactly as applications/v16_external.R (survival)
  L <- immune_library(key = "symbol"); EXT <- file.path(G, "external", "cohorts")
  man <- read.csv(file.path(EXT, "manifest.csv"), stringsAsFactors = FALSE); man <- man[as.logical(man$included), ]
  co <- lapply(man$cohort, function(id) {
    e <- read.delim(gzfile(file.path(EXT, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE)
    M <- as.matrix(e); M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
    cl <- read.delim(file.path(EXT, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
    dfdr_cohort(M[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample),
                name = id)
  })
  lb <- L$lib
  nmeas <- sapply(lb$sets, function(g) min(vapply(co, function(c) sum(g %in% rownames(c$Z)), 1L)))
  keep <- names(nmeas)[nmeas >= 3]; wm <- keep[lb$engine[keep] == "wmean"]
  lib <- dfdr_library(lb$sets[keep], lb$engine[keep], lb$weights[wm], lb$direction[keep]); print(lib)
  dec <- readRDS(file.path(G, "icb", "v16_external", "external_decoys.rds"))
  units <- read.csv(gzfile(file.path(G, "icb", "v16_external", "external_surv_units.csv.gz")), stringsAsFactors = FALSE)
  # the steps of decoyfdr::fit_units (seed 20261102, n_draw 1000), bank draws made and not scored
  n_draw <- 1000L
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261102L)
  streams <- vector("list", length(co)); s <- .Random.seed
  for (i in seq_along(co)) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }
  DEC <- dec$decoys
  loadNamespace("GSVA")
  score_one <- function(ci) {
    assign(".Random.seed", streams[[ci]], envir = globalenv())
    cc <- co[[ci]]; present <- rownames(cc$Z)
    tg <- names(lib$sets)
    tset <- lapply(lib$sets, function(g) g[g %in% present])
    tw <- lapply(tg, function(s) if (lib$engine[s] == "wmean") lib$weights[[s]][lib$sets[[s]] %in% present] else NULL)
    names(tw) <- tg
    tg <- tg[vapply(tset[tg], length, 1L) >= 3L]
    obs_t <- score_library(cc, dfdr_library(tset[tg], lib$engine[tg], tw[tg], lib$direction[tg]))
    for (s in tg) {
      k <- length(tset[[s]]); eng <- lib$engine[[s]]
      sets <- replicate(n_draw, present[sample.int(length(present), k)], simplify = FALSE)
      if (eng == "wmean") w <- lapply(seq_len(n_draw), function(i) sample(tw[[s]]))
    }
    dn <- names(DEC)[vapply(DEC, function(d) d$target %in% tg, TRUE)]
    obs_d <- matrix(NA_real_, length(dn), cc$n, dimnames = list(dn, cc$samples))
    for (eng in unique(lib$engine[tg])) {
      nm <- dn[vapply(DEC[dn], function(d) lib$engine[[d$target]] == eng, TRUE)]
      if (!length(nm)) next
      sets <- lapply(DEC[nm], function(d) d$genes[d$genes %in% present])
      w <- if (eng == "wmean") lapply(nm, function(n_) { ws <- tw[[DEC[[n_]]$target]]; L <- length(sets[[n_]])
        sample(ws, L, replace = length(ws) != L) }) else NULL
      if (!is.null(w)) names(w) <- nm
      obs_d[nm, ] <- score_sets(cc, sets, eng, w)
    }
    list(targets = obs_t, decoys = obs_d)
  }
  sf <- file.path(OUT, "external_scores.rds")
  if (file.exists(sf)) { scores <- readRDS(sf); cat("scores reused from", sf, "\n") } else {
    invisible(gc())
    scores <- parallel::mclapply(seq_along(co), score_one, mc.cores = NC, mc.preschedule = FALSE)
    bad <- vapply(scores, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
    if (any(bad)) stop("scoring failed: ", paste(man$cohort[bad], collapse = ", "))
    cat("scored", mins(), "\n")
    saveRDS(scores, sf)
  }
  if (app == "external") run_splits(co, scores, dec, units, 20261106L, "external") else
    run_boot(co, scores, dec, 20261108L, "external")
}
