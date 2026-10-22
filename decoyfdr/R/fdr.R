#' Target-decoy false discovery rate of a library claim
#'
#' The claim is "these R of the m signatures are called". Each target's decoys estimate the probability
#' that a module with the target's structure and no signature gene is called; the expected number of
#' false calls is the sum of those probabilities over targets, estimated as m times their mean, and the
#' FDR is that number over R. Targets without an eligible decoy drop out of both terms. The sum runs over
#' every target, true ones included, so this is the estimate with pi0 = 1. The pi0-corrected estimate
#' multiplies it by the share of targets estimated to be null. Each target's statistic, the sign of its
#' median delta times -log10 of its P value, is ranked among its eligible decoys': q is one plus the number
#' of decoys at or above it, over one plus the number of decoys; pi0 is the number of targets with q above
#' `lambda` over m (1 - lambda), capped at 1 (Storey 2002, with the decoys as the null reference instead of
#' a uniform P value).
#'
#' @param calls output of [call_units()].
#' @param stat `"fold"` or `"dir"`.
#' @param decoy_class decoy class to use (`"coex"`, `"rand"`, or the class labels of a custom table).
#' @param same_direction for the directed statistic, use only decoys declared in their target's
#'   direction (as in the paper). Ignored for the folded statistic.
#' @param B bootstrap replicates.
#' @param clusters optional named vector mapping targets to clusters; the bootstrap then resamples
#'   clusters, each with its targets and their decoys, instead of targets.
#' @param seed bootstrap seed.
#' @param decoy_filter optional logical vector over the rows of `calls` restricting the eligible decoys
#'   further (e.g. the closest-matched ones).
#' @param lambda the pi0 threshold.
#' @return A one-row data frame: `stat, class, m, called, decoy_rate, exp_false, exp_true, fdr, fdr_lo,
#'   fdr_hi, n_decoys, n_sets, fdr_set_min, fdr_set_max, pi0, fdr_pi0, fdr_pi0_lo, fdr_pi0_hi` (rates and
#'   FDRs in percent). With several decoy sets the estimate pools them, which equals the mean of the
#'   per-set estimates; `fdr_set_min` and `fdr_set_max` give their range.
#' @examples
#' U <- read.csv(system.file("extdata", "response_units_N0.csv.gz", package = "decoyfdr"))
#' calls <- call_units(U)
#' decoy_fdr(calls, "fold", "coexA", B = 200)
#' @export
decoy_fdr <- function(calls, stat = "fold", decoy_class = "coex", same_direction = TRUE, B = 6000L,
                      clusters = NULL, seed = 20261019L, decoy_filter = NULL, lambda = 0.5) {
  C <- calls[calls$stat == stat, , drop = FALSE]
  if (!is.null(decoy_filter)) decoy_filter <- decoy_filter[calls$stat == stat]
  tt <- C[C$kind == "target", , drop = FALSE]
  dsel <- C$kind == "decoy" & C$class == decoy_class
  if (!is.null(decoy_filter)) dsel <- dsel & decoy_filter
  dd <- C[dsel, , drop = FALSE]
  if (stat == "dir" && same_direction) dd <- dd[dd$dirn == tt$dirn[match(dd$target, tt$target)], , drop = FALSE]
  if (is.null(dd$set)) dd$set <- rep(1L, nrow(dd))
  pit <- tapply(dd$call, dd$target, mean)[tt$target]
  call <- as.numeric(tt$call)
  sstat <- function(x) sign(x$med_delta) * -log10(pmax(x$p, 1e-300))
  st_t <- sstat(tt); st_d <- split(sstat(dd), dd$target)
  q <- vapply(seq_len(nrow(tt)), function(i) { x <- st_d[[tt$target[i]]]
    if (is.null(x)) NA_real_ else (1 + sum(x >= st_t[i] - 1e-12)) / (1 + length(x)) }, 1.0)
  est <- function(ix) {
    ok <- !is.na(pit[ix]); Rn <- sum(call[ix][ok])
    E <- if (any(ok)) mean(pit[ix][ok]) * sum(ok) else NA_real_
    p0 <- if (any(ok)) min(1, sum(q[ix][ok] > lambda) / (sum(ok) * (1 - lambda))) else NA_real_
    c(R = Rn, pi = if (any(ok)) mean(pit[ix][ok]) else NA_real_, E = E, F = if (Rn > 0) E / Rn else NA_real_,
      pi0 = p0, F0 = if (Rn > 0) p0 * E / Rn else NA_real_)
  }
  pt <- est(seq_along(call))
  sets <- sort(unique(dd$set))
  fset <- vapply(sets, function(s) { pk <- tapply(dd$call[dd$set == s], dd$target[dd$set == s], mean)[tt$target]
    ok <- !is.na(pk); Rn <- sum(call[ok]); if (Rn > 0) mean(pk[ok]) * sum(ok) / Rn else NA_real_ }, 1.0)
  set.seed(seed)
  if (is.null(clusters)) {
    bs <- replicate(B, est(sample.int(length(call), replace = TRUE)))
  } else {
    cl <- clusters[tt$target]; if (anyNA(cl)) stop("every target needs a cluster")
    groups <- split(seq_along(call), cl)
    bs <- replicate(B, est(unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE)))
  }
  q95 <- apply(bs, 1, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(stat = stat, class = decoy_class, m = sum(!is.na(pit)), called = pt[["R"]],
             decoy_rate = 100 * pt[["pi"]], exp_false = pt[["E"]], exp_true = pt[["R"]] - pt[["E"]],
             fdr = 100 * pt[["F"]], fdr_lo = 100 * q95[1, "F"], fdr_hi = 100 * q95[2, "F"], n_decoys = nrow(dd),
             n_sets = length(sets), fdr_set_min = 100 * min(fset, na.rm = TRUE), fdr_set_max = 100 * max(fset, na.rm = TRUE),
             pi0 = pt[["pi0"]], fdr_pi0 = 100 * pt[["F0"]], fdr_pi0_lo = 100 * q95[1, "F0"], fdr_pi0_hi = 100 * q95[2, "F0"],
             row.names = NULL)
}

#' Where each target sits among its own decoys
#'
#' @inheritParams decoy_fdr
#' @return One row per target: its median delta, its Wilcoxon P and call, the number of its decoys with
#'   a median delta at least as large, and the decoy-calibrated P `(1 + that number) / (1 + decoys)`.
#' @export
target_vs_decoys <- function(calls, stat = "fold", decoy_class = "coex", same_direction = FALSE) {
  C <- calls[calls$stat == stat, , drop = FALSE]
  tt <- C[C$kind == "target", , drop = FALSE]
  dd <- C[C$kind == "decoy" & C$class == decoy_class, , drop = FALSE]
  if (stat == "dir" && same_direction) dd <- dd[dd$dirn == tt$dirn[match(dd$target, tt$target)], , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(tt)), function(i) {
    x <- dd$med_delta[dd$target == tt$target[i]]
    data.frame(target = tt$target[i], med_delta = tt$med_delta[i], p = tt$p[i], call = tt$call[i],
               decoys = length(x), decoys_at_or_above = sum(x >= tt$med_delta[i]),
               p_decoy = (1 + sum(x >= tt$med_delta[i])) / (1 + length(x)))
  }))
}
