#' The whole procedure in one call
#'
#' Grows decoys ([grow_decoys()]), scores targets, decoys and null banks ([fit_units()]), calls every
#' unit ([call_units()]) and estimates the FDR of the library claim for each statistic ([decoy_fdr()]).
#'
#' @param cohorts list of `dfdr_cohort` objects.
#' @param lib a `dfdr_library`.
#' @param n_coex,n_rand,n_sets,tol,exclude passed to [grow_decoys()]; the default of 10 decoy sets is the
#'   paper's.
#' @param n_draw,keep_draws passed to [fit_units()].
#' @param stats statistics; `"dir"` needs declared directions in `lib`.
#' @param seed one seed for decoys, nulls and bootstrap.
#' @param cores forked workers.
#' @param B bootstrap replicates for the FDR interval.
#' @return A list of class `dfdr_result`: `decoys`, `fit`, `calls` and `fdr` (one row per statistic and
#'   decoy class).
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy()
#' res <- run_decoy_fdr(toy$cohorts, toy$lib, n_coex = 4, n_rand = 1, n_sets = 2, n_draw = 50, seed = 3, B = 100)
#' res$fdr[, c("stat", "class", "called", "fdr")]
#' @export
run_decoy_fdr <- function(cohorts, lib, n_coex = 32L, n_rand = 8L, n_sets = 10L, tol = 0.02, exclude = character(0),
                          n_draw = 1000L, keep_draws = 0L, stats = c("fold", "dir"), seed = 1L, cores = 1L,
                          B = 6000L) {
  if ("dir" %in% stats && all(is.na(lib$direction))) stats <- setdiff(stats, "dir")
  dec <- grow_decoys(cohorts, lib, exclude = exclude, n_coex = n_coex, n_rand = n_rand, n_sets = n_sets, tol = tol,
                     seed = seed, cores = cores)
  fit <- fit_units(cohorts, lib, dec, n_draw = n_draw, keep_draws = keep_draws, seed = seed + 1L, cores = cores)
  calls <- call_units(fit$units, stats = stats)
  fdr <- do.call(rbind, lapply(stats, function(st) do.call(rbind, lapply(intersect(c("coex", "rand"), dec$info$class),
    function(cl) decoy_fdr(calls, st, cl, B = B, seed = seed)))))
  structure(list(decoys = dec, fit = fit, calls = calls, fdr = fdr), class = "dfdr_result")
}

#' @export
print.dfdr_result <- function(x, ...) {
  cat("decoyfdr result\n"); print(x$decoys)
  f <- x$fdr
  for (i in seq_len(nrow(f)))
    cat(sprintf("  %-4s %-4s decoys: %2d of %d called, decoy call rate %.1f%%, expected false %.1f, FDR %.1f%% (%.1f-%.1f); pi0 %.2f, FDR %.1f%%\n",
                f$stat[i], f$class[i], as.integer(f$called[i]), f$m[i], f$decoy_rate[i], f$exp_false[i], f$fdr[i],
                f$fdr_lo[i], f$fdr_hi[i], f$pi0[i], f$fdr_pi0[i]))
  invisible(x)
}
