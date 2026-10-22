#' Calibrate the estimate under the global null
#'
#' Permutes the outcome within each cohort `B` times, which makes every target null and every call
#' false, and reruns the whole procedure on each permutation from the scores kept by [fit_units()]
#' (`keep_draws > 0`): signed statistics of targets, decoys and the kept bank draws, bank medians, deltas,
#' the signed-rank calls and the decoy estimate. The ratio of the mean estimated to the mean realized
#' number of false calls is 1 for a calibrated estimate.
#'
#' @param fit a `dfdr_fit` with kept scores.
#' @param B number of permutations.
#' @param stats statistics, any of `"fold"` and `"dir"`.
#' @param decoy_class decoy class for the estimate.
#' @param alpha,min_cohorts the call rule of [call_units()].
#' @param seed permutation seed (one L'Ecuyer-CMRG stream per cohort).
#' @param cores forked workers over cohorts.
#' @param nboot bootstrap replicates (over permutations) for the interval of the ratio.
#' @return A list: `summary` (per statistic: mean realized and estimated false calls, their ratio with a
#'   95% interval, and the per-unit call rates of targets and of each decoy class), `per_perm` (R and E for
#'   every permutation) and `per_target` (each target's call rate and its decoys' mean call rate).
#' @export
calibrate_null <- function(fit, B = 200L, stats = c("fold", "dir"), decoy_class = "coex", alpha = 0.05,
                           min_cohorts = 6L, seed = 1L, cores = 1L, nboot = 2000L) {
  if (any(vapply(fit$scores, is.null, TRUE))) stop("fit_units() must be run with keep_draws > 0")
  lib <- fit$lib; DEC <- fit$decoys$decoys
  targets <- names(lib$sets); decoys <- names(DEC)
  dtarget <- vapply(DEC, `[[`, "", "target"); ddirn <- vapply(DEC, `[[`, 1L, "dirn")
  dclass <- vapply(DEC, `[[`, "", "class")
  tdirn <- ifelse(is.na(lib$direction), 1L, lib$direction)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  streams <- vector("list", length(fit$cohorts)); s <- .Random.seed
  for (i in seq_along(fit$cohorts)) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }

  one <- function(ci) {
    assign(".Random.seed", streams[[ci]], envir = globalenv())
    co <- fit$cohorts[[ci]]
    perm_deltas(co, fit$scores[[co$name]], replicate(B, sample.int(co$n)), tdirn, ddirn, dtarget)
  }
  if (cores > 1L) invisible(gc())   # a forked worker copies whatever garbage its own collections touch
  res <- if (cores > 1L) parallel::mclapply(seq_along(fit$cohorts), one, mc.cores = cores, mc.preschedule = FALSE) else
    lapply(seq_along(fit$cohorts), one)
  err <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)   # NULL: a forked worker died (often memory)
  if (any(err)) stop("calibration failed: ", if (is.null(res[[which(err)[1]]]))
    "a worker returned no result (out of memory? use fewer cores)" else as.character(res[[which(err)[1]]]))

  nc <- length(res)
  arr <- function(units, key, which) {
    X <- array(NA_real_, c(nc, length(units), B))
    for (ci in seq_len(nc)) { r <- res[[ci]]; nm <- r[[which]]; X[ci, match(nm, units), ] <- r[[key]] }
    X
  }
  callfun <- function(X) {
    U <- dim(X)[2]; out <- matrix(FALSE, U, B)
    for (u in seq_len(U)) {
      M <- X[, u, , drop = TRUE]; if (is.null(dim(M))) M <- matrix(M, nrow = nc)
      ok <- rowSums(is.finite(M)) == B
      if (sum(ok) < min_cohorts) next
      M <- M[ok, , drop = FALSE]
      p <- signed_rank_p(M); m <- matrixStats::colMedians(M)
      p[colSums(M != 0) == 0] <- 1
      out[u, ] <- p < alpha & m > 0
    }
    out
  }
  summ <- list(); perb <- list(); pert <- list()
  for (st in stats) {
    kt <- if (st == "fold") "tf" else "td"; kd <- if (st == "fold") "df" else "dd"
    cT <- callfun(arr(targets, kt, "tg")); cD <- callfun(arr(decoys, kd, "dn"))
    elig <- dclass == decoy_class
    if (st == "dir") elig <- elig & ddirn == tdirn[dtarget]
    num <- rowsum(cD[elig, , drop = FALSE] * 1, dtarget[elig])
    den <- as.numeric(table(dtarget[elig])[rownames(num)])
    pit <- num / den
    okt <- targets %in% rownames(pit)
    Rb <- colSums(cT[okt, , drop = FALSE]); Eb <- colMeans(pit) * sum(okt)
    rho <- mean(Eb) / mean(Rb)
    bs <- replicate(nboot, { i <- sample.int(B, replace = TRUE); mean(Eb[i]) / mean(Rb[i]) })
    rate <- function(m) 100 * mean(m)
    summ[[st]] <- data.frame(stat = st, perms = B, realized = mean(Rb), estimated = mean(Eb), ratio = rho,
      ratio_lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)), ratio_hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)),
      rate_targets = rate(cT[okt, ]), rate_coex = rate(cD[dclass == "coex", ]),
      rate_rand = if (any(dclass == "rand")) rate(cD[dclass == "rand", ]) else NA_real_)
    perb[[st]] <- data.frame(stat = st, perm = seq_len(B), R = Rb, E = Eb)
    pert[[st]] <- data.frame(stat = st, target = targets[okt], rate = rowMeans(cT[okt, , drop = FALSE]),
                             decoy_rate = rowMeans(pit[targets[okt], , drop = FALSE]))
  }
  list(summary = do.call(rbind, summ), per_perm = do.call(rbind, perb), per_target = do.call(rbind, pert))
}

# Folded and directed deltas of every target and decoy of one cohort under each column of perms (a matrix of
# permutations of the samples), read against the medians of the kept bank draws under the same permutation.
perm_deltas <- function(co, sc, perms, tdirn, ddirn, dtarget) {
  tg <- rownames(sc$targets); dn <- rownames(sc$decoys)
  kb <- vapply(sc$bank[tg], nrow, 1L)
  stack <- rbind(sc$targets, sc$decoys, do.call(rbind, sc$bank[tg]))
  B <- ncol(perms)
  A <- if (co$type == "binary") {
    R <- matrixStats::rowRanks(stack, ties.method = "average")
    Y <- matrix(co$y[perms], co$n, B)
    n1 <- sum(co$y); n0 <- co$n - n1
    (R %*% Y - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else {
    vapply(seq_len(B), function(b) cindex_rows(stack, co$time[perms[, b]], co$event[perms[, b]]),
           numeric(nrow(stack)))
  }
  A <- matrix(A, nrow(stack), B)
  nt <- length(tg); nd <- length(dn)
  aT <- A[seq_len(nt), , drop = FALSE]; aD <- A[nt + seq_len(nd), , drop = FALSE]
  off <- nt + nd + c(0, cumsum(kb))
  mf <- mp <- matrix(NA_real_, nt, B, dimnames = list(tg, NULL))
  for (i in seq_len(nt)) {
    Bi <- A[(off[i] + 1):off[i + 1], , drop = FALSE]
    mf[i, ] <- matrixStats::colMedians(fold(Bi)); mp[i, ] <- matrixStats::colMedians(Bi)
  }
  j <- match(dtarget[dn], tg)
  st <- ifelse(tdirn[tg] == 1L, 1, -1); sd_ <- ifelse(ddirn[dn] == 1L, 1, -1)
  # directed: a - m in the declared direction, (1 - a) - (1 - m) = -(a - m) against it
  list(tg = tg, dn = dn,
       tf = fold(aT) - mf, td = (aT - mp) * st,
       df = fold(aD) - mf[j, , drop = FALSE], dd = (aD - mp[j, , drop = FALSE]) * sd_)
}
