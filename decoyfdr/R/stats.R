#' Signed outcome statistics for many scores at once
#'
#' `auc_rows()` gives, for each row of a score matrix, the probability that a responder scores above a
#' non-responder (ties one half). `cindex_rows()` gives Harrell's C with a larger score going with longer
#' survival, with the conventions of `survival::concordance()`. Both equal 0.5 for an uninformative score;
#' a statistic below 0.5 means the score runs the other way. `fold()` maps a statistic to `max(s, 1 - s)`.
#'
#' @param S numeric matrix, one row per score, one column per sample.
#' @param y 0/1 vector of length `ncol(S)` (1 = responder, or the favorable class).
#' @param time,event survival time and event indicator (1 = event), of length `ncol(S)`.
#' @return A numeric vector with one value per row of `S`.
#' @examples
#' S <- rbind(a = c(1, 2, 3, 4, 5, 6), b = c(6, 1, 5, 2, 4, 3))
#' auc_rows(S, c(0, 0, 0, 1, 1, 1))
#' fold(auc_rows(S, c(0, 0, 0, 1, 1, 1)))
#' cindex_rows(S, time = c(5, 4, 6, 2, 1, 3), event = c(1, 1, 0, 1, 1, 1))
#' @export
auc_rows <- function(S, y) {
  S <- as_row_matrix(S)
  y <- as.integer(y)
  if (length(y) != ncol(S) || anyNA(y) || !all(y %in% 0:1)) stop("y must be 0/1, one entry per column")
  n1 <- sum(y); n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) stop("y needs both classes")
  R <- matrixStats::rowRanks(S, ties.method = "average")
  (as.numeric(R %*% y) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @rdname auc_rows
#' @export
cindex_rows <- function(S, time, event) {
  S <- as_row_matrix(S)
  event <- as.integer(event)
  if (length(time) != ncol(S) || length(event) != ncol(S)) stop("time and event need one entry per column")
  cindex_rows_cpp(S, as.numeric(time), event)
}

#' @rdname auc_rows
#' @param s a signed statistic.
#' @export
fold <- function(s) pmax(s, 1 - s)

as_row_matrix <- function(S) {
  if (is.null(dim(S))) S <- matrix(S, nrow = 1L)
  storage.mode(S) <- "double"
  S
}

#' Two-sided Wilcoxon signed-rank P values, column by column
#'
#' Exact P values as SciPy computes them: the exact null distribution of the signed-rank statistic when
#' no difference is zero and none is tied; with zeros (dropped) or ties (mid-ranks), the exact
#' distribution over all sign patterns of the non-zero ranks, for up to `max_exact` non-zero
#' differences; beyond that, the normal approximation with the tie correction and no continuity
#' correction. A column whose differences are all zero gets P = 1.
#'
#' @param D numeric matrix, units in columns, cohorts in rows (or a vector for one unit).
#' @param max_exact largest number of non-zero differences for the exact distributions.
#' @return A vector of P values, one per column.
#' @examples
#' D <- cbind(up = c(0.04, 0.02, 0.05, 0.01, 0.03, 0.06), none = c(0.02, -0.03, 0.01, -0.01, 0.04, -0.02))
#' signed_rank_p(D)
#' @export
signed_rank_p <- function(D, max_exact = 20L) {
  if (is.null(dim(D))) D <- matrix(D, ncol = 1L)
  apply(D, 2L, signed_rank_p1, max_exact = max_exact)
}

signed_rank_p1 <- function(d, max_exact = 20L) {
  d <- d[is.finite(d)]
  nz <- d[d != 0]
  k <- length(nz)
  if (k == 0L) return(1)
  r <- rank(abs(nz))
  w <- sum(r[nz > 0])
  ties <- anyDuplicated(abs(nz)) > 0
  if (k == length(d) && !ties && k <= 50L) {
    cdf <- signed_rank_cdf(k)
    tt <- min(w, k * (k + 1) / 2 - w)
    return(min(1, 2 * cdf[round(tt) + 1L]))
  }
  if (k <= max_exact) {
    tt <- min(w, sum(r) - w)
    W <- sign_patterns(k) %*% r
    return(min(1, 2 * mean(W <= tt + 1e-9)))
  }
  m <- k * (k + 1) / 4
  tie_tab <- table(r)
  v <- k * (k + 1) * (2 * k + 1) / 24 - sum(tie_tab^3 - tie_tab) / 48
  z <- (w - m) / sqrt(v)
  min(1, 2 * stats::pnorm(-abs(z)))
}

signed_rank_cdf <- local({
  cache <- list()
  function(n) {
    key <- as.character(n)
    if (!is.null(cache[[key]])) return(cache[[key]])
    cnt <- numeric(n * (n + 1) / 2 + 1); cnt[1] <- 1
    for (r in seq_len(n)) {
      nxt <- cnt
      idx <- (r + 1):length(cnt)
      nxt[idx] <- nxt[idx] + cnt[idx - r]
      cnt <- nxt
    }
    out <- cumsum(cnt) / 2^n
    cache[[key]] <<- out
    out
  }
})

sign_patterns <- local({
  cache <- list()
  function(k) {
    key <- as.character(k)
    if (!is.null(cache[[key]])) return(cache[[key]])
    M <- outer(0:(2^k - 1), 0:(k - 1), function(a, b) bitwAnd(bitwShiftR(a, b), 1L))
    storage.mode(M) <- "double"
    cache[[key]] <<- M
    M
  }
})
