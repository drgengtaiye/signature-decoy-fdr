#' Prepare one cohort
#'
#' Keeps the samples that have an outcome, the genes that are finite and non-constant across them (the
#' first row of a duplicated gene identifier), and computes gene-wise z-scores and the expression decile
#' of every gene (deciles of the gene means within the cohort).
#'
#' @param expr numeric matrix, genes in rows (row names are gene identifiers), samples in columns
#'   (column names are sample identifiers).
#' @param outcome either a 0/1 vector (1 = responder or favorable class) or a data frame / matrix with
#'   columns `time` and `event` (1 = event). Named by sample, or in the column order of `expr`.
#' @param name cohort name.
#' @return An object of class `dfdr_cohort`.
#' @export
dfdr_cohort <- function(expr, outcome, name = "cohort") {
  expr <- as.matrix(expr)
  if (is.null(rownames(expr)) || is.null(colnames(expr))) stop("expr needs row and column names")
  if (is.data.frame(outcome) || is.matrix(outcome)) {
    outcome <- as.data.frame(outcome)
    if (!all(c("time", "event") %in% names(outcome))) stop("a survival outcome needs columns time and event")
    sid <- if (!is.null(rownames(outcome)) && all(colnames(expr) %in% rownames(outcome)) &&
               !identical(rownames(outcome), as.character(seq_len(nrow(outcome))))) rownames(outcome) else colnames(expr)
    if (nrow(outcome) != length(sid)) stop("outcome rows must match the samples")
    tm <- stats::setNames(as.numeric(outcome$time), sid)
    ev <- stats::setNames(as.integer(outcome$event), sid)
    keep <- intersect(colnames(expr), sid[is.finite(tm) & tm > 0 & !is.na(ev)])
    type <- "survival"
    out <- list(time = unname(tm[keep]), event = unname(ev[keep]))
  } else {
    y <- outcome
    if (is.null(names(y))) { if (length(y) != ncol(expr)) stop("unnamed outcome must match ncol(expr)"); names(y) <- colnames(expr) }
    keep <- intersect(colnames(expr), names(y)[!is.na(y)])
    type <- "binary"
    out <- list(y = as.integer(y[keep]))
    if (!all(out$y %in% 0:1)) stop("a binary outcome must be 0/1")
  }
  M <- expr[, keep, drop = FALSE]
  M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & matrixStats::rowSds(M) > 0, , drop = FALSE]
  s <- matrixStats::rowSds(M)
  Z <- (M - rowMeans(M)) / s
  mu <- rowMeans(M)
  brk <- stats::quantile(mu, probs = seq(0, 1, by = 0.1), na.rm = TRUE); brk[1] <- -Inf; brk[11] <- Inf
  dec <- stats::setNames(cut(mu, breaks = brk, labels = FALSE, include.lowest = TRUE), rownames(M))
  structure(c(list(name = name, type = type, expr = M, Z = Z, decile = dec, n = ncol(M), samples = keep), out),
            class = "dfdr_cohort")
}

#' @export
print.dfdr_cohort <- function(x, ...) {
  ev <- if (x$type == "binary") sprintf("%d responders", sum(x$y)) else sprintf("%d events", sum(x$event))
  cat(sprintf("dfdr_cohort '%s': %d genes x %d samples, %s outcome, %s\n", x$name, nrow(x$Z), x$n, x$type, ev))
  invisible(x)
}

#' Shuffle the outcome of a cohort
#'
#' Returns the cohort with its outcome permuted across samples (for survival, (time, event) pairs are
#' permuted together), which makes every score null.
#' @param cohort a `dfdr_cohort`.
#' @param perm a permutation of `seq_len(cohort$n)`; random if omitted.
#' @export
permute_outcome <- function(cohort, perm = sample.int(cohort$n)) {
  if (cohort$type == "binary") cohort$y <- cohort$y[perm] else {
    cohort$time <- cohort$time[perm]; cohort$event <- cohort$event[perm] }
  cohort
}

#' The signed statistic of every score in a cohort
#'
#' AUC for a binary outcome, Harrell's C (larger score, longer survival) for survival.
#' @param cohort a `dfdr_cohort`.
#' @param S score matrix, rows = scores, columns = the cohort's samples.
#' @export
outcome_stat <- function(cohort, S) {
  if (cohort$type == "binary") auc_rows(S, cohort$y) else cindex_rows(S, cohort$time, cohort$event)
}
