#' A signature library
#'
#' @param sets named list of gene identifier vectors.
#' @param engine scoring engine per set (recycled): `"mean"` (mean z), `"wmean"` (weighted mean of z,
#'   divided by the sum of absolute weights), `"gsva"`, `"ssgsea"` (both from the GSVA package) or `"pc1"`
#'   (first principal component of the z-scores, signed to correlate positively with the mean z).
#' @param weights optional named list of weight vectors, aligned with `sets`, for `"wmean"` sets.
#' @param direction optional named vector, +1 when a high score is declared to go with the favorable
#'   outcome (response, longer survival), -1 for the unfavorable one, NA when nothing is declared.
#' @return An object of class `dfdr_library`.
#' @examples
#' dfdr_library(list(a = c("G1", "G2", "G3"), b = c("G4", "G5", "G6")), engine = "mean",
#'              direction = c(a = 1, b = -1))
#' @export
dfdr_library <- function(sets, engine = "mean", weights = NULL, direction = NULL) {
  if (is.null(names(sets)) || anyDuplicated(names(sets))) stop("sets must be a uniquely named list")
  nm <- names(sets)
  engine <- stats::setNames(rep_len(engine, length(sets)), nm)
  bad <- setdiff(unique(engine), c("mean", "wmean", "gsva", "ssgsea", "pc1"))
  if (length(bad)) stop("unknown engine: ", paste(bad, collapse = ", "))
  w <- stats::setNames(vector("list", length(sets)), nm)
  if (!is.null(weights)) for (s in intersect(names(weights), nm)) w[[s]] <- as.numeric(weights[[s]])
  for (s in nm[engine == "wmean"]) {
    if (is.null(w[[s]])) w[[s]] <- rep(1, length(sets[[s]]))
    if (length(w[[s]]) != length(sets[[s]])) stop("weights of ", s, " do not match its genes")
  }
  dirn <- stats::setNames(rep(NA_integer_, length(sets)), nm)
  if (!is.null(direction)) dirn[intersect(names(direction), nm)] <- as.integer(direction[intersect(names(direction), nm)])
  structure(list(sets = lapply(sets, as.character), engine = engine, weights = w, direction = dirn),
            class = "dfdr_library")
}

#' @export
print.dfdr_library <- function(x, ...) {
  cat(sprintf("dfdr_library: %d sets; engines %s\n", length(x$sets),
              paste(sprintf("%s %d", names(table(x$engine)), table(x$engine)), collapse = ", ")))
  invisible(x)
}

#' Score gene sets in a cohort
#'
#' @param cohort a `dfdr_cohort`.
#' @param sets named list of gene vectors (only genes measured in the cohort are used; sets with fewer
#'   than `min_genes` measured genes give a row of NA).
#' @param engine one engine for all sets (see [dfdr_library()]).
#' @param weights for `"wmean"`: a named list of weight vectors aligned with `sets` (NULL = unit weights).
#' @param min_genes smallest number of measured genes for a score.
#' @return Score matrix, sets in rows, the cohort's samples in columns.
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy(n_cohorts = 1)
#' score_sets(toy$cohorts[[1]], toy$lib$sets, "mean")[, 1:4]
#' @export
score_sets <- function(cohort, sets, engine, weights = NULL, min_genes = 3L) {
  stopifnot(length(engine) == 1L)
  out <- matrix(NA_real_, length(sets), cohort$n, dimnames = list(names(sets), cohort$samples))
  if (!length(sets)) return(out)
  present <- rownames(cohort$Z)
  gs <- lapply(sets, function(g) g[g %in% present])
  ok <- vapply(gs, length, 1L) >= min_genes
  if (!any(ok)) return(out)
  Z <- cohort$Z
  if (engine == "mean") {
    for (i in which(ok)) out[i, ] <- colMeans(Z[gs[[i]], , drop = FALSE])
  } else if (engine == "wmean") {
    for (i in which(ok)) {
      g <- sets[[i]]; w <- if (is.null(weights[[i]])) rep(1, length(g)) else weights[[i]]
      keep <- g %in% present; g <- g[keep]; w <- w[keep]
      if (sum(abs(w)) == 0) next
      out[i, ] <- as.numeric(crossprod(Z[g, , drop = FALSE], w)) / sum(abs(w))
    }
  } else if (engine == "pc1") {
    for (i in which(ok)) out[i, ] <- pc1_scores(Z[gs[[i]], , drop = FALSE])
  } else {
    ok_ns <- tryCatch({ loadNamespace("GSVA"); TRUE }, error = function(e) conditionMessage(e))
    if (!isTRUE(ok_ns)) stop("engine ", engine, " needs the GSVA package: ", ok_ns)
    par <- if (engine == "ssgsea") GSVA::ssgseaParam(cohort$expr, gs[ok], minSize = min_genes) else
      GSVA::gsvaParam(cohort$expr, gs[ok], minSize = min_genes)
    rr <- GSVA::gsva(par, verbose = FALSE)
    out[rownames(rr), ] <- rr
  }
  out
}

#' Score every set of a library in a cohort, each by its own engine
#' @param cohort a `dfdr_cohort`.
#' @param lib a `dfdr_library`.
#' @return Score matrix, the library's sets in rows, the cohort's samples in columns.
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy(n_cohorts = 1)
#' score_library(toy$cohorts[[1]], toy$lib)[, 1:4]
#' @export
score_library <- function(cohort, lib) {
  out <- matrix(NA_real_, length(lib$sets), cohort$n, dimnames = list(names(lib$sets), cohort$samples))
  for (eng in unique(lib$engine)) {
    nm <- names(lib$sets)[lib$engine == eng]
    out[nm, ] <- score_sets(cohort, lib$sets[nm], eng, lib$weights[nm])
  }
  out
}

# First principal component of the z-scores X (genes x samples; genes already centered), as sample scores, signed
# to correlate positively with the mean z. The leading eigenvector of the smaller Gram matrix gives it: X X' when
# there are fewer genes than samples (scores X' v), X' X otherwise; RSpectra's Lanczos solver is used for large
# matrices when installed, with no random start drawn from R's stream. Random gene sets often have a small gap
# between their first two components, which makes power iteration slow; the Gram route does not depend on it.
pc1_scores <- function(X) {
  k <- nrow(X); n <- ncol(X)
  lead <- function(G) {
    if (nrow(G) > 200L && requireNamespace("RSpectra", quietly = TRUE)) {
      v <- tryCatch(RSpectra::eigs_sym(G, 1L, which = "LA")$vectors[, 1], error = function(e) NULL)
      if (!is.null(v) && all(is.finite(v))) return(v)
    }
    eigen(G, symmetric = TRUE)$vectors[, 1]
  }
  sc <- if (k <= n) as.numeric(crossprod(X, lead(tcrossprod(X)))) else lead(crossprod(X))
  if (stats::cor(sc, colMeans(X)) < 0) sc <- -sc
  sc
}
