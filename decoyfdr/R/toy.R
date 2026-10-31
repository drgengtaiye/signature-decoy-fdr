#' A small simulated library for examples
#'
#' Draws from the current random number stream; set a seed first for a reproducible library.
#'
#' Cohorts that share two co-expressed programs, the first of which drives a binary outcome, and a library
#' of four signatures: two covering the response program, one sampling the other program and one of noise
#' genes (the setting of the vignette, smaller).
#'
#' @param n_cohorts number of cohorts.
#' @param n samples per cohort.
#' @param n_genes genes (at least 300).
#' @return A list with `cohorts` (a list of `dfdr_cohort` objects) and `lib` (a `dfdr_library`).
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy()
#' toy$cohorts[[1]]
#' toy$lib
#' @export
dfdr_toy <- function(n_cohorts = 6L, n = 40L, n_genes = 300L) {
    if (n_genes < 300L) stop("n_genes must be at least 300")
    genes <- sprintf("G%03d", seq_len(n_genes))
    cohorts <- lapply(seq_len(n_cohorts), function(i) {
        f <- stats::rnorm(n)
        X <- matrix(stats::rnorm(n_genes * n), n_genes, n,
                    dimnames = list(genes, sprintf("c%d_%02d", i, seq_len(n))))
        X[1:40, ] <- X[1:40, ] + matrix(f, 40, n, byrow = TRUE)
        X[101:250, ] <- X[101:250, ] + matrix(stats::rnorm(n), 150, n, byrow = TRUE)
        y <- stats::rbinom(n, 1, stats::plogis(1.2 * f))
        names(y) <- colnames(X)
        dfdr_cohort(X + 5, y, name = paste0("C", i))
    })
    lib <- dfdr_library(list(prog_a = genes[1:20], prog_b = genes[21:40], other = genes[101:115],
                             noise = genes[270:285]),
                        engine = "mean", direction = c(prog_a = 1, prog_b = 1, other = 1, noise = 1))
    list(cohorts = cohorts, lib = lib)
}
