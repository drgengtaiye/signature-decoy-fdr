#' Grow decoy modules matched to each signature
#'
#' For every target in the library, grows `n_coex` co-expressed decoys and draws `n_rand` random-gene
#' decoys from a pool of genes measured in every cohort and in no signature. A co-expressed decoy has the
#' target's size and expression-decile composition (the median decile of each gene across cohorts) and is
#' grown one gene at a time: each new gene is drawn from the decile still needed, among the candidates
#' whose mean correlation with the genes already chosen lies within `tol` of the best match to the
#' target's own mean pairwise correlation, the match being the root mean square over cohorts. This is the
#' seed-and-grow of the accompanying paper (calib_decoys_v3.R).
#'
#' @param cohorts list of `dfdr_cohort` objects.
#' @param lib a `dfdr_library`.
#' @param exclude extra genes to keep out of the pool (e.g. genes of related signatures not in `lib`).
#' @param n_coex,n_rand decoys per target and set.
#' @param n_sets number of independent decoy sets. Set 1 is drawn from the target's own L'Ecuyer-CMRG
#'   stream, set `s > 1` from substream `s - 1` of it, so set 1 does not depend on `n_sets`; averaging the
#'   estimate over sets reduces its decoy-draw variability (Keich et al. 2019).
#' @param tol width of the near-best band.
#' @param seed seed; each target draws from its own L'Ecuyer-CMRG stream, so results do not depend on
#'   `cores`.
#' @param cores number of forked workers (parallel::mclapply).
#' @param precompute keep each cohort's pool correlation matrix in memory (fast; needs about
#'   8 x pool^2 bytes per cohort). `NULL` decides from `max_bytes`.
#' @param max_bytes memory allowed for the precomputed matrices.
#' @param coherence `"signed"` matches the mean pairwise correlation of the genes (the coherence a mean or
#'   weighted-mean score depends on); `"absolute"` matches the mean absolute pairwise correlation, for scores that
#'   do not depend on the sign of the correlations, such as the first principal component.
#' @return An object of class `dfdr_decoys`: a list with `decoys` (one entry per decoy: name, target,
#'   class `"coex"` or `"rand"`, set, genes, k, engine, dirn, the target's and the achieved internal
#'   correlation per cohort), `info` (one row per decoy) and `pool`. Decoys are ordered set by set.
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy()
#' grow_decoys(toy$cohorts, toy$lib, n_coex = 4, n_rand = 1, seed = 1)
#' @export
grow_decoys <- function(cohorts, lib, exclude = character(0), n_coex = 32L, n_rand = 8L, n_sets = 1L, tol = 0.02,
                        seed = 1L, cores = 1L, precompute = NULL, max_bytes = 48e9,
                        coherence = c("signed", "absolute")) {
  n_sets <- as.integer(n_sets); stopifnot(n_sets >= 1L)
  coherence <- match.arg(coherence)
  tf <- if (coherence == "absolute") abs else identity
  stopifnot(inherits(lib, "dfdr_library"), length(cohorts) >= 1L)
  CN <- vapply(cohorts, `[[`, "", "name")
  univ <- Reduce(intersect, lapply(cohorts, function(c) rownames(c$Z)))
  pool <- setdiff(univ, c(unlist(lib$sets), exclude))
  if (length(pool) < 100L) stop("the decoy pool has fewer than 100 genes")
  DM <- vapply(cohorts, function(c) as.numeric(c$decile[univ]), numeric(length(univ)))
  refdec <- stats::setNames(as.integer(round(matrixStats::rowMedians(matrix(DM, nrow = length(univ))))), univ)
  Zp <- lapply(cohorts, function(c) c$Z[pool, , drop = FALSE])
  if (is.null(precompute)) precompute <- 8 * length(pool)^2 * length(cohorts) <= max_bytes
  Cl <- if (precompute) lapply(Zp, function(Z) tf(tcrossprod(Z) / (ncol(Z) - 1))) else NULL
  corr_col <- function(j, pick) if (precompute) Cl[[j]][, pick] else
    tf(as.numeric(Zp[[j]] %*% Zp[[j]][pick, ]) / (ncol(Zp[[j]]) - 1))
  mean_r <- function(Z) { if (nrow(Z) < 2) return(NA_real_); C <- tf(tcrossprod(Z) / (ncol(Z) - 1))
    (sum(C) - sum(diag(C))) / (nrow(Z) * (nrow(Z) - 1)) }
  ibA <- split(seq_along(pool), refdec[pool])

  targets <- names(lib$sets)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  streams <- vector("list", length(targets)); s <- .Random.seed
  for (i in seq_along(targets)) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }

  one <- function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    tg <- targets[i]; g <- lib$sets[[tg]]
    rt <- vapply(cohorts, function(c) mean_r(c$Z[intersect(g, rownames(c$Z)), , drop = FALSE]), 1.0)
    kk <- vapply(cohorts, function(c) sum(g %in% rownames(c$Z)), 1L)
    k <- as.integer(round(stats::median(kk)))
    if (k < 3L) return(list())                     # list(), not NULL: NULL marks a dead worker
    dprof <- refdec[intersect(g, univ)]
    if (length(dprof) < k) dprof <- c(dprof, dprof[sample.int(length(dprof), k - length(dprof), replace = TRUE)])
    if (length(dprof) > k) dprof <- dprof[sample.int(length(dprof), k)]
    dn <- table(factor(dprof, levels = 1:10))
    use_c <- which(is.finite(rt))
    grow <- function() {
      need <- dn[dn > 0]; S <- integer(0)
      cs <- matrix(0, length(use_c), length(pool))
      while (length(S) < k) {
        av <- names(need)[need > 0]
        if (!length(av)) { av <- as.character(1:10); need[av[!av %in% names(need)]] <- 0L }
        d <- if (length(av) == 1L) av else sample(av, 1L, prob = pmax(need[av], 1e-9))
        cand <- setdiff(ibA[[d]], S)
        if (!length(cand)) { need[d] <- 0L; next }
        pick <- if (!length(S) || !length(use_c)) cand[sample.int(length(cand), 1L)] else {
          L <- sqrt(colMeans((cs[, cand, drop = FALSE] / length(S) - rt[use_c])^2))
          okk <- which(L <= min(L) + tol); cand[okk[sample.int(length(okk), 1L)]] }
        S <- c(S, pick)
        for (j in seq_along(use_c)) cs[j, ] <- cs[j, ] + corr_col(use_c[j], pick)
        need[d] <- max(0L, need[d] - 1L)
      }
      S
    }
    mk <- function(idx, cls, b, set) {
      gg <- pool[idx]
      nm <- if (set == 1L) sprintf("dcy_%s_%s_%02d", cls, tg, b) else sprintf("dcy_%s_%s_s%02d_%02d", cls, tg, set, b)
      list(name = nm, target = tg, class = cls, set = set, genes = gg, k = length(gg),
           engine = unname(lib$engine[tg]), dirn = sample(c(1L, -1L), 1L), r_target = rt,
           r_achieved = vapply(seq_along(cohorts), function(j) mean_r(Zp[[j]][idx, , drop = FALSE]), 1.0))
    }
    out <- list(); sub <- streams[[i]]
    for (set in seq_len(n_sets)) {
      if (set > 1L) { sub <- parallel::nextRNGSubStream(sub); assign(".Random.seed", sub, envir = globalenv()) }
      for (b in seq_len(n_coex)) out[[length(out) + 1L]] <- mk(grow(), "coex", b, set)
      for (b in seq_len(n_rand)) out[[length(out) + 1L]] <- mk(sample.int(length(pool), k), "rand", b, set)
    }
    out
  }
  if (cores > 1L) invisible(gc())   # a forked worker copies whatever garbage its own collections touch
  res <- if (cores > 1L) parallel::mclapply(seq_along(targets), one, mc.cores = cores, mc.preschedule = FALSE) else
    lapply(seq_along(targets), one)
  err <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)   # NULL: a forked worker died (often memory)
  if (any(err)) stop("decoy growth failed for ", paste(targets[err], collapse = ", "), ": ",
                     if (is.null(res[[which(err)[1]]])) "a worker returned no result (out of memory? use fewer cores)" else
                       as.character(res[[which(err)[1]]]))
  DEC <- unlist(res, recursive = FALSE)
  if (!length(DEC)) stop("no target has 3 or more genes in the cohorts; do the gene identifiers match?")
  DEC <- DEC[order(vapply(DEC, `[[`, 1L, "set"))]              # set by set; stable within a set
  names(DEC) <- vapply(DEC, `[[`, "", "name")
  info <- do.call(rbind, lapply(DEC, function(d) data.frame(decoy = d$name, target = d$target, class = d$class,
    set = d$set, k = d$k, engine = d$engine, dirn = d$dirn, r_target_med = stats::median(d$r_target, na.rm = TRUE),
    r_achieved_med = stats::median(d$r_achieved, na.rm = TRUE),
    genes = paste(d$genes, collapse = ";"), stringsAsFactors = FALSE)))
  rownames(info) <- NULL
  structure(list(decoys = DEC, info = info, pool = pool, cohorts = CN, coherence = coherence), class = "dfdr_decoys")
}

#' @export
print.dfdr_decoys <- function(x, ...) {
  i <- x$info
  ns <- if (is.null(i$set)) 1L else length(unique(i$set))
  cat(sprintf("dfdr_decoys: %d decoys for %d targets in %d set%s (pool %d genes, %d cohorts)\n", nrow(i),
              length(unique(i$target)), ns, if (ns == 1L) "" else "s", length(x$pool), length(x$cohorts)))
  for (cl in unique(i$class)) {
    j <- i$class == cl
    cat(sprintf("  %-4s internal r, median over decoys: target %.3f, achieved %.3f\n", cl,
                stats::median(i$r_target_med[j], na.rm = TRUE), stats::median(i$r_achieved_med[j], na.rm = TRUE)))
  }
  invisible(x)
}
