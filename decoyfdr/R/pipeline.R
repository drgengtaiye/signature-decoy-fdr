#' Score targets, decoys and matched null banks, cohort by cohort
#'
#' For every cohort: each target is scored by its engine; `n_draw` random gene sets of the target's
#' measured size are drawn from all measured genes and scored the same way (a weighted-mean target's own
#' weights, permuted, go to each draw, so weights of both signs are kept); each decoy is scored by its
#' target's engine, with its target's weights permuted; and every unit's signed statistic is read
#' against its target's bank under the folded and the directed transformation.
#'
#' @param cohorts list of `dfdr_cohort` objects (all binary or all survival).
#' @param lib a `dfdr_library`.
#' @param decoys a `dfdr_decoys` object from [grow_decoys()].
#' @param n_draw null draws per target and cohort.
#' @param keep_draws number of draws whose scores are kept for [calibrate_null()] (0 keeps none).
#' @param seed seed; each cohort draws from its own L'Ecuyer-CMRG stream.
#' @param cores forked workers over cohorts.
#' @return An object of class `dfdr_fit`: `units` (one row per unit and cohort: the signed statistic and
#'   the folded and directed observed value, bank median and 95th percentile) and, if `keep_draws > 0`,
#'   the scores needed for calibration.
#' @export
fit_units <- function(cohorts, lib, decoys, n_draw = 1000L, keep_draws = 0L, seed = 1L, cores = 1L) {
  types <- unique(vapply(cohorts, `[[`, "", "type"))
  if (length(types) != 1L) stop("cohorts must share one outcome type")
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  streams <- vector("list", length(cohorts)); s <- .Random.seed
  for (i in seq_along(cohorts)) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }
  DEC <- decoys$decoys
  # load GSVA once here: workers loading it at the same time race to create HDF5Array's dump directory
  if (any(lib$engine %in% c("gsva", "ssgsea"))) loadNamespace("GSVA")
  one <- function(ci) {
    assign(".Random.seed", streams[[ci]], envir = globalenv())
    co <- cohorts[[ci]]; present <- rownames(co$Z)
    tg <- names(lib$sets)
    tset <- lapply(lib$sets, function(g) g[g %in% present])
    tw <- lapply(tg, function(s) if (lib$engine[s] == "wmean") lib$weights[[s]][lib$sets[[s]] %in% present] else NULL)
    names(tw) <- tg
    tg <- tg[vapply(tset[tg], length, 1L) >= 3L]
    obs_t <- score_library(co, dfdr_library(tset[tg], lib$engine[tg], tw[tg], lib$direction[tg]))
    a_t <- outcome_stat(co, obs_t)
    # null banks: all draws are made first (target by target, in library order), then every engine's
    # draws are scored in one call, which keeps the GSVA kernel step to one per engine and cohort
    bank <- list(); keepS <- list(); bsets <- list(); bw <- list()
    for (s in tg) {
      k <- length(tset[[s]]); eng <- lib$engine[[s]]
      sets <- replicate(n_draw, present[sample.int(length(present), k)], simplify = FALSE)
      names(sets) <- sprintf("%s__n%04d", s, seq_len(n_draw))
      bsets[[s]] <- sets
      if (eng == "wmean") { w <- lapply(seq_len(n_draw), function(i) sample(tw[[s]])); names(w) <- names(sets); bw[[s]] <- w }
    }
    for (eng in unique(lib$engine[tg])) {
      ss <- tg[lib$engine[tg] == eng]
      S <- score_sets(co, unlist(unname(bsets[ss]), recursive = FALSE), eng,
                      if (eng == "wmean") unlist(unname(bw[ss]), recursive = FALSE) else NULL)
      A <- outcome_stat(co, S)
      for (s in ss) {
        i <- match(names(bsets[[s]]), rownames(S))
        bank[[s]] <- A[i]
        if (keep_draws > 0L) keepS[[s]] <- S[i[seq_len(min(keep_draws, n_draw))], , drop = FALSE]
      }
      rm(S)
    }
    # decoys
    dn <- names(DEC)[vapply(DEC, function(d) d$target %in% tg, TRUE)]
    obs_d <- matrix(NA_real_, length(dn), co$n, dimnames = list(dn, co$samples))
    for (eng in unique(lib$engine[tg])) {
      nm <- dn[vapply(DEC[dn], function(d) lib$engine[[d$target]] == eng, TRUE)]
      if (!length(nm)) next
      sets <- lapply(DEC[nm], function(d) d$genes[d$genes %in% present])
      w <- if (eng == "wmean") lapply(nm, function(n_) { ws <- tw[[DEC[[n_]]$target]]; L <- length(sets[[n_]])
        sample(ws, L, replace = length(ws) != L) }) else NULL
      if (!is.null(w)) names(w) <- nm
      obs_d[nm, ] <- score_sets(co, sets, eng, w)
    }
    a_d <- outcome_stat(co, obs_d)
    rd <- function(a, b, dirn) {
      dirn <- ifelse(is.na(dirn), 1L, dirn)
      bf <- fold(b); bd <- if (dirn == 1L) b else 1 - b
      c(stat = a, fold_obs = fold(a), fold_med = stats::median(bf), fold_p95 = unname(stats::quantile(bf, 0.95)),
        dir_obs = if (dirn == 1L) a else 1 - a, dir_med = stats::median(bd), dir_p95 = unname(stats::quantile(bd, 0.95)))
    }
    rows_t <- lapply(tg, function(s) c(rd(a_t[[match(s, tg)]], bank[[s]], lib$direction[[s]])))
    rows_d <- lapply(seq_along(dn), function(i) rd(a_d[[i]], bank[[DEC[[dn[i]]]$target]], DEC[[dn[i]]]$dirn))
    U <- data.frame(cohort = co$name, n = co$n,
      unit = c(tg, dn), kind = c(rep("target", length(tg)), rep("decoy", length(dn))),
      target = c(tg, vapply(DEC[dn], `[[`, "", "target")),
      class = c(rep("signature", length(tg)), vapply(DEC[dn], `[[`, "", "class")),
      set = c(rep(NA_integer_, length(tg)), vapply(DEC[dn], function(d) if (is.null(d$set)) 1L else as.integer(d$set), 1L)),
      k = c(vapply(tset[tg], length, 1L), vapply(DEC[dn], function(d) sum(d$genes %in% present), 1L)),
      dirn = c(unname(lib$direction[tg]), vapply(DEC[dn], `[[`, 1L, "dirn")),
      do.call(rbind, c(rows_t, rows_d)), stringsAsFactors = FALSE)
    rownames(U) <- NULL
    keep <- if (keep_draws > 0L) list(targets = obs_t, decoys = obs_d, bank = keepS) else NULL
    list(units = U, scores = keep)
  }
  if (cores > 1L) invisible(gc())   # a forked worker copies whatever garbage its own collections touch
  res <- if (cores > 1L) parallel::mclapply(seq_along(cohorts), one, mc.cores = cores, mc.preschedule = FALSE) else
    lapply(seq_along(cohorts), one)
  err <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)   # NULL: a forked worker died (often memory)
  if (any(err)) stop("scoring failed in ", paste(vapply(cohorts[err], `[[`, "", "name"), collapse = ", "), ": ",
                     if (is.null(res[[which(err)[1]]])) "a worker returned no result (out of memory? use fewer cores)" else
                       as.character(res[[which(err)[1]]]))
  structure(list(units = do.call(rbind, lapply(res, `[[`, "units")),
                 scores = stats::setNames(lapply(res, `[[`, "scores"), vapply(cohorts, `[[`, "", "name")),
                 cohorts = cohorts, lib = lib, decoys = decoys, type = types),
            class = "dfdr_fit")
}

#' Call every unit: Wilcoxon signed-rank over cohorts
#'
#' For each unit (target or decoy) and statistic, delta = observed - bank median in each cohort; the
#' unit is called when the two-sided signed-rank P over cohorts is below `alpha` and the median delta is
#' positive. Units tested in fewer than `min_cohorts` cohorts get P = 1.
#'
#' @param units a units table (the `units` element of a `dfdr_fit`, or a table with the same columns:
#'   `unit, kind, target, class, dirn` and `<stat>_obs, <stat>_med`).
#' @param stats statistics to call, any of `"fold"` and `"dir"`.
#' @param alpha,min_cohorts the call rule.
#' @param exclude optional data frame with columns `target` and `cohort`: the target's tests in those
#'   cohorts, and its decoys', are dropped (e.g. the cohorts a signature was developed in).
#' @return One row per unit and statistic: `unit, kind, target, class, set, dirn, stat, ncoh, med_delta, p,
#'   call, neg` (`set` is the decoy set, NA for targets; absent if `units` has no `set` column).
#' @examples
#' U <- read.csv(system.file("extdata", "response_units_N0.csv.gz", package = "decoyfdr"))
#' calls <- call_units(U)
#' table(calls$stat[calls$kind == "target"], calls$call[calls$kind == "target"])
#' @export
call_units <- function(units, stats = c("fold", "dir"), alpha = 0.05, min_cohorts = 6L, exclude = NULL) {
  U <- units
  if (!"class" %in% names(U) && all(c("pool", "structure") %in% names(U)))
    U$class <- ifelse(U$kind == "target", "signature", paste0(U$structure, U$pool))
  if (!is.null(exclude) && nrow(exclude)) {
    drop <- paste(U$target, U$cohort) %in% paste(exclude$target, exclude$cohort)
    U <- U[!drop, , drop = FALSE]
  }
  out <- list()
  for (st in stats) {
    d <- U[[paste0(st, "_obs")]] - U[[paste0(st, "_med")]]
    ok <- is.finite(d)
    sp <- split(d[ok], U$unit[ok])
    first <- U[!duplicated(U$unit), intersect(c("unit", "kind", "target", "class", "set", "dirn"), names(U))]
    first <- first[match(names(sp), first$unit), ]
    ncoh <- vapply(sp, length, 1L)
    p <- vapply(sp, function(x) if (length(x) >= min_cohorts && any(x != 0)) signed_rank_p1(x) else 1, 1.0)
    m <- vapply(sp, stats::median, 1.0)
    out[[st]] <- data.frame(first, stat = st, ncoh = ncoh, med_delta = m, p = p,
                            call = p < alpha & m > 0, neg = p < alpha & m < 0, row.names = NULL, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
