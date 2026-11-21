#' Check prognostic risk models against matched decoys
#'
#' The audit of published prognostic risk models of the accompanying paper, for one model or many. A model is a vector
#' of coefficients named by gene; its score in a cohort is the sum of coefficient times the gene's within-cohort
#' z-score, and a high score is taken to mark shorter survival, as for a Cox risk score. Each model of three or more
#' genes receives `n_coex` co-expressed and `n_rand` random-gene decoys ([grow_decoys()]): gene sets of its size and
#' expression-level composition, the co-expressed ones also matched on internal co-expression, drawn from genes in no
#' model. A decoy carries the model's coefficients in a random order and the direction of its own association with
#' survival in the training cohort. Model and decoys then go through the same test in every validation cohort:
#' rule (L), the score split at its median, log-rank P < `alpha` with the high group faring worse.
#'
#' Genes. A model keeps the genes measured in every cohort of `match`; if some are missing there it keeps the rest when
#' they are at least `min_measured` of its printed genes (`n_printed`) and at least two, and is dropped otherwise. In a
#' validation cohort outside `match`, a model or decoy is scored when at least `min_measured` of its genes are
#' measured there, on those genes.
#'
#' Reading the result. `d_coex` is how often a model with this model's structure, built from genes unrelated to it,
#' passes the same test in the same cohort. Summed over a set of models and divided by the number validated, it
#' estimates the share of validations that are false ([validation_fdr()]). `rho` is the share of co-expressed decoys
#' the model beats on the directed C-index: it ranks models, it does not sort true from false. In the paper, among
#' published models validated in one cohort, those with rho of at least 0.95 there replicated in another independent
#' cohort in 30.5% of pairs and the others in 21.9%. Validation rates also depend on the cohort ([cohort_fvr()]).
#'
#' @param models a numeric vector of coefficients named by gene (one model), a named list of them, or a data frame with
#'   columns `model`, `gene` and `coef`.
#' @param train the `dfdr_cohort` (survival outcome) the models were fitted in, such as TCGA of the cancer type; it sets
#'   the decoys' directions.
#' @param validate a `dfdr_cohort` or a list of them (survival outcomes): the cohorts the models are tested in.
#' @param match cohorts the decoys are matched in (`train` is added first if missing). Default: `train` and every
#'   validation cohort.
#' @param n_printed number of genes of each model as printed (named by model); default the number given.
#' @param exclude genes kept out of the decoy pool besides the models' own (for instance the gene sets the models were
#'   selected from).
#' @param n_coex,n_rand decoys per model.
#' @param min_measured smallest share of a model's genes that must be measured.
#' @param alpha log-rank threshold of rule (L).
#' @param seed seed of the decoys and of the coefficient shuffles; with the same models in the same order, the same
#'   cohorts and the same seed, the result is identical.
#' @param cores forked workers.
#' @param max_bytes memory allowed for the precomputed correlation matrices of [grow_decoys()].
#' @return An object of class `dfdr_check`, a list with
#'   \describe{
#'   \item{units}{one row per model and validation cohort: genes used, `scorable`, the directed C-index `cdir` (the
#'     probability that of two patients the one with the higher score dies first), log-rank `p`, `worse` (high group
#'     fares worse), `validated`, the training-cohort statistics, and, for models with decoys, `n_coex`, `d_coex`,
#'     `n_rand`, `d_rand` and `rho`.}
#'   \item{scores}{one row per model or decoy and cohort (training cohort included).}
#'   \item{decoys}{one row per decoy: model, class, genes, shuffled coefficients, direction.}
#'   \item{models}{one row per model: genes given, printed and used, and whether it was kept.}
#'   }
#'   and the names of the cohorts and the seed.
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy_survival()
#' chk <- check_models(toy$models, toy$train, toy$validate, n_coex = 8, n_rand = 2, seed = 1)
#' chk
#' validation_fdr(chk, B = 200)
#' @export
check_models <- function(models, train, validate, match = NULL, n_printed = NULL, exclude = character(0),
                         n_coex = 32L, n_rand = 8L, min_measured = 0.8, alpha = 0.05, seed = 1L, cores = 1L,
                         max_bytes = 48e9) {
  models <- as_models(models)
  if (inherits(validate, "dfdr_cohort")) validate <- list(validate)
  if (is.null(match)) match <- c(list(train), validate)
  nm_of <- function(L) vapply(L, `[[`, "", "name")
  if (!train$name %in% nm_of(match)) match <- c(list(train), match)
  for (co in c(list(train), validate, match)) {
    if (!inherits(co, "dfdr_cohort")) stop("train, validate and match take dfdr_cohort objects")
    if (co$type != "survival") stop("cohort '", co$name, "' has no survival outcome")
  }
  if (anyDuplicated(nm_of(validate))) stop("validation cohorts need distinct names")
  all_co <- c(list(train), validate, match); all_co <- all_co[!duplicated(nm_of(all_co))]
  CN <- nm_of(all_co); MN <- nm_of(match); VN <- nm_of(validate)
  if (is.null(n_printed)) n_printed <- lengths(models)
  n_printed <- stats::setNames(as.integer(n_printed[names(models)]), names(models))
  if (anyNA(n_printed)) stop("n_printed lacks a model")

  old <- rng_save(); on.exit(rng_restore(old), add = TRUE)

  # genes measured in every matching cohort
  univ <- Reduce(intersect, lapply(match, function(c) rownames(c$Z)))
  sets <- lapply(models, names); W <- lapply(models, unname)
  status <- stats::setNames(rep("kept", length(models)), names(models))
  for (m in names(sets)) {
    k <- sets[[m]] %in% univ
    if (all(k)) next
    if (sum(k) / n_printed[[m]] < min_measured || sum(k) < 2) {
      status[[m]] <- sprintf("dropped: %d of %d printed genes measured in every matching cohort", sum(k), n_printed[[m]])
      sets[[m]] <- NULL; W[[m]] <- NULL
    } else { sets[[m]] <- sets[[m]][k]; W[[m]] <- W[[m]][k] }
  }
  MD <- data.frame(model = names(models), genes_given = lengths(models), genes_printed = n_printed,
                   genes_used = vapply(names(models), function(m) if (is.null(sets[[m]])) 0L else length(sets[[m]]), 1L),
                   status = unname(status), stringsAsFactors = FALSE, row.names = NULL)
  if (!length(sets)) stop("no model has enough genes measured in the matching cohorts")

  # decoys (models of three or more genes), their coefficients and directions
  big <- names(sets)[lengths(sets) >= 3]
  dec <- NULL; wd <- list(); ddir <- integer(0)
  if (length(big)) {
    lib <- dfdr_library(sets[big], engine = "wmean", weights = W[big])
    dec <- grow_decoys(match, lib, exclude = exclude, n_coex = n_coex, n_rand = n_rand, n_sets = 1L, seed = seed,
                       cores = cores, max_bytes = max_bytes)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
    wd <- lapply(dec$decoys, function(d) { w <- W[[d$target]]; stopifnot(length(w) == d$k); w[sample.int(length(w))] })
    ddir <- vapply(seq_along(dec$decoys), function(i) {
      st <- lin_score(train, dec$decoys[[i]]$genes, wd[[i]])
      if (cindex_rows(st, train$time, train$event) <= 0.5) 1L else -1L     # C <= 1/2: higher score, shorter survival
    }, 1L)
  }
  U <- c(lapply(names(sets), function(m) list(kind = "model", id = m, model = m, class = NA_character_,
                                              genes = sets[[m]], w = W[[m]], given = names(models[[m]]),
                                              wg = unname(models[[m]]), dirn = 1L)),
         if (!is.null(dec)) lapply(seq_along(dec$decoys), function(i) {
           d <- dec$decoys[[i]]
           list(kind = "decoy", id = d$name, model = d$target, class = d$class, genes = d$genes, w = wd[[i]],
                given = d$genes, wg = wd[[i]], dirn = ddir[i]) }))

  # every unit in every cohort
  one_cohort <- function(co) {
    in_match <- co$name %in% MN
    r <- parallel::mclapply(U, function(u) {
      g <- if (in_match) u$genes else u$given
      w <- if (in_match) u$w else u$wg
      k <- g %in% rownames(co$Z)
      base <- data.frame(cohort = co$name, kind = u$kind, id = u$id, model = u$model, class = u$class,
                         genes = length(g), measured = sum(k), stringsAsFactors = FALSE)
      if (mean(k) < min_measured || sum(k) < 2)
        return(cbind(base, scorable = FALSE, cdir = NA_real_, p = NA_real_, worse = NA_real_, validated = NA_real_))
      s <- lin_score(co, g[k], w[k])
      ce <- cindex_rows(s, co$time, co$event); kv <- km_valid(u$dirn * s, co$time, co$event)
      cbind(base, scorable = TRUE, cdir = if (u$dirn > 0) 1 - ce else ce, p = kv[["p"]], worse = kv[["worse"]],
            validated = as.numeric(isTRUE(kv[["p"]] < alpha) && isTRUE(kv[["worse"]] == 1)))
    }, mc.cores = cores, mc.preschedule = TRUE)
    bad <- vapply(r, function(x) is.null(x) || inherits(x, "try-error"), TRUE)
    if (any(bad)) stop("scoring failed in cohort ", co$name, if (is.null(r[[which(bad)[1]]])) " (a worker died)" else
      paste0(": ", as.character(r[[which(bad)[1]]])))
    do.call(rbind, r)
  }
  S <- do.call(rbind, lapply(all_co[CN %in% c(train$name, VN)], one_cohort))
  rownames(S) <- NULL

  # per model and validation cohort
  tr <- S[S$cohort == train$name & S$kind == "model", c("id", "cdir", "p", "worse", "validated")]
  names(tr) <- c("model", "train_cdir", "train_p", "train_worse", "train_validated")
  units <- lapply(VN, function(v) {
    x <- S[S$cohort == v, ]
    t <- x[x$kind == "model", c("cohort", "id", "genes", "measured", "scorable", "cdir", "p", "worse", "validated")]
    names(t)[2] <- "model"
    d <- x[x$kind == "decoy" & x$scorable, ]
    agg <- lapply(t$model, function(m) {
      dc <- d[d$model == m & d$class == "coex", ]; dr <- d[d$model == m & d$class == "rand", ]
      ct <- t$cdir[t$model == m]
      data.frame(n_coex = nrow(dc), d_coex = if (nrow(dc)) mean(dc$validated) else NA_real_,
                 n_rand = nrow(dr), d_rand = if (nrow(dr)) mean(dr$validated) else NA_real_,
                 rho = if (nrow(dc) && is.finite(ct)) mean((dc$cdir < ct) + 0.5 * (dc$cdir == ct)) else NA_real_)
    })
    cbind(t, do.call(rbind, agg))
  })
  units <- merge(do.call(rbind, units), tr, by = "model", sort = FALSE)
  units <- units[order(match(units$cohort, VN), match(units$model, names(models))), ]
  rownames(units) <- NULL
  DI <- if (is.null(dec)) NULL else data.frame(decoy = dec$info$decoy, model = dec$info$target, class = dec$info$class,
    k = dec$info$k, genes = dec$info$genes,
    coefs = vapply(wd, function(w) paste(signif(w, 8), collapse = ";"), ""), dirn = ddir,
    r_target_med = dec$info$r_target_med, r_achieved_med = dec$info$r_achieved_med, stringsAsFactors = FALSE)
  structure(list(units = units, scores = S, decoys = DI, models = MD, train = train$name, match = MN, validate = VN,
                 seed = seed, alpha = alpha, pool = if (is.null(dec)) 0L else length(dec$pool)),
            class = "dfdr_check")
}

#' @export
print.dfdr_check <- function(x, ...) {
  u <- x$units
  cat(sprintf("dfdr_check: %d model%s, trained in '%s', decoys matched in %s (pool %d genes)\n", nrow(x$models),
              if (nrow(x$models) == 1) "" else "s", x$train, paste0("'", x$match, "'", collapse = ", "), x$pool))
  drop <- x$models$status != "kept"
  if (any(drop)) cat(sprintf("  %d dropped for lack of measured genes\n", sum(drop)))
  for (v in x$validate) {
    y <- u[u$cohort == v, ]
    cat(sprintf("\nValidation cohort '%s': %d of %d scorable models validated (log-rank P < %.2g, high group worse)\n",
                v, sum(y$validated, na.rm = TRUE), sum(y$scorable), x$alpha))
    f <- tryCatch(cohort_fvr(v), error = function(e) NULL)
    if (!is.null(f) && nrow(f) == 1)
      cat(sprintf("  this cohort (%s) in the paper: %.0f%% of null models fitted in TCGA %s validated here\n",
                  f$accession, 100 * f$fvr, f$type))
    show <- data.frame(model = y$model, genes = paste0(y$measured, "/", y$genes),
                       validated = ifelse(is.na(y$validated), "-", ifelse(y$validated == 1, "yes", "no")),
                       C = sprintf("%.3f", y$cdir), P = formatC(y$p, digits = 2, format = "g"),
                       decoys_validated = ifelse(is.na(y$d_coex), "-", sprintf("%.0f%% of %d", 100 * y$d_coex, y$n_coex)),
                       rho = ifelse(is.na(y$rho), "-", sprintf("%.2f", y$rho)))
    if (nrow(show) > 20) { print(show[1:20, ], row.names = FALSE); cat(sprintf("  ... %d more\n", nrow(show) - 20)) }
    else print(show, row.names = FALSE)
  }
  if (nrow(u) && any(!is.na(u$d_coex)))
    cat("\ndecoys_validated: share of co-expressed decoys passing the same test; rho: share of them the model beats.",
        "\nSee ?check_models for how to read them, and validation_fdr() for the estimate over a set of models.\n")
  invisible(x)
}

#' Decoy-estimated false discovery rate of validation over a set of models
#'
#' The estimator of the accompanying paper: the decoys' validation shares summed over the models, divided by the
#' number of models validated, in one validation cohort or pooled over several. Only models with decoys and scorable
#' in the cohort count. The interval resamples models (within `strata` when given) with replacement.
#'
#' @param x a `dfdr_check`.
#' @param cohort validation cohort(s) to pool; default all.
#' @param models models to include; default all.
#' @param decoys `"coex"` (co-expressed decoys, the paper's primary estimate) or `"rand"`.
#' @param B bootstrap draws (0 for none).
#' @param strata optional vector of strata named by model (the paper resampled within cancer type).
#' @param seed seed of the bootstrap; the user's random number state is left as it was.
#' @return A one-row data frame: `units` (model-cohort pairs), `validated`, `validated_share`, `decoy_share`, `fdr`,
#'   and the 2.5% and 97.5% bootstrap quantiles `fdr_lo`, `fdr_hi`.
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy_survival()
#' chk <- check_models(toy$models, toy$train, toy$validate, n_coex = 8, n_rand = 2, seed = 1)
#' validation_fdr(chk, B = 200)
#' @export
validation_fdr <- function(x, cohort = NULL, models = NULL, decoys = c("coex", "rand"), B = 2000L, strata = NULL,
                           seed = 1L) {
  stopifnot(inherits(x, "dfdr_check"))
  decoys <- match.arg(decoys)
  dcol <- paste0("d_", decoys)
  u <- x$units[x$units$scorable & !is.na(x$units[[dcol]]), ]
  if (!is.null(cohort)) u <- u[u$cohort %in% cohort, ]
  if (!is.null(models)) u <- u[u$model %in% models, ]
  if (!nrow(u)) stop("no scorable model with decoys")
  est <- function(v) c(V = sum(u$validated[v]), E = sum(u[[dcol]][v]), M = length(v))
  e <- est(seq_len(nrow(u)))
  lo <- hi <- NA_real_
  if (B > 0) {
    old <- rng_save(); on.exit(rng_restore(old), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
    grp <- if (is.null(strata)) rep(1L, nrow(u)) else {
      g <- strata[u$model]; if (anyNA(g)) stop("strata lacks a model"); g }
    idx <- split(seq_len(nrow(u)), grp)
    f <- vapply(seq_len(B), function(b) {
      v <- unlist(lapply(idx, function(i) i[sample.int(length(i), length(i), replace = TRUE)]), use.names = FALSE)
      r <- est(v); if (r[["V"]] > 0) r[["E"]] / r[["V"]] else NA_real_ }, 1.0)
    q <- stats::quantile(f, c(0.025, 0.975), na.rm = TRUE, names = FALSE); lo <- q[1]; hi <- q[2]
  }
  data.frame(decoys = decoys, units = e[["M"]], validated = e[["V"]], validated_share = e[["V"]] / e[["M"]],
             decoy_share = e[["E"]] / e[["M"]], fdr = if (e[["V"]] > 0) e[["E"]] / e[["V"]] else NA_real_,
             fdr_lo = lo, fdr_hi = hi)
}

#' False-validation rates of the paper's validation cohorts
#'
#' For each validation cohort built in the accompanying paper, the share of null models that validated
#' there under rule (L). A null model is the usual pipeline (univariate Cox screen at P < 0.05, then LASSO Cox) run in
#' TCGA of the cohort's cancer type on a decoy pool: genes matched in number, and for co-expressed pools in
#' co-expression, to one of 15 published theme gene sets (ferroptosis, hypoxia and others), drawn from genes in no theme
#' set and no published model of the type. The table also gives, for published models of the type, the validation rate among those whose paper did not
#' name the cohort and among those whose paper did, and their co-expressed decoys' rate. A cohort sharing patients with
#' another carries that cohort's rates (`note`).
#'
#' @param cohort optional: cohort names or accessions. A name of the table (`"LUAD_GSE41271"`) selects its row;
#'   otherwise each is matched by the accession it contains (`"GSE39582"`, `"my_GSE39582_subset"`). GSE41271 holds
#'   two cohorts, its adenocarcinomas and its squamous carcinomas.
#' @param type optional: TCGA cancer type codes.
#' @return A data frame, one row per cohort; with neither argument, the whole table.
#' @examples
#' cohort_fvr("GSE39582")
#' head(cohort_fvr(type = "BRCA"))
#' @export
cohort_fvr <- function(cohort = NULL, type = NULL) {
  f <- system.file("extdata", "validation_cohorts.csv", package = "decoyfdr", mustWork = TRUE)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (!is.null(type)) d <- d[d$type %in% type, ]
  if (!is.null(cohort)) {
    keys <- toupper(d$accession)
    hit <- unique(unlist(lapply(toupper(cohort), function(s) {
      exact <- which(toupper(d$cohort) == s)                  # a cohort name of the table, e.g. LUAD_GSE41271
      if (length(exact)) return(exact)
      which(vapply(keys, function(k) grepl(paste0("(^|[^A-Z0-9])", k, "($|[^A-Z0-9])"), s), TRUE)) })))
    d <- d[sort(hit), ]
  }
  rownames(d) <- NULL
  d
}

#' A small simulated survival setting for `check_models()` examples
#'
#' Draws from the current random number stream; set a seed first. A training and two validation cohorts share a
#' co-expressed program of genes 1-10 that shortens survival, and a second program of genes 101-250 that does not.
#' Three models: six genes of the survival program, five of the other program, four noise genes. The survival program
#' is small, so the first model's co-expressed decoys come mostly from the other program.
#'
#' @param n samples per cohort.
#' @param n_genes genes (at least 300).
#' @return A list with `train` and `validate` (a list of two `dfdr_cohort` objects) and `models` (a named list of
#'   coefficient vectors).
#' @examples
#' set.seed(1)
#' toy <- dfdr_toy_survival()
#' toy$models
#' @export
dfdr_toy_survival <- function(n = 120L, n_genes = 400L) {
  if (n_genes < 300L) stop("n_genes must be at least 300")
  genes <- sprintf("G%03d", seq_len(n_genes))
  mk <- function(nm) {
    f <- stats::rnorm(n)
    X <- matrix(stats::rnorm(n_genes * n), n_genes, n, dimnames = list(genes, sprintf("%s_%03d", nm, seq_len(n))))
    X[1:10, ] <- X[1:10, ] + matrix(f, 10, n, byrow = TRUE)
    X[101:250, ] <- X[101:250, ] + matrix(stats::rnorm(n), 150, n, byrow = TRUE)
    tm <- stats::rexp(n, exp(0.9 * f) / 5); cens <- stats::runif(n, 2, 12)
    dfdr_cohort(X + 5, data.frame(time = pmin(tm, cens), event = as.integer(tm <= cens), row.names = colnames(X)),
                name = nm)
  }
  co <- lapply(c("train", "valid_1", "valid_2"), mk)
  models <- list(program = stats::setNames(c(0.4, 0.3, 0.5, 0.2, 0.3, 0.4), genes[c(1, 2, 4, 5, 7, 9)]),
                 other = stats::setNames(c(0.3, -0.2, 0.4, 0.3, 0.2), genes[c(105, 130, 160, 190, 222)]),
                 noise = stats::setNames(c(0.2, 0.3, -0.4, 0.3), genes[c(271, 280, 301, 333)]))
  list(train = co[[1]], validate = co[2:3], models = models)
}

# ------------------------------------------------------------------ internal

as_models <- function(models) {
  if (is.data.frame(models)) {
    if (!all(c("model", "gene", "coef") %in% names(models))) stop("a model table needs columns model, gene and coef")
    ord <- unique(as.character(models$model))
    models <- lapply(split(models, factor(models$model, levels = ord)),
                     function(d) stats::setNames(as.numeric(d$coef), as.character(d$gene)))
  } else if (is.numeric(models)) models <- list(model = models)
  if (!is.list(models) || is.null(names(models)) || anyDuplicated(names(models)) || any(names(models) == ""))
    stop("models must be a uniquely named list of coefficient vectors")
  for (m in names(models)) {
    x <- models[[m]]
    if (!is.numeric(x) || is.null(names(x)) || any(names(x) == "") || anyNA(names(x)))
      stop("model ", m, ": coefficients must be numeric and named by gene")
    if (any(!is.finite(x))) stop("model ", m, ": coefficients must be finite")
    if (anyDuplicated(names(x))) stop("model ", m, ": a gene appears twice")
    if (length(x) < 2) stop("model ", m, ": a model needs two or more genes")
  }
  models
}

# sum of coefficient x within-cohort z-score over the genes measured
lin_score <- function(co, genes, w) {
  k <- genes %in% rownames(co$Z)
  as.numeric(crossprod(co$Z[genes[k], , drop = FALSE], w[k]))
}

# rule (L): split at the median (ties to the lower group), log-rank P, and whether the high group fares worse
km_valid <- function(s, time, event) {
  g <- as.integer(s > stats::median(s))
  if (length(unique(g)) < 2) return(c(p = NA_real_, worse = NA_real_))
  d <- survival::survdiff(survival::Surv(time, event) ~ g)
  c(p = 1 - stats::pchisq(d$chisq, 1), worse = as.numeric(d$obs[2] > d$exp[2]))
}

rng_save <- function() list(kind = RNGkind(), seed = if (exists(".Random.seed", globalenv(), inherits = FALSE))
  get(".Random.seed", globalenv(), inherits = FALSE) else NULL)

rng_restore <- function(s) {
  RNGkind(s$kind[1], s$kind[2], s$kind[3])
  if (is.null(s$seed)) { if (exists(".Random.seed", globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()) }
  else assign(".Random.seed", s$seed, envir = globalenv())
}
