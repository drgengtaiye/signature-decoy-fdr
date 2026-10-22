sim_cohorts <- function(nc = 6, n = 60, p = 400, type = "binary", seed = 4) {
  set.seed(seed)
  genes <- sprintf("G%03d", seq_len(p))
  prog <- genes[1:40]                                     # a co-expressed program
  lapply(seq_len(nc), function(i) {
    f <- rnorm(n)
    X <- matrix(rnorm(p * n), p, n, dimnames = list(genes, sprintf("c%d_s%02d", i, seq_len(n))))
    X[prog, ] <- X[prog, ] + matrix(f, length(prog), n, byrow = TRUE)
    X <- X + 5
    if (type == "binary") {
      y <- rbinom(n, 1, plogis(1.2 * f)); names(y) <- colnames(X)
      dfdr_cohort(X, y, name = sprintf("C%d", i))
    } else {
      tm <- rexp(n, exp(-0.8 * f)); ev <- rbinom(n, 1, 0.8)
      dfdr_cohort(X, data.frame(time = tm, event = ev, row.names = colnames(X)), name = sprintf("C%d", i))
    }
  })
}

test_that("the whole procedure runs and calls the planted program, for both outcome types", {
  for (type in c("binary", "survival")) {
    co <- sim_cohorts(type = type)
    lib <- dfdr_library(list(prog_a = sprintf("G%03d", 1:10), prog_b = sprintf("G%03d", 11:25),
                             noise_a = sprintf("G%03d", 301:312), noise_b = sprintf("G%03d", 320:340)),
                        engine = c("mean", "wmean", "pc1", "mean"), direction = c(prog_a = 1, prog_b = 1, noise_a = 1, noise_b = -1))
    res <- run_decoy_fdr(co, lib, n_coex = 6, n_rand = 2, n_sets = 1, n_draw = 100, keep_draws = 50, seed = 7, B = 100)
    tc <- res$calls[res$calls$kind == "target" & res$calls$stat == "dir", ]
    expect_true(all(tc$call[tc$target %in% c("prog_a", "prog_b")]))
    expect_equal(nrow(res$decoys$info), 4 * 8)
    expect_true(all(is.finite(res$fdr$fdr) | res$fdr$called == 0))
    cal <- calibrate_null(res$fit, B = 20, seed = 8)
    expect_equal(nrow(cal$summary), 2)
  }
})

test_that("decoys match their target's size and are drawn from outside the library", {
  co <- sim_cohorts(nc = 3)
  lib <- dfdr_library(list(a = sprintf("G%03d", 1:12)))
  dec <- grow_decoys(co, lib, n_coex = 4, n_rand = 2, seed = 9)
  expect_true(all(dec$info$k == 12))
  expect_false(any(unlist(lapply(dec$decoys, `[[`, "genes")) %in% lib$sets$a))
  dec2 <- grow_decoys(co, lib, n_coex = 4, n_rand = 2, seed = 9, precompute = FALSE)
  expect_identical(dec$info$genes, dec2$info$genes)
})

test_that("calibration deltas at the identity permutation equal the fitted deltas, both statistics", {
  for (type in c("binary", "survival")) {
    co <- sim_cohorts(type = type, nc = 6)
    lib <- dfdr_library(list(a = sprintf("G%03d", 1:10), b = sprintf("G%03d", 301:315)), engine = "mean",
                        direction = c(a = 1, b = -1))
    dec <- grow_decoys(co, lib, n_coex = 4, n_rand = 2, seed = 3)
    fit <- fit_units(co, lib, dec, n_draw = 60, keep_draws = 60, seed = 4)
    DEC <- dec$decoys
    for (ci in seq_along(co)) {
      pd <- decoyfdr:::perm_deltas(co[[ci]], fit$scores[[co[[ci]]$name]], matrix(seq_len(co[[ci]]$n)),
                                   ifelse(is.na(lib$direction), 1L, lib$direction),
                                   vapply(DEC, `[[`, 1L, "dirn"), vapply(DEC, `[[`, "", "target"))
      U <- fit$units[fit$units$cohort == co[[ci]]$name, ]
      ut <- U[match(pd$tg, U$unit), ]; ud <- U[match(pd$dn, U$unit), ]
      expect_equal(as.numeric(pd$tf), ut$fold_obs - ut$fold_med, tolerance = 1e-12)
      expect_equal(as.numeric(pd$td), ut$dir_obs - ut$dir_med, tolerance = 1e-12)
      expect_equal(as.numeric(pd$df), ud$fold_obs - ud$fold_med, tolerance = 1e-12)
      expect_equal(as.numeric(pd$dd), ud$dir_obs - ud$dir_med, tolerance = 1e-12)
    }
  }
})

test_that("absolute coherence matches the mean absolute correlation of a mixed-sign signature", {
  co <- sim_cohorts(nc = 3)
  for (i in seq_along(co)) co[[i]]$Z[21:40, ] <- -co[[i]]$Z[21:40, ]   # genes 21-40 anticorrelated with 1-20
  lib <- dfdr_library(list(mix = sprintf("G%03d", c(1:10, 21:30))), engine = "pc1")
  ds <- grow_decoys(co, lib, n_coex = 4, n_rand = 0, seed = 2, coherence = "signed")
  da <- grow_decoys(co, lib, n_coex = 4, n_rand = 0, seed = 2, coherence = "absolute")
  expect_gt(median(da$info$r_achieved_med), 0.2)              # random 20-gene sets reach about 0.1 here
  expect_gt(median(da$info$r_target_med), median(ds$info$r_target_med))
})

test_that("set 1 of several decoy sets is the single-set draw, scores and weights included", {
  co <- sim_cohorts(nc = 6)
  lib <- dfdr_library(list(a = sprintf("G%03d", 1:10), w = sprintf("G%03d", 301:312)), engine = c("mean", "wmean"),
                      weights = list(w = c(rep(1, 6), rep(-1, 6))), direction = c(a = 1, w = 1))
  d1 <- grow_decoys(co, lib, n_coex = 4, n_rand = 2, seed = 11)
  d3 <- grow_decoys(co, lib, n_coex = 4, n_rand = 2, n_sets = 3, seed = 11)
  expect_equal(nrow(d3$info), 3 * nrow(d1$info))
  expect_identical(d3$info[d3$info$set == 1, c("decoy", "genes", "dirn")], d1$info[, c("decoy", "genes", "dirn")],
                   ignore_attr = TRUE)
  expect_false(any(d3$info$genes[d3$info$set == 2] %in% d1$info$genes))
  f1 <- fit_units(co, lib, d1, n_draw = 50, seed = 12)
  f3 <- fit_units(co, lib, d3, n_draw = 50, seed = 12)
  u1 <- f1$units; u3 <- f3$units[f3$units$unit %in% u1$unit, ]
  expect_equal(u3[, c("unit", "cohort", "fold_obs", "fold_med", "dir_obs")], u1[, c("unit", "cohort", "fold_obs", "fold_med", "dir_obs")],
               ignore_attr = TRUE)
})

test_that("the estimate pooled over sets is the mean of the per-set estimates, and pi0 behaves", {
  co <- sim_cohorts(nc = 6)
  lib <- dfdr_library(list(prog_a = sprintf("G%03d", 1:10), prog_b = sprintf("G%03d", 11:25),
                           noise_a = sprintf("G%03d", 301:312), noise_b = sprintf("G%03d", 320:340)), engine = "mean")
  res <- run_decoy_fdr(co, lib, n_coex = 6, n_rand = 2, n_sets = 3, n_draw = 100, stats = "fold", seed = 7, B = 50)
  f <- res$fdr[res$fdr$class == "coex", ]
  per <- vapply(1:3, function(s) decoy_fdr(res$calls, "fold", "coex", B = 2,
                                            decoy_filter = is.na(res$calls$set) | res$calls$set == s)$fdr, 1.0)
  if (f$called > 0) expect_equal(f$fdr, mean(per), tolerance = 1e-10)
  expect_equal(f$n_sets, 3)
  expect_true(f$pi0 >= 0 && f$pi0 <= 1)
  if (f$called > 0) expect_equal(f$fdr_pi0, f$pi0 * f$fdr, tolerance = 1e-10)
  # every decoy made a target, its siblings its decoys: all targets are null, pi0 near 1
  dd <- res$calls[res$calls$stat == "fold" & res$calls$kind == "decoy" & res$calls$class == "coex", ]
  fake <- do.call(rbind, lapply(seq_len(nrow(dd)), function(i) {
    sib <- dd[dd$target == dd$target[i] & seq_len(nrow(dd)) != i, ]; sib$target <- dd$unit[i]
    tg <- dd[i, ]; tg$kind <- "target"; tg$class <- "signature"; tg$set <- NA; tg$target <- dd$unit[i]
    rbind(tg, sib) }))
  expect_gt(decoy_fdr(fake, "fold", "coex", B = 2)$pi0, 0.7)
})
