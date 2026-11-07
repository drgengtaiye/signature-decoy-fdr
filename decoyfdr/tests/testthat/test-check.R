toy_check <- function(...) {
  set.seed(1)
  toy <- dfdr_toy_survival()
  list(toy = toy, chk = check_models(toy$models, toy$train, toy$validate, n_coex = 8, n_rand = 2, seed = 3, ...))
}

test_that("check_models scores every model and decoy in every validation cohort", {
  x <- toy_check(); chk <- x$chk
  expect_s3_class(chk, "dfdr_check")
  expect_equal(nrow(chk$units), 3 * 2)
  expect_equal(nrow(chk$decoys), 3 * (8 + 2))
  expect_true(all(chk$units$n_coex == 8) && all(chk$units$n_rand == 2))
  expect_equal(nrow(chk$scores), (3 + 30) * 3)                 # training cohort and two validation cohorts
  p <- chk$units[chk$units$model == "program", ]
  expect_true(all(p$validated == 1))
  expect_true(all(p$rho >= 0.9))
  expect_true(all(chk$units$validated[chk$units$model == "other"] == 0))
  d <- chk$scores[chk$scores$kind == "decoy" & chk$scores$cohort == "valid_1" & chk$scores$model == "program" &
                    chk$scores$class == "coex", ]
  expect_equal(p$d_coex[p$cohort == "valid_1"], mean(d$validated))
  ct <- p$cdir[p$cohort == "valid_1"]
  expect_equal(p$rho[p$cohort == "valid_1"], mean((d$cdir < ct) + 0.5 * (d$cdir == ct)))
  expect_output(print(chk), "Validation cohort 'valid_1'")
})

test_that("decoys are drawn outside the models, keep their model's coefficients and size", {
  x <- toy_check(); chk <- x$chk
  genes <- unlist(lapply(x$toy$models, names))
  dg <- unlist(strsplit(chk$decoys$genes, ";"))
  expect_length(intersect(dg, genes), 0)
  for (i in seq_len(nrow(chk$decoys))) {
    m <- chk$decoys$model[i]
    w <- as.numeric(strsplit(chk$decoys$coefs[i], ";")[[1]])
    expect_equal(sort(w), sort(signif(unname(x$toy$models[[m]]), 8)))
  }
})

test_that("the result depends on the seed only, and the caller's random state is left alone", {
  set.seed(1); toy <- dfdr_toy_survival()
  RNGkind("Mersenne-Twister"); set.seed(42); s0 <- .Random.seed
  a <- check_models(toy$models, toy$train, toy$validate, n_coex = 6, n_rand = 2, seed = 5)
  expect_identical(.Random.seed, s0)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  b <- check_models(toy$models, toy$train, toy$validate, n_coex = 6, n_rand = 2, seed = 5)
  expect_identical(a$units, b$units)
  if (.Platform$OS.type == "unix") {
    c2 <- check_models(toy$models, toy$train, toy$validate, n_coex = 6, n_rand = 2, seed = 5, cores = 2)
    expect_identical(a$units, c2$units)
  }
  tab <- do.call(rbind, lapply(names(toy$models), function(m)
    data.frame(model = m, gene = names(toy$models[[m]]), coef = unname(toy$models[[m]]))))
  d <- check_models(tab, toy$train, toy$validate, n_coex = 6, n_rand = 2, seed = 5)
  expect_identical(a$units, d$units)
})

test_that("genes missing from a cohort follow the 80% rules", {
  set.seed(1); toy <- dfdr_toy_survival()
  v1 <- toy$validate[[1]]
  keep <- setdiff(rownames(v1$expr), c("G001", "G002"))       # two of the six genes of 'program'
  v1b <- dfdr_cohort(v1$expr[keep, ], data.frame(time = v1$time, event = v1$event, row.names = v1$samples),
                     name = "valid_short")
  # outside the matching cohorts: 4 of 6 genes (67%) measured, so 'program' is not scorable there
  chk <- check_models(toy$models, toy$train, list(v1b), match = list(toy$train), n_coex = 4, n_rand = 1, seed = 2)
  u <- chk$units[chk$units$model == "program", ]
  expect_false(u$scorable); expect_true(is.na(u$validated))
  # inside them: kept on its 4 measured genes only if 4 >= 80% of the printed count
  chk2 <- check_models(toy$models, toy$train, list(v1b), n_coex = 4, n_rand = 1, seed = 2)
  expect_match(chk2$models$status[chk2$models$model == "program"], "^dropped")
  chk3 <- check_models(toy$models, toy$train, list(v1b), n_printed = c(program = 5, other = 5, noise = 4),
                       n_coex = 4, n_rand = 1, seed = 2)            # 4 of 5 printed genes: 80%, kept
  expect_equal(chk3$models$status[chk3$models$model == "program"], "kept")
  expect_equal(chk3$models$genes_used[chk3$models$model == "program"], 4L)
  expect_equal(chk3$units$measured[chk3$units$model == "program"], 4L)
})

test_that("validation_fdr is the decoys' validation shares over the number validated", {
  chk <- toy_check()$chk
  f <- validation_fdr(chk, B = 0)
  expect_equal(f$fdr, sum(chk$units$d_coex) / sum(chk$units$validated))
  expect_equal(f$units, 6)
  f1 <- validation_fdr(chk, cohort = "valid_1", B = 100, seed = 1)
  expect_true(f1$fdr_lo <= f1$fdr && f1$fdr <= f1$fdr_hi)
  expect_identical(validation_fdr(chk, B = 100, seed = 1), validation_fdr(chk, B = 100, seed = 1))
})

test_that("inputs are checked", {
  set.seed(1); toy <- dfdr_toy_survival(); bin <- dfdr_toy(n_cohorts = 1)$cohorts[[1]]
  expect_error(check_models(toy$models, toy$train, bin), "no survival outcome")
  expect_error(check_models(list(a = c(G001 = 1, G001 = 2)), toy$train, toy$validate), "appears twice")
  expect_error(check_models(list(a = c(1, 2)), toy$train, toy$validate), "named by gene")
  expect_error(check_models(list(a = c(G001 = 1)), toy$train, toy$validate), "two or more genes")
})

test_that("cohort_fvr finds cohorts by name or accession", {
  all <- cohort_fvr()
  expect_equal(nrow(all), 43)
  expect_true(all(all$fvr >= 0 & all$fvr <= 1))
  expect_equal(all$null_validated, round(all$null_models * all$fvr))
  expect_equal(cohort_fvr("GSE39582")$cohort, "COAD_GSE39582")
  expect_equal(cohort_fvr("my_GSE39582_subset")$cohort, "COAD_GSE39582")
  expect_equal(sort(cohort_fvr("GSE41271")$cohort), c("LUAD_GSE41271", "LUSC_GSE41271"))
  expect_equal(cohort_fvr("LUSC_GSE41271")$cohort, "LUSC_GSE41271")
  expect_equal(nrow(cohort_fvr("GSE1")), 0)
  expect_true(all(cohort_fvr(type = "BRCA")$type == "BRCA"))
})
