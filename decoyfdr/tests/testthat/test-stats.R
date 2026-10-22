test_that("auc_rows equals the Mann-Whitney AUC, ties one half", {
  set.seed(1)
  S <- matrix(round(rnorm(5 * 40), 1), 5)
  y <- rep(0:1, 20)
  ref <- apply(S, 1, function(v) { a <- v[y == 1]; b <- v[y == 0]
    mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="), na.rm = TRUE) })
  expect_equal(auc_rows(S, y), ref, tolerance = 1e-12)
})

test_that("cindex_rows equals survival::concordance, with ties in time and score", {
  skip_if_not_installed("survival")
  set.seed(2)
  n <- 120
  tm <- sample(1:40, n, replace = TRUE); ev <- rbinom(n, 1, 0.6)
  S <- rbind(rnorm(n), round(rnorm(n)), rep(1:4, length.out = n), -tm + rnorm(n, 0, 5))
  ref <- apply(S, 1, function(x) unname(survival::concordance(survival::Surv(tm, ev) ~ x)$concordance))
  expect_equal(cindex_rows(S, tm, ev), ref, tolerance = 1e-12)
})

test_that("signed_rank_p matches stats::wilcox.test where both are exact", {
  set.seed(3)
  for (n in c(6, 9, 12, 20)) {
    d <- rnorm(n) + 0.3
    expect_equal(signed_rank_p(d), wilcox.test(d, exact = TRUE)$p.value, tolerance = 1e-12)
  }
  expect_equal(signed_rank_p(rep(0, 8)), 1)
})

test_that("signed_rank_p handles zeros and ties by exact enumeration", {
  d <- c(0, 0.1, -0.1, 0.2, 0.2, 0.3, 0.4, -0.05)
  nz <- d[d != 0]; r <- rank(abs(nz)); w <- sum(r[nz > 0])
  pats <- as.matrix(expand.grid(rep(list(0:1), length(nz))))
  W <- pats %*% r
  ref <- min(1, 2 * mean(W <= min(w, sum(r) - w) + 1e-9))
  expect_equal(signed_rank_p(d), ref)
})

test_that("pc1 scores by power iteration agree with the SVD", {
  set.seed(5)
  f <- rnorm(200); X <- matrix(rnorm(120 * 200), 120) + outer(rnorm(120, 1, .3), f)
  X <- (X - rowMeans(X)) / apply(X, 1, sd)
  ref <- svd(X)$v[, 1]; if (cor(ref, colMeans(X)) < 0) ref <- -ref
  got <- decoyfdr:::pc1_scores(X)
  expect_gt(cor(got, ref), 1 - 1e-10)
  expect_equal(rank(got), rank(ref))
})

test_that("pc1 scores agree with the SVD when genes outnumber samples and with a large gene set", {
  set.seed(6)
  for (dims in list(c(300, 80), c(400, 600))) {
    X <- matrix(rnorm(prod(dims)), dims[1]) + outer(rnorm(dims[1], 1, .3), rnorm(dims[2]))
    X <- (X - rowMeans(X)) / apply(X, 1, sd)
    ref <- svd(X)$v[, 1]; if (cor(ref, colMeans(X)) < 0) ref <- -ref
    expect_gt(cor(decoyfdr:::pc1_scores(X), ref), 1 - 1e-8)
  }
})
