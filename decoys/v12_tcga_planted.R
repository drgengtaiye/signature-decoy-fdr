#!/usr/bin/env Rscript
# Analyses 2, 3 and 6 of validation/PLAN_round7_estimand_and_bias.md: planted signal in TCGA overall survival with new
# causal gene sets, read against the call-probability truth.
#   Rscript v12_tcga_planted.R <setting 0-5> [replicates] [M] [cores]   -> icb/v12_truth/tcga_planted_<setting>.csv
# Settings: T1 library (2 targets), T2 outside (1 Hallmark set, its decoy-pool genes), T3 mixed (1 target + 1 Hallmark
# set), each at beta 0.3 and 0.6 (log hazard ratio per SD of the causal score; a higher score, longer survival).
# Every stored row of the round-6 fit (52 targets, 32 co-expressed and 8 random decoys per target, 200 kept null draws
# per target and type) is read by Harrell's C under each simulated outcome; folded deltas against the kept bank's
# folded median; signed-rank calls over the 20 types. Seeds: causal set and outcome 720000 + 1000 x setting +
# replicate, truth draws 730000 + 1000 x setting + replicate, bootstrap 732000 + ... (arbitrary labels).
suppressMessages({library(decoyfdr, lib.loc = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v12_pkg/Rlib5")
  library(matrixStats); library(survival); library(parallel); library(Rcpp)})
args <- commandArgs(trailingOnly = TRUE)
SET <- as.integer(args[1]); NREP <- if (length(args) > 1) as.integer(args[2]) else 100L
M <- if (length(args) > 2) as.integer(args[3]) else 100L; NC <- if (length(args) > 3) as.integer(args[4]) else 1L
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"; OUTD <- file.path(D, "v12_truth"); dir.create(OUTD, FALSE)
SETTINGS <- data.frame(kind = rep(c("T1 library", "T2 outside", "T3 mixed"), each = 2), beta = rep(c(0.3, 0.6), 3))
kind <- SETTINGS$kind[SET + 1]; beta <- SETTINGS$beta[SET + 1]
LAMBDA <- 0.5

sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace Rcpp;
// Harrell C (decoyfdr::cindex_rows conventions) for every unit (column of RK, dense ranks 1..U[j]) under every outcome
// (columns of TM, EV); units with U[j] == 0 give NA.
// [[Rcpp::export]]
NumericMatrix cidx_multi(IntegerMatrix RK, IntegerVector U, NumericMatrix TM, IntegerMatrix EV) {
  const int n = RK.nrow(), m = RK.ncol(), D = TM.ncol();
  NumericMatrix out(m, D);
  std::vector<int> ord(n); std::vector<double> bit(n + 2);
  for (int d = 0; d < D; ++d) {
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      if (TM(a, d) != TM(b, d)) return TM(a, d) > TM(b, d);
      return EV(a, d) < EV(b, d); });
    std::vector<int> gstart; for (int p = 0; p < n; ++p) if (p == 0 || TM(ord[p], d) != TM(ord[p - 1], d)) gstart.push_back(p);
    gstart.push_back(n);
    for (int j = 0; j < m; ++j) {
      const int u = U[j]; if (u <= 0) { out(j, d) = NA_REAL; continue; }
      std::fill(bit.begin(), bit.begin() + u + 1, 0.0);
      const int* rk = &RK(0, j);
      double conc = 0, disc = 0, tie = 0, ins = 0;
      for (size_t g = 0; g + 1 < gstart.size(); ++g) {
        int pos = gstart[g], end = gstart[g + 1];
        for (int p = pos; p < end; ++p) if (EV(ord[p], d) == 0) { for (int i = rk[ord[p]]; i <= u; i += i & (-i)) bit[i] += 1.0; ins += 1; }
        for (int p = pos; p < end; ++p) if (EV(ord[p], d) == 1) {
          int ri = rk[ord[p]]; double less = 0, leq = 0;
          for (int i = ri - 1; i > 0; i -= i & (-i)) less += bit[i];
          for (int i = ri; i > 0; i -= i & (-i)) leq += bit[i];
          conc += ins - leq; disc += less; tie += leq - less; }
        for (int p = pos; p < end; ++p) if (EV(ord[p], d) == 1) { for (int i = rk[ord[p]]; i <= u; i += i & (-i)) bit[i] += 1.0; ins += 1; }
      }
      double den = conc + disc + tie;
      out(j, d) = den > 0 ? (conc + 0.5 * tie) / den : NA_REAL;
    }
  }
  return out;
}')

t0 <- Sys.time()
fit <- readRDS(file.path(D, "v11_pkg", "tcga_OS_fit.rds"))
DEC <- fit$decoys$decoys; targets <- names(fit$lib$sets)
dn_all <- names(DEC); dtarget <- vapply(DEC, `[[`, "", "target"); dclass <- vapply(DEC, `[[`, "", "class")
coex <- dclass == "coex"
CLf <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
CL <- CLf$cluster_fixed[match(targets, CLf$signature)]
# per type: dense ranks of every stored row (targets, decoys, kept draws), the z-scores, the censoring distribution
prep <- lapply(fit$cohorts, function(co) {
  sc <- fit$scores[[co$name]]
  S <- rbind(sc$targets, sc$decoys, do.call(rbind, sc$bank[rownames(sc$targets)]))
  RK <- rowRanks(S, ties.method = "dense"); U <- rowMaxs(RK); U[is.na(U)] <- 0L; RK[is.na(RK)] <- 1L
  storage.mode(RK) <- "integer"
  km <- survfit(Surv(co$time, 1 - co$event) ~ 1)
  pj <- -diff(c(1, km$surv)); ct <- c(km$time, max(co$time)); pc <- c(pj, tail(km$surv, 1))
  list(name = co$name, n = co$n, Z = co$Z, time = co$time, event = co$event, RKt = t(RK), U = as.integer(U),
       tg = rownames(sc$targets), dn = rownames(sc$decoys), nb = vapply(sc$bank[rownames(sc$targets)], nrow, 1L),
       lam = sum(co$event) / sum(co$time), ct = ct, pc = pc / sum(pc))
})
rm(fit); invisible(gc())
# check: under the observed outcome the C of every target and decoy equals the fit's
u0 <- read.csv(gzfile(file.path(D, "v11_pkg", "tcga_OS_units.csv.gz")), stringsAsFactors = FALSE)
for (p in prep[1:3]) {
  C0 <- cidx_multi(p$RKt, p$U, matrix(p$time), matrix(as.integer(p$event)))[, 1]
  ref <- u0[u0$cohort == p$name, ]; ref <- ref$stat[match(c(p$tg, p$dn), ref$unit)]
  stopifnot(isTRUE(all.equal(C0[seq_along(ref)], ref, tolerance = 1e-12)))
}
cat("C check passed;", length(prep), "types; prepared in", format(Sys.time() - t0), "\n")

suppressMessages(library(msigdbr))
H <- msigdbr(species = "Homo sapiens", collection = "H")
pool <- readRDS(file.path(D, "v11_pkg", "tcga_decoys.rds"))$pool
Hs <- lapply(split(H$gene_symbol, H$gs_name), function(g) intersect(unique(g), pool))
Hs <- Hs[vapply(Hs, length, 1L) >= 15]
cat(sprintf("msigdbr %s: %d Hallmark sets with at least 15 decoy-pool genes\n", packageVersion("msigdbr"), length(Hs)))
Sys.setenv(DFDR_LIB = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v12_pkg/Rlib5")
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R", local = TRUE)
LIB <- immune_library(key = "symbol")$lib
stopifnot(identical(names(LIB$sets), targets))

cdf_tab <- lapply(1:25, function(k) decoyfdr:::signed_rank_cdf(k))
# Wilcoxon signed-rank P and median for every column of X (types x columns); exact, as decoyfdr::signed_rank_p
wsr <- function(X) {
  ok <- is.finite(X); nn <- colSums(ok)
  A <- abs(X); A[!ok] <- NA
  R <- colRanks(A, ties.method = "average", preserveShape = TRUE)
  W <- colSums(R * (X > 0 & ok), na.rm = TRUE)
  p <- rep(1, ncol(X)); med <- colMedians(X, na.rm = TRUE)
  fast <- nn >= 6 & colSums(X == 0 & ok) == 0 & !colAnyDups(A)
  Tm <- pmin(W, nn * (nn + 1) / 2 - W)
  for (k in unique(nn[fast])) { j <- which(fast & nn == k); p[j] <- pmin(1, 2 * cdf_tab[[k]][round(Tm[j]) + 1]) }
  for (j in which(!fast & nn >= 6)) p[j] <- decoyfdr:::signed_rank_p1(X[ok[, j], j])
  p[nn < 6] <- 1; med[nn < 6] <- 0
  list(p = p, med = med, call = p < 0.05 & med > 0)
}
colAnyDups <- function(A) apply(A, 2, function(x) anyDuplicated(x[!is.na(x)]) > 0)
local({                                   # wsr against the package's signed-rank P, with NAs, zeros and ties
  set.seed(1); X <- matrix(round(rnorm(20 * 300), 1), 20); X[sample(length(X), 200)] <- NA
  ref <- apply(X, 2, function(x) if (sum(is.finite(x)) >= 6) decoyfdr:::signed_rank_p1(x[is.finite(x)]) else 1)
  stopifnot(isTRUE(all.equal(wsr(X)$p, ref, tolerance = 1e-12)))
})

# folded deltas of targets and decoys (types x units x draws) under outcome matrices (per type: n x draws)
deltas <- function(G, gen) {
  nd <- NULL; Dt <- NULL; Dd <- NULL
  for (ci in seq_along(prep)) {
    p <- prep[[ci]]; o <- gen(p)
    C <- cidx_multi(p$RKt, p$U, o$tm, o$ev); nd <- ncol(C)
    if (is.null(Dt)) { Dt <- array(NA_real_, c(length(prep), length(targets), nd)); Dd <- array(NA_real_, c(length(prep), length(dn_all), nd)) }
    nt <- length(p$tg); ndc <- length(p$dn)
    Fo <- pmax(C, 1 - C)
    off <- nt + ndc + c(0, cumsum(p$nb))
    med <- t(vapply(seq_len(nt), function(i) colMedians(Fo[(off[i] + 1):off[i + 1], , drop = FALSE]), numeric(nd)))
    if (nd == 1) med <- matrix(med, ncol = 1)
    it <- match(p$tg, targets); id <- match(p$dn, dn_all)
    Dt[ci, it, ] <- Fo[seq_len(nt), , drop = FALSE] - med
    Dd[ci, id, ] <- Fo[nt + seq_len(ndc), , drop = FALSE] - med[match(dtarget[p$dn], p$tg), , drop = FALSE]
  }
  list(t = Dt, d = Dd)
}
outcome_gen <- function(G, nd) function(p) {
  gi <- intersect(G, rownames(p$Z))
  b <- if (length(gi) >= 3) beta else 0
  u <- if (length(gi) >= 3) colMeans(p$Z[gi, , drop = FALSE]) else rep(0, p$n)
  u <- if (sd(u) > 0) (u - mean(u)) / sd(u) else u * 0
  h <- p$lam * exp(-b * u) / mean(exp(-b * u))
  Tt <- matrix(rexp(p$n * nd, rate = rep(h, nd)), p$n, nd)
  Cc <- matrix(sample(p$ct, p$n * nd, replace = TRUE, prob = p$pc), p$n, nd)
  list(tm = pmin(Tt, Cc), ev = matrix(as.integer(Tt <= Cc), p$n, nd))
}
calls_of <- function(DL) {
  nd <- dim(DL$t)[3]
  Xt <- matrix(DL$t, dim(DL$t)[1]); Xd <- matrix(DL$d, dim(DL$d)[1])          # types x (units * draws)
  wt <- wsr(Xt); wd <- wsr(Xd)
  sh <- function(x, nu) matrix(x, nu, nd)
  list(tp = sh(wt$p, length(targets)), tm = sh(wt$med, length(targets)), tc = sh(wt$call, length(targets)),
       dp = sh(wd$p, length(dn_all)), dm = sh(wd$med, length(dn_all)), dc = sh(wd$call, length(dn_all)))
}
pit_of <- function(dc) {                                  # per-target mean call of its co-expressed decoys, per draw
  r <- rowsum(dc[coex, , drop = FALSE] * 1, dtarget[coex]) / as.numeric(table(dtarget[coex]))
  r[targets, , drop = FALSE]
}

one <- function(r) {
  set.seed(720000 + 1000 * SET + r)
  G <- switch(kind,
    "T1 library" = unique(unlist(LIB$sets[sample(targets, 2)])),
    "T2 outside" = Hs[[sample(names(Hs), 1)]],
    "T3 mixed"   = unique(c(LIB$sets[[sample(targets, 1)]], Hs[[sample(names(Hs), 1)]])))
  obs <- calls_of(deltas(G, outcome_gen(G, 1L)))
  set.seed(730000 + 1000 * SET + r)
  tru <- calls_of(deltas(G, outcome_gen(G, M)))
  tc <- obs$tc[, 1]; pit <- pit_of(obs$dc)[, 1]; ok <- !is.na(pit)
  s_t <- sign(obs$tm[, 1]) * -log10(pmax(obs$tp[, 1], 1e-300)); s_d <- sign(obs$dm[, 1]) * -log10(pmax(obs$dp[, 1], 1e-300))
  q <- vapply(seq_along(targets), function(i) { x <- s_d[coex & dtarget == targets[i]]
    if (!length(x)) NA_real_ else (1 + sum(x >= s_t[i] - 1e-12)) / (1 + length(x)) }, 1.0)
  Rn <- sum(tc[ok]); E <- sum(pit[ok]); pi0 <- min(1, sum(q[ok] > LAMBDA) / (sum(ok) * (1 - LAMBDA)))
  ft <- rowMeans(tru$tc); kt <- rowSums(tru$tc); pid <- rowMeans(pit_of(tru$dc))
  null_cp <- ok & ft <= pid
  true_bin <- ok & pbinom(kt - 1, M, pmin(1, pid), lower.tail = FALSE) < 0.05
  set.seed(732000 + 1000 * SET + r)
  v <- cbind(ifelse(ok, tc, 0), ifelse(ok, pit, 0), ok, ok & q > LAMBDA)
  bt <- function(Wm) { R <- Wm[, 1]; F <- ifelse(R > 0, Wm[, 2] / R, NA); F0 <- pmin(1, Wm[, 4] / (Wm[, 3] * (1 - LAMBDA))) * F
    c(quantile(F, c(.025, .975), na.rm = TRUE), quantile(F0, c(.025, .975), na.rm = TRUE)) }
  Wt <- t(rmultinom(1000, length(targets), rep(1, length(targets)))) %*% v
  uc <- sort(unique(CL)); Sc <- rowsum(v, CL)[as.character(uc), ]
  Wc <- t(rmultinom(1000, length(uc), rep(1, length(uc)))) %*% Sc
  bi <- c(bt(Wt), bt(Wc))
  data.frame(setting = SET, kind = kind, beta = beta, rep = r, n_causal = length(G), R = Rn, E = E, pi0 = pi0,
             V_cp = sum(tc & null_cp), V_bin = sum(tc & ok & !true_bin), n_true_cp = sum(ok & !null_cp),
             n_true_bin = sum(true_bin), exp_V = sum(ft[null_cp]), pi_null = sum(pid[null_cp]),
             pi_true = sum(pid[ok & !null_cp]), excess_null = sum((pid - ft)[null_cp]),
             F_lo_t = bi[1], F_hi_t = bi[2], F0_lo_t = bi[3], F0_hi_t = bi[4],
             F_lo_c = bi[5], F_hi_c = bi[6], F0_lo_c = bi[7], F0_hi_c = bi[8],
             decoy_rate_obs = mean(obs$dc[coex, 1]), row.names = NULL)
}
t1 <- Sys.time(); x <- one(0L); cat("replicate 0 in", format(Sys.time() - t1), "\n"); print(x)
res <- mclapply(seq_len(NREP - 1L), function(r) tryCatch(one(r), error = function(e) conditionMessage(e)),
                mc.cores = NC, mc.preschedule = FALSE)
bad <- !vapply(res, is.data.frame, TRUE)
if (any(bad)) { cat("failed replicates:", which(bad), "\n"); print(unique(unlist(res[bad]))) }
out <- do.call(rbind, c(list(x), res[!bad]))
write.csv(out, file.path(OUTD, sprintf("tcga_planted_%d.csv", SET)), row.names = FALSE)
cat(kind, beta, ":", nrow(out), "replicates in", format(Sys.time() - t0), "\n")
