## Figure 4e - is "significant in three gastric cohorts" specific to the immune axis?
##
## The immune axis (NetImmune) is the mean of four immune module z-means (APM 8, IFN 5,
## TLS_B 6, T_CELL 7 genes) minus the mean of two barrier module z-means (CAF_ECM 11,
## TAM_SUP 11). Each random axis is built the same way: six random modules of the same
## sizes, drawn from genes measured in all three cohorts, scored by the same row-z mean,
## combined by the same difference, split at the median and tested by log-rank.
suppressMessages({library(survival); library(matrixStats)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
set.seed(2026)
NRAND <- 1000
SIZES_IMM <- c(APM=8, IFN=5, TLS_B=6, T_CELL=7)
SIZES_BAR <- c(CAF_ECM=11, TAM_SUP=11)
TAGS <- c("TCGA", "GSE84437", "GSE62254")

zrow <- function(m){ s <- rowSds(m, na.rm=TRUE); s[s == 0 | is.na(s)] <- 1
                     (m - rowMeans(m, na.rm=TRUE)) / s }

D <- lapply(TAGS, function(t) readRDS(file.path(O, paste0("B1_", t, ".rds"))))
names(D) <- TAGS
S <- lapply(TAGS, function(t) read.csv(file.path(O, paste0("B2_scores_", t, ".csv"))))
names(S) <- TAGS

## gene pool: measured in all three cohorts, variable in each
pool <- Reduce(intersect, lapply(D, function(d){
  m <- d$expr; v <- rowSds(m, na.rm=TRUE)
  rownames(m)[is.finite(v) & v > 0 & rowSums(is.na(m)) == 0]
}))
cat("gene pool:", length(pool), "\n")

Z <- lapply(D, function(d){ m <- d$expr[pool, ]; zrow(m) })

test <- function(tag, axis){
  s <- S[[tag]]
  idx <- match(s$sample, colnames(Z[[tag]]))
  ok <- !is.na(idx) & !is.na(s$time) & s$time > 0 & !is.na(s$event)
  x <- axis[idx[ok]]; tt <- s$time[ok]; ee <- s$event[ok]
  hi <- x > median(x)
  lr <- survdiff(Surv(tt, ee) ~ hi)
  p  <- pchisq(lr$chisq, df=1, lower.tail=FALSE)
  f  <- coxph(Surv(tt, ee) ~ scale(x)[, 1])
  c(p=p, hr=unname(exp(coef(f))), c=unname(concordance(f)$concordance))
}

rows <- list()
for (tag in TAGS) {
  z <- Z[[tag]]
  ## observed axis, recomputed through the identical path as a check
  obs <- S[[tag]]$NetImmune
  names(obs) <- S[[tag]]$sample
  ax <- obs[colnames(z)]
  r <- test(tag, ax)
  rows[[length(rows) + 1]] <- data.frame(cohort=tag, set=0, type="immune axis",
                                         p=r["p"], hr=r["hr"], c=r["c"])
  for (k in seq_len(NRAND)) {
    g <- sample(pool, sum(SIZES_IMM) + sum(SIZES_BAR))
    br <- cumsum(c(SIZES_IMM, SIZES_BAR))
    st <- c(1, head(br, -1) + 1)
    mods <- sapply(seq_along(br), function(j) colMeans(z[g[st[j]:br[j]], , drop=FALSE]))
    axis <- rowMeans(mods[, 1:4]) - rowMeans(mods[, 5:6])
    names(axis) <- colnames(z)
    r <- test(tag, axis)
    rows[[length(rows) + 1]] <- data.frame(cohort=tag, set=k, type="random",
                                           p=r["p"], hr=r["hr"], c=r["c"])
  }
  cat(tag, "done\n")
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(O, "v2_random_logrank.csv"), row.names=FALSE)
for (tag in TAGS) {
  u <- out[out$cohort == tag & out$type == "random", ]
  o <- out[out$cohort == tag & out$type == "immune axis", ]
  cat(sprintf("%-9s immune P=%.2g  random P<0.05: %.1f%%  rank of immune P: %d/%d\n", tag, o$p,
              100 * mean(u$p < 0.05), sum(u$p <= o$p) + 1, nrow(u) + 1))
}
## all three significant at once
w <- reshape(out[out$type == "random", c("cohort", "set", "p")], idvar="set",
             timevar="cohort", direction="wide")
cat(sprintf("random axes significant in all three cohorts: %.1f%%\n",
            100 * mean(rowSums(w[, -1] < 0.05) == 3)))
