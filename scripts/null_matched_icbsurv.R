#!/usr/bin/env Rscript
# Survival under checkpoint blockade, against the same size- and engine-matched null.
#
# Same instrument as null_matched_v2.R (response) and null_matched_prog.R (TCGA survival):
# each signature is scored by its publication's own engine and compared with 200 random gene
# sets of the nearest size on the ladder, scored by the same engine in the same cohort.
# Endpoints are the collection's own overall and progression-free survival, and performance is
# the folded concordance index. Eligibility, endpoints and the metric were fixed in
# validation/PLAN_icb_survival.md (SHA-256 905669928b69dd2f…) before this was run.
suppressMessages({library(GSVA); library(matrixStats); library(MultiAssayExperiment)
                  library(SummarizedExperiment); library(survival); library(parallel)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
NRAND <- 200
LADDER <- c(4, 7, 12, 20, 32, 50, 80, 130, 260, 440)
MINN <- 25; MINEV <- 10
NCORE <- as.integer(Sys.getenv("NCORE", "7"))
OUTDIR <- file.path(D, "null_icbsurv_parts"); dir.create(OUTDIR, showWarnings = FALSE)

meta <- read.csv(file.path(D, "sig_methods_clean.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(meta) == 59)
SPEC <- meta$signature[grepl("Specific", meta$method)]
load_sig <- function(nm, key) {
  f <- file.path(D, "SignatureSets-main/data-raw/IO", paste0(nm, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "latin1")
  g <- if (key == "ensembl" && "gene_id" %in% names(d)) sub("\\..*$", "", d$gene_id) else d$gene_name
  w <- if ("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_, nrow(d))
  ok <- !is.na(g) & g != ""
  list(genes = g[ok], weights = w[ok])
}
cidx <- function(v, y) { ok <- is.finite(v); if (length(unique(y[ok, 1])) < 2) return(NA_real_)
  c_ <- as.numeric(concordance(y[ok] ~ I(-as.numeric(v[ok])))$concordance); max(c_, 1 - c_) }
engine <- function(m) if (grepl("Weighted", m)) "wmean" else if (grepl("ssGSEA", m)) "ssgsea" else "gsva"

run_cohort <- function(f) {
  nm <- sub("\\.rds$", "", basename(f))
  if (nm == "ICB_VanDenEnde") return(NULL)                       # excluded, as elsewhere
  part <- file.path(OUTDIR, paste0(nm, ".csv")); if (file.exists(part)) return(nm)
  set.seed(sum(utf8ToInt(nm)))                                   # same seed as the response run
  x <- try(suppressWarnings(readRDS(f)), silent = TRUE); if (inherits(x, "try-error")) return(NULL)
  ex <- try(experiments(x), silent = TRUE); if (inherits(ex, "try-error")) return(NULL)
  an <- names(ex); pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  if (is.na(pick)) pick <- grep("^expr", an, value = TRUE)[1]
  if (is.na(pick)) return(NULL)
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  eps <- list()
  for (ep in c("OS", "PFS")) {
    tc <- paste0("survival_time_", tolower(ep)); ec <- paste0("event_occurred_", tolower(ep))
    if (!all(c(tc, ec) %in% colnames(cd))) next
    tm <- suppressWarnings(as.numeric(cd[[tc]])); ev <- suppressWarnings(as.numeric(cd[[ec]]))
    keep <- intersect(colnames(M), rownames(cd)[is.finite(tm) & tm > 0 & is.finite(ev)])
    if (length(keep) >= MINN && sum(ev[match(keep, rownames(cd))], na.rm = TRUE) >= MINEV)
      eps[[ep]] <- list(keep = keep, y = Surv(tm[match(keep, rownames(cd))], ev[match(keep, rownames(cd))]))
  }
  if (!length(eps)) return(NULL)
  keep_all <- sort(unique(unlist(lapply(eps, `[[`, "keep"))))
  M <- M[, keep_all, drop = FALSE]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  if (nrow(M) < 3000) return(NULL)
  key <- if (mean(grepl("^ENSG", rownames(M))) > 0.5) "ensembl" else "symbol"
  Z <- { s <- rowSds(M); s[s == 0] <- 1; (M - rowMeans(M)) / s }
  pool <- rownames(M)

  obs <- list(); use <- list(); gsva_sets <- list(); ssg_sets <- list()
  for (i in seq_len(nrow(meta))) {
    s_ <- meta$signature[i]; mth <- meta$method[i]
    if (s_ %in% SPEC) next
    S <- load_sig(s_, key); if (is.null(S)) next
    g <- intersect(S$genes, pool); if (length(g) < 3) next
    if (engine(mth) == "wmean") {
      w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g))
      w[is.na(w)] <- 0; if (sum(abs(w)) == 0) next
      obs[[s_]] <- setNames(as.numeric(crossprod(Z[g, , drop = FALSE], w)) / sum(abs(w)), colnames(Z))
    } else if (engine(mth) == "ssgsea") ssg_sets[[s_]] <- g else gsva_sets[[s_]] <- g
    use[[s_]] <- list(k = length(g), eng = engine(mth))
  }
  if (length(gsva_sets)) { rr <- gsva(gsvaParam(M, gsva_sets, minSize = 3), verbose = FALSE)
                           for (n_ in rownames(rr)) obs[[n_]] <- setNames(as.numeric(rr[n_, ]), colnames(rr)) }
  if (length(ssg_sets))  { rr <- gsva(ssgseaParam(M, ssg_sets, minSize = 3), verbose = FALSE)
                           for (n_ in rownames(rr)) obs[[n_]] <- setNames(as.numeric(rr[n_, ]), colnames(rr)) }

  # random-set scores once per engine and rung, then evaluated on each endpoint
  scores <- list()
  for (eng in sort(unique(sapply(use[names(obs)], `[[`, "eng")))) {
    for (k in LADDER) {
      if (k > length(pool)) next
      sets <- lapply(seq_len(NRAND), function(i) sample(pool, k)); names(sets) <- sprintf("r%04d", seq_along(sets))
      scores[[paste(eng, k)]] <- if (eng == "wmean") t(sapply(sets, function(g) colMeans(Z[g, , drop = FALSE])))
        else gsva(if (eng == "ssgsea") ssgseaParam(M, sets, minSize = 3) else gsvaParam(M, sets, minSize = 3), verbose = FALSE)
    }
  }
  rows <- list()
  for (ep in names(eps)) {
    kp <- eps[[ep]]$keep; y <- eps[[ep]]$y
    nulls <- lapply(scores, function(S_) apply(S_[, kp, drop = FALSE], 1, function(z) cidx(z, y)))
    for (s_ in names(obs)) {
      u <- use[[s_]]; k <- LADDER[which.min(abs(LADDER - u$k))]
      v <- nulls[[paste(u$eng, k)]]; if (is.null(v)) next
      rows[[paste(ep, s_)]] <- data.frame(cohort = nm, endpoint = ep, n = length(kp),
        events = sum(y[, 2]), signature = s_, n_genes = u$k, engine = u$eng, null_size = k,
        C = cidx(obs[[s_]][kp], y), rand_med = median(v, na.rm = TRUE),
        rand_p95 = quantile(v, .95, na.rm = TRUE, names = FALSE),
        rand_mean = mean(v, na.rm = TRUE), rand_sd = sd(v, na.rm = TRUE))
    }
  }
  out <- do.call(rbind, rows); write.csv(out, part, row.names = FALSE)
  cat(sprintf("ICBSURV %-18s endpoints=%s n=%s sigs=%d\n", nm, paste(names(eps), collapse = "/"),
              paste(sapply(eps, function(e) length(e$keep)), collapse = "/"), nrow(out))); nm
}

files <- list.files(D, pattern = "^ICB_.*\\.rds$", full.names = TRUE)
invisible(mclapply(files, function(f) tryCatch(run_cohort(f), error = function(e) {
  cat("FAILED", basename(f), conditionMessage(e), "\n"); NULL }), mc.cores = NCORE, mc.preschedule = FALSE))
R <- do.call(rbind, lapply(list.files(OUTDIR, pattern = "\\.csv$", full.names = TRUE), read.csv))
write.csv(R, file.path(D, "native_icbsurv_matched.csv"), row.names = FALSE)
cat("cohorts:", length(unique(R$cohort)), " tests:", nrow(R), "\n")
