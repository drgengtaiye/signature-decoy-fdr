# Shared loaders for the decoyfdr runs (validation/PLAN_round6_applications.md; round 7 sets DFDR_LIB to its frozen
# install of decoyfdr 0.2.0).
RUNLIB <- Sys.getenv("DFDR_LIB", "/path/to/analysis_root/icb/v11_pkg/Rlib3")   # a frozen install (Rlib3: calibrate_null directed deltas fixed; pc1 by the Gram matrix): development installs leave running jobs alone
suppressMessages({library(decoyfdr, lib.loc = RUNLIB); library(matrixStats)})
G   <- "/path/to/analysis_root"
D   <- file.path(G, "icb")
P   <- file.path(G, "pancan")

# The 52 immune signatures: genes, weights, engine and declared direction, as in null_matched_v3_coexpr.R
immune_library <- function(key = "symbol") {
  meta <- read.csv(file.path(D, "sig_methods_clean.csv"), stringsAsFactors = FALSE)
  meta$signature <- trimws(meta$signature); meta$method <- trimws(meta$method)
  meta$association <- trimws(meta$association)
  SPEC <- meta$signature[grepl("Specific", meta$method)]
  load_sig <- function(nm) {
    f <- file.path(D, "SignatureSets-main/data-raw/IO", paste0(nm, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "latin1")
    g <- if (key == "ensembl" && "gene_id" %in% names(d)) sub("\\..*$", "", d$gene_id) else d$gene_name
    w <- if ("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_, nrow(d))
    ok <- !is.na(g) & g != ""
    list(genes = g[ok], weights = w[ok])
  }
  engine <- function(m) if (grepl("Weighted", m)) "wmean" else if (grepl("ssGSEA", m)) "ssgsea" else "gsva"
  targets <- sort(unique(read.csv(file.path(D, "null_v3_tidy.csv"))$signature))
  stopifnot(length(targets) == 52)
  sets <- list(); eng <- c(); w <- list(); dirn <- c()
  for (s in targets) {
    S <- load_sig(s); i <- match(s, meta$signature)
    keep <- !duplicated(S$genes)
    sets[[s]] <- S$genes[keep]; eng[s] <- engine(meta$method[i])
    if (eng[s] == "wmean") { ww <- S$weights[keep]; if (all(is.na(ww))) ww <- rep(1, length(ww)); ww[is.na(ww)] <- 0; w[[s]] <- ww }
    dirn[s] <- if (meta$association[i] == "sensitive") 1L else -1L
  }
  allsig <- unique(unlist(lapply(setdiff(meta$signature, SPEC), function(s) load_sig(s)$genes)))
  list(lib = dfdr_library(sets, eng, w, dirn), exclude = allsig)
}

# The 12 response cohorts with the filters of null_matched_v3_coexpr.R
response_cohorts <- function() {
  suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment)})
  COH <- c("ICB_Mariathasan", "ICB_Braun", "ICB_Gide", "ICB_Liu", "ICB_Kim", "ICB_Jung",
           "ICB_Hugo", "ICB_Riaz", "ICB_Van_Allen", "ICB_Miao1", "ICB_Puch", "ICB_Limagne1")
  lapply(COH, function(co) {
    x <- readRDS(file.path(D, paste0(co, ".rds")))
    ex <- experiments(x); an <- names(ex)
    pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
    if (is.na(pick)) pick <- grep("^expr", an, value = TRUE)[1]
    M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
    r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
    keep <- intersect(colnames(M), names(y)[!is.na(y)])
    M <- M[, keep, drop = FALSE]
    rownames(M) <- sub("\\..*$", "", rownames(M))
    dfdr_cohort(M, y[keep], name = co)
  })
}

# TCGA primary tumors with TCGA-CDR endpoints, as null_matched_prog.R
tcga_cohorts <- function(endpoint = "OS", types = NULL) {
  surv <- read.delim(file.path(P, "Survival_SupplementalTable_S1_20171025_xena_sp"), check.names = FALSE,
                     stringsAsFactors = FALSE)
  cancers <- sub("\\.gz$", "", list.files(file.path(P, "expr"), pattern = "\\.gz$"))
  if (!is.null(types)) cancers <- intersect(cancers, types)
  out <- list()
  for (ca in cancers) {
    e <- as.matrix(read.delim(gzfile(file.path(P, "expr", paste0(ca, ".gz"))), row.names = 1, check.names = FALSE))
    colnames(e) <- gsub("\\.", "-", colnames(e)); e <- e[, substr(colnames(e), 14, 15) == "01", drop = FALSE]
    s <- surv[match(colnames(e), surv$sample), ]
    ok <- !is.na(s$OS) & !is.na(s$OS.time) & s$OS.time > 0          # the benchmark's sample set
    e <- e[, ok, drop = FALSE]; s <- s[ok, ]
    e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
    if (ncol(e) < 80 || sum(s$OS) < 25) next
    ev <- suppressWarnings(as.numeric(s[[endpoint]])); tm <- suppressWarnings(as.numeric(s[[paste0(endpoint, ".time")]]))
    k <- which(!is.na(ev) & !is.na(tm) & tm > 0)
    if (length(k) < 80 || sum(ev[k]) < 25) next
    out[[ca]] <- dfdr_cohort(e[, k, drop = FALSE], data.frame(time = tm[k], event = ev[k], row.names = colnames(e)[k]),
                             name = ca)
  }
  out
}

# The breast cohorts prepared by breast/prep_breast_cohorts.py (Entrez identifiers; one endpoint per cohort)
BREAST <- file.path(G, "breast")
breast_cohorts <- function() {
  man <- read.csv(file.path(BREAST, "cohorts", "manifest.csv"), stringsAsFactors = FALSE)
  inc <- man$cohort[as.logical(man$included)]
  lapply(inc, function(id) {
    e <- read.delim(gzfile(file.path(BREAST, "cohorts", paste0(id, "_expr.tsv.gz"))), check.names = FALSE)
    M <- as.matrix(e[, -1, drop = FALSE]); rownames(M) <- as.character(e$entrez)
    cl <- read.delim(file.path(BREAST, "cohorts", paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
    dfdr_cohort(M, data.frame(time = cl$time_months, event = cl$event, row.names = cl$sample), name = id)
  })
}
breast_sets <- function() {
  s <- read.csv(file.path(BREAST, "library", "venet_cancer_signatures.csv"), stringsAsFactors = FALSE)
  lapply(split(as.character(s$entrez), s$signature), unique)
}
