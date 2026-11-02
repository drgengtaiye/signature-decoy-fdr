# Round 13 (validation/PLAN_round13_published_risk_models.md): shared pieces of the published-model audit (part A,
# v19_published_models.R) and of the pipeline run on theme pools (part B, v19_pipeline.R).
source("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/common.R")
suppressMessages({library(survival); library(glmnet)})
OUT <- file.path(G, "icb", "v19_audit"); dir.create(OUT, showWarnings = FALSE)
TYPES <- c("STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA")   # plan order = index
EXTID <- c(STAD = "STAD_GSE84437", SKCM = "SKCM_GSE65904", LUAD = "LUAD_GSE72094", LUSC = "LUSC_GSE157011",
           COAD = "COAD_GSE39582", CESC = "CESC_GSE44001", SARC = "SARC_GSE21050", HNSC = "HNSC_GSE65858",
           BLCA = "BLCA_GSE13507", OV = "OV_GSE32062", BRCA = "BRCA_GSE96058")
MODELS <- file.path(G, "signature-decoy-fdr", "validation", "published_models.csv")   # frozen by amendment

# TCGA primary tumors of the type (overall survival, the benchmark's sample set) and its round-11 external cohort
load_pair <- function(tp) {
  tc <- tcga_cohorts("OS", types = tp)[[tp]]
  id <- EXTID[[tp]]; EXT <- file.path(G, "external", "cohorts")
  e <- as.matrix(read.delim(gzfile(file.path(EXT, paste0(id, "_expr.tsv.gz"))), row.names = 1, check.names = FALSE))
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  cl <- read.delim(file.path(EXT, paste0(id, "_clin.tsv")), stringsAsFactors = FALSE)
  ex <- dfdr_cohort(e[, cl$sample, drop = FALSE], data.frame(time = cl$time, event = cl$event, row.names = cl$sample),
                    name = id)
  list(tcga = tc, ext = ex)
}

# The theme pools of plan B1 (msigdbr 26.1.0 and three stated lists)
theme_pools <- function() {
  m <- msigdbr::msigdbr(species = "Homo sapiens"); s <- split(m$gene_symbol, m$gs_name)
  P <- list(ferroptosis = s$WP_FERROPTOSIS, pyroptosis = s$REACTOME_PYROPTOSIS,
            necroptosis = s$GOBP_NECROPTOTIC_PROCESS, anoikis = s$GOBP_REGULATION_OF_ANOIKIS,
            autophagy = s$KEGG_REGULATION_OF_AUTOPHAGY, mitophagy = s$REACTOME_MITOPHAGY,
            senescence = s$FRIDMAN_SENESCENCE_UP, hypoxia = s$HALLMARK_HYPOXIA, glycolysis = s$HALLMARK_GLYCOLYSIS,
            emt = s$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, inflammatory = s$HALLMARK_INFLAMMATORY_RESPONSE,
            lactate = s$GOBP_LACTATE_METABOLIC_PROCESS,
            m6a = c("METTL3", "METTL14", "METTL16", "WTAP", "VIRMA", "RBM15", "RBM15B", "ZC3H13", "FTO", "ALKBH5",
                    "YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", "IGF2BP1", "IGF2BP2", "IGF2BP3", "HNRNPA2B1",
                    "HNRNPC", "FMR1", "LRPPRC", "ELAVL1"),
            cuproptosis = c("NFE2L2", "NLRP3", "ATP7B", "ATP7A", "SLC31A1", "FDX1", "LIAS", "LIPT1", "LIPT2", "DLD",
                            "DLAT", "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A", "DBT", "GCSH", "DLST"),
            disulfidptosis = c("GYS1", "NDUFS1", "OXSM", "LRPPRC", "NDUFA11", "NUBPL", "NCKAP1", "RPN1", "SLC3A2",
                               "SLC7A11"))
  stopifnot(!any(vapply(P, is.null, TRUE)))
  lapply(P, unique)
}

# The frozen published models of a type: list(sets, weights) on current symbols, as validation/published_models.csv
published_models <- function(tp = NULL) {
  if (!file.exists(MODELS)) return(NULL)
  d <- read.csv(MODELS, stringsAsFactors = FALSE)
  d <- d[as.logical(d$scored), ]
  if (!is.null(tp)) d <- d[d$type == tp, ]
  d
}

# (L) of plan A8: score split at its median (ties to the lower group), log-rank P, and whether the high group fares worse
km_valid <- function(s, time, event) {
  g <- as.integer(s > stats::median(s))
  if (length(unique(g)) < 2) return(c(p = NA_real_, worse = NA_real_))
  d <- survdiff(Surv(time, event) ~ g)
  c(p = 1 - stats::pchisq(d$chisq, 1), worse = as.numeric(d$obs[2] > d$exp[2]))
}

# Linear score of plan A6: sum of coefficient x within-cohort z-score, on the genes measured
lin_score <- function(co, genes, w) {
  k <- genes %in% rownames(co$Z)
  as.numeric(crossprod(co$Z[genes[k], , drop = FALSE], w[k]))
}

# Statistics of one score in the two cohorts of a type; dirn = +1 when higher scores mark shorter survival
score_stats <- function(st, se, P, dirn = 1) {
  ct <- cindex_rows(st, P$tcga$time, P$tcga$event); ce <- cindex_rows(se, P$ext$time, P$ext$event)
  kt <- km_valid(dirn * st, P$tcga$time, P$tcga$event); ke <- km_valid(dirn * se, P$ext$time, P$ext$event)
  c(tcga_cdir = if (dirn > 0) 1 - ct else ct, ext_cdir = if (dirn > 0) 1 - ce else ce,
    tcga_p = kt[["p"]], tcga_worse = kt[["worse"]], ext_p = ke[["p"]], ext_worse = ke[["worse"]],
    validated = as.numeric(isTRUE(ke[["p"]] < 0.05) && isTRUE(ke[["worse"]] == 1)))
}
