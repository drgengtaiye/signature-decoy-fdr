# Round 13 part A (validation/PLAN_round13_published_risk_models.md, A3-A5): the eligible models of the screening tables
# validation/round13_search/screen_<TYPE>.csv, with genes and coefficients parsed, printed symbols mapped through Entrez
# Gene (org.Hs.eg.db symbols and aliases) to the symbols measured in both cohorts of the type, and the plan's 80% rule
# applied -> validation/published_models.csv (frozen with its SHA-256 by amendment before any model is scored).
# Only gene identifiers of the cohorts are read here, no expression value or outcome.
suppressMessages(library(org.Hs.eg.db))
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
V <- file.path(G, "signature-decoy-fdr", "validation")
TYPES <- c("STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA")
EXTID <- c(STAD = "STAD_GSE84437", SKCM = "SKCM_GSE65904", LUAD = "LUAD_GSE72094", LUSC = "LUSC_GSE157011",
           COAD = "COAD_GSE39582", CESC = "CESC_GSE44001", SARC = "SARC_GSE21050", HNSC = "HNSC_GSE65858",
           BLCA = "BLCA_GSE13507", OV = "OV_GSE32062", BRCA = "BRCA_GSE96058")
# accessions related to each external cohort in GEO (Series_relation, read 2026-11-02): a paper naming any of them
# anywhere in its text counts as using the external cohort (plan A10, as amended)
RELATED <- list(STAD = c("GSE84426", "GSE84433", "GSE147163"), LUSC = c("GSE157009", "GSE157010"), COAD = "GSE40967",
                BRCA = "GSE81540")
rowids <- function(f) { x <- readLines(gzfile(f)); sub("\t.*$", "", x[-1]) }
sym <- toTable(org.Hs.egSYMBOL); ali <- toTable(org.Hs.egALIAS2EG)
eg_of <- function(s) unique(c(sym$gene_id[sym$symbol == s], ali$gene_id[ali$alias_symbol == s]))
names_of <- function(eg) unique(c(sym$symbol[sym$gene_id %in% eg], ali$alias_symbol[ali$gene_id %in% eg]))
num <- function(x) {
  for (m in c("−", "–", "—", "˗", "‐", "‑", "‒", "﹣", "－"))
    x <- gsub(m, "-", x, fixed = TRUE, useBytes = TRUE)
  x <- gsub("[[:space:]]+", "", x, useBytes = TRUE)
  x <- sub("^[[(](.*)[])]$", "\\1", x, useBytes = TRUE)
  x <- gsub("×10^", "e", x, fixed = TRUE, useBytes = TRUE)
  suppressWarnings(as.numeric(x))
}
out <- list()
for (tp in TYPES) {
  f <- file.path(V, "round13_search", sprintf("screen_%s.csv", tp))
  if (!file.exists(f)) { cat("no screening table for", tp, "\n"); next }
  s <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  s <- s[toupper(as.character(s$eligible)) == "TRUE", ]
  if (!nrow(s)) next
  tc <- rowids(file.path(G, "pancan", "expr", paste0(tp, ".gz")))
  ex <- rowids(file.path(G, "external", "cohorts", paste0(EXTID[[tp]], "_expr.tsv.gz")))
  univ <- intersect(tc, ex)
  gse <- c(sub("^[A-Z]+_", "", EXTID[[tp]]), RELATED[[tp]])   # with its GEO SuperSeries / SubSeries
  for (i in seq_len(nrow(s))) {
    g <- trimws(strsplit(s$genes[i], ";")[[1]]); w <- num(strsplit(s$coefs[i], ";")[[1]])
    why <- ""
    if (length(g) != length(w) || anyNA(w)) why <- "genes and coefficients do not parse to equal lengths"
    cur <- vapply(g, function(x) {
      if (x %in% univ) return(x)
      if (toupper(x) %in% univ) return(toupper(x))
      cand <- intersect(names_of(eg_of(x)), univ)
      if (!length(cand)) cand <- intersect(names_of(eg_of(toupper(x))), univ)
      if (length(cand) == 1) cand else NA_character_ }, "")
    ok <- !is.na(cur) & !duplicated(cur)
    frac <- mean(ok)
    if (why == "" && frac < 0.8) why <- sprintf("%.0f%% of genes measured in both cohorts", 100 * frac)
    if (why == "" && sum(ok) < 2) why <- "fewer than two genes measured"
    vc <- as.character(s$validation_cohorts[i]); vc[is.na(vc)] <- ""
    out[[length(out) + 1]] <- data.frame(
      type = tp, model_id = paste0(tp, "_", s$pmcid[i]), pmcid = s$pmcid[i], pmid = s$pmid[i], year = s$year[i],
      theme = s$theme[i], endpoint = s$endpoint[i], n_genes = length(g), genes_printed = paste(g, collapse = ";"),
      coefs_printed = s$coefs[i], genes_current = if (why == "") paste(cur[ok], collapse = ";") else "",
      coefs_current = if (why == "") paste(w[ok], collapse = ";") else "", n_measured = sum(ok),
      frac_measured = round(frac, 3), unmapped = paste(g[!ok], collapse = ";"),
      uses_external = any(vapply(gse, function(a) grepl(paste0("\\b", a, "\\b"), vc), TRUE)), validation_cohorts = vc,
      scored = why == "", reason_not_scored = why, stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, out)
write.csv(R, file.path(V, "published_models.csv"), row.names = FALSE)
cat(nrow(R), "eligible models;", sum(R$scored), "scored;", sum(R$scored & R$n_measured >= 3), "with three or more genes;",
    sum(R$scored & !R$uses_external), "scored and independent of the external cohort\n")
print(table(R$type, R$scored))
