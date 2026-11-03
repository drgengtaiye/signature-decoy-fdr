# Round 15 part A (validation/PLAN_round15_literature_cohorts_controls.md A1-A3, AMENDMENT_1_round15): the table of all
# eligible models, round 13's frozen rows unchanged and the new ones of validation/round15_search/screen_r15_final.csv
# parsed and mapped exactly as v19_freeze_models.R did (org.Hs.eg.db symbols and aliases to the symbols measured in both
# cohorts of the type, the 80% rule), with the amendment's explicit symbol readings -> validation/published_models_r15.csv
# (frozen with its SHA-256 by amendment before any new model is scored). Only gene identifiers of the cohorts are read.
suppressMessages(library(org.Hs.eg.db))
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
V <- file.path(G, "signature-decoy-fdr", "validation")
TYPES <- c("STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA")
EXTID <- c(STAD = "STAD_GSE84437", SKCM = "SKCM_GSE65904", LUAD = "LUAD_GSE72094", LUSC = "LUSC_GSE157011",
           COAD = "COAD_GSE39582", CESC = "CESC_GSE44001", SARC = "SARC_GSE21050", HNSC = "HNSC_GSE65858",
           BLCA = "BLCA_GSE13507", OV = "OV_GSE32062", BRCA = "BRCA_GSE96058")
RELATED <- list(STAD = c("GSE84426", "GSE84433", "GSE147163"), LUSC = c("GSE157009", "GSE157010"), COAD = "GSE40967",
                BRCA = "GSE81540")
# amendment 1, rule 9: printed names whose alias is ambiguous or absent, read as the paper itself describes them
READ_AS <- c("PD-1" = "PDCD1", "NOV" = "CCN3", "HSP27" = "HSPB1", "LOC158830" = "CXorf65")
rowids <- function(f) { x <- readLines(gzfile(f)); sub("\t.*$", "", x[-1]) }
sym <- toTable(org.Hs.egSYMBOL); ali <- toTable(org.Hs.egALIAS2EG)
eg_of <- function(s) unique(c(sym$gene_id[sym$symbol == s], ali$gene_id[ali$alias_symbol == s]))
names_of <- function(eg) unique(c(sym$symbol[sym$gene_id %in% eg], ali$alias_symbol[ali$gene_id %in% eg]))
num <- function(x) {
  for (m in c("−", "–", "—", "˗", "‐", "‑", "‒", "﹣", "－"))
    x <- gsub(m, "-", x, fixed = TRUE, useBytes = TRUE)
  x <- gsub("[[:space:]]+", "", x, useBytes = TRUE)
  x <- gsub("^\\+", "", x)
  x <- sub("^[[(](.*)[])]$", "\\1", x, useBytes = TRUE)
  x <- gsub("(×|x|\\*)10\\^?", "e", x, useBytes = TRUE)
  suppressWarnings(as.numeric(x))
}
old <- read.csv(file.path(V, "published_models.csv"), stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
old$round <- 13
s_all <- read.csv(file.path(V, "round15_search", "screen_r15_final.csv"), stringsAsFactors = FALSE, check.names = FALSE,
                  encoding = "UTF-8")
out <- list()
for (tp in TYPES) {
  s <- s_all[s_all$type == tp & toupper(as.character(s_all$eligible)) == "TRUE", ]
  if (!nrow(s)) next
  tc <- rowids(file.path(G, "pancan", "expr", paste0(tp, ".gz")))
  ex <- rowids(file.path(G, "external", "cohorts", paste0(EXTID[[tp]], "_expr.tsv.gz")))
  univ <- intersect(tc, ex)
  gse <- c(sub("^[A-Z]+_", "", EXTID[[tp]]), RELATED[[tp]])
  for (i in seq_len(nrow(s))) {
    g <- trimws(strsplit(s$genes[i], ";")[[1]]); w <- num(trimws(strsplit(s$coefs[i], ";")[[1]]))
    why <- ""
    if (length(g) != length(w) || anyNA(w)) why <- "genes and coefficients do not parse to equal lengths"
    cur <- vapply(g, function(x) {
      if (x %in% names(READ_AS)) x <- READ_AS[[x]]
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
      scored = why == "", reason_not_scored = why, round = 15, stringsAsFactors = FALSE)
  }
}
R <- rbind(old, do.call(rbind, out))
hr <- s_all$pmcid[s_all$hr_weights == "TRUE"]
R$hr_weights <- paste(R$type, R$pmcid) %in% paste(s_all$type, s_all$pmcid)[s_all$hr_weights == "TRUE"]
stopifnot(!anyDuplicated(R$model_id), !anyDuplicated(R$pmcid))
write.csv(R, file.path(V, "published_models_r15.csv"), row.names = FALSE)
cat(nrow(R), "eligible models (round 13:", sum(R$round == 13), "; round 15:", sum(R$round == 15), ");", sum(R$scored),
    "scored;", sum(R$scored & R$n_measured >= 3), "with three or more genes;",
    sum(R$scored & !R$uses_external & R$n_measured >= 3), "scored, three or more genes, independent of the external cohort\n")
print(table(R$type, R$scored))
print(R[R$round == 15 & !R$scored, c("model_id", "reason_not_scored", "unmapped")])
