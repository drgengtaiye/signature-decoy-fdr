# Round 15 B2 (amendment 3, decision 7): ICGC OV-AU and CGCI-HTMCP-CC from their open files into the cohort format of
# round 11 (<cohort>_expr.tsv.gz: symbol x sample, log2(x + 1), the feature with the largest mean kept per symbol;
# <cohort>_clin.tsv: sample, time in months, event = death). One primary-tumour specimen per patient (the first by
# identifier when there are several). Downloaded files are read as text only.
suppressMessages({library(org.Hs.eg.db); library(data.table)})
R <- "/path/to/analysis_root/external/raw_r15"
O <- "/path/to/analysis_root/external/cohorts_r15"
write_cohort <- function(id, X, clin) {
  X <- log2(pmax(X, 0) + 1)
  m <- rowMeans(X); X <- X[order(-m), , drop = FALSE]; X <- X[!duplicated(rownames(X)), , drop = FALSE]
  X <- X[rowSums(is.finite(X)) == ncol(X), clin$sample, drop = FALSE]
  fwrite(data.table(symbol = rownames(X), X), file.path(O, paste0(id, "_expr.tsv.gz")), sep = "\t")
  fwrite(clin, file.path(O, paste0(id, "_clin.tsv")), sep = "\t")
  cat(id, ":", nrow(X), "genes,", nrow(clin), "patients,", sum(clin$event), "events\n")
}

# ---- ICGC OV-AU: exp_seq normalized read counts of primary solid tumour specimens; OS from the donor records
hdr <- function(f) strsplit(readLines(gzfile(file.path(R, "ICGC_OV-AU", f)), 1), "\t")[[1]]
rd <- function(dir, h, kind) rbindlist(lapply(grep(paste0("/", kind, "/"), list.files(file.path(R, "ICGC_OV-AU", dir), "part-",
                                                         recursive = TRUE, full.names = TRUE), value = TRUE),
                                         function(f) fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = FALSE, col.names = h,
                                                           colClasses = "character", quote = "")))
don <- rd("donor", hdr("donor_header.tsv.gz"), "donor")
spe <- rd("spec", hdr("specimen_header.tsv.gz"), "specimen")
spe <- spe[spe$specimen_type == "Primary tumour - solid tissue"]
ex <- rd("exp_seq", hdr("exp_seq_header.tsv.gz"), "exp_seq")
ex <- ex[ex$icgc_specimen_id %in% spe$icgc_specimen_id]
first <- ex[, .(spec = sort(unique(icgc_specimen_id))[1]), by = icgc_donor_id]
ex <- ex[paste(ex$icgc_donor_id, ex$icgc_specimen_id) %in% paste(first$icgc_donor_id, first$spec)]
sym <- mapIds(org.Hs.eg.db, unique(sub("\\..*$", "", ex$gene_id)), "SYMBOL", "ENSEMBL", multiVals = "first")
ex$symbol <- sym[sub("\\..*$", "", ex$gene_id)]
ex <- ex[!is.na(ex$symbol)]
W <- dcast(ex, symbol + gene_id ~ icgc_donor_id, value.var = "normalized_read_count", fun.aggregate = function(x) as.numeric(x[1]))
X <- as.matrix(W[, -(1:2)]); rownames(X) <- W$symbol; X[is.na(X)] <- 0
d <- don[don$icgc_donor_id %in% colnames(X)]
clin <- data.frame(sample = d$icgc_donor_id, time = as.numeric(d$donor_survival_time) / 30.4375,
                   event = as.integer(d$donor_vital_status == "deceased"))
clin <- clin[is.finite(clin$time) & clin$time > 0, ]
write_cohort("OV_OV-AU", X, clin)

# ---- CGCI-HTMCP-CC: GDC STAR counts (TPM, unstranded) of primary tumours; OS from vital status and days
fl <- fread(file.path(R, "GDC_HTMCP-CC", "files.tsv"))
fl <- fl[fl[["cases.0.samples.0.sample_type"]] == "Primary Tumor"][order(file_id)]
fl <- fl[!duplicated(fl[["cases.0.case_id"]])]
M <- lapply(fl$file_id, function(f) {
  x <- fread(file.path(R, "GDC_HTMCP-CC", paste0(f, ".tsv")), skip = 1, sep = "\t")
  x <- x[grepl("^ENSG", x$gene_id)]; setNames(x$tpm_unstranded, x$gene_name) })
X <- do.call(cbind, M); colnames(X) <- fl[["cases.0.case_id"]]
cs <- jsonlite::fromJSON(file.path(R, "GDC_HTMCP-CC", "cases.json"))$data$hits
vs <- cs$demographic$vital_status; dd <- as.numeric(cs$demographic$days_to_death)
fu <- vapply(cs$diagnoses, function(d) if (is.null(d) || !nrow(d) || is.null(d$days_to_last_follow_up)) NA_real_ else
  suppressWarnings(max(as.numeric(d$days_to_last_follow_up), na.rm = TRUE)), 1.0)
fu[!is.finite(fu)] <- NA
clin <- data.frame(sample = cs$case_id, time = ifelse(vs == "Dead", dd, fu) / 30.4375, event = as.integer(vs == "Dead"))
clin <- clin[clin$sample %in% colnames(X) & is.finite(clin$time) & clin$time > 0, ]
write_cohort("CESC_HTMCP-CC", X, clin)
