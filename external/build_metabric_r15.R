# Round 15 B2 (amendment 3, decision 7): METABRIC from the round-7 breast preparation (cBioPortal Illumina HT-12,
# rows keyed by Entrez Gene) with rows mapped to current symbols (org.Hs.eg.db; the row with the largest mean kept per
# symbol) and overall survival from the cBioPortal patient file (written by the round-15 build: BRCA_METABRIC_clin.tsv).
suppressMessages({library(org.Hs.eg.db); library(data.table)})
B <- "/groups/adv2105_gp/yichen/Yi/com/gastric/breast/cohorts"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/external/cohorts_r15"
x <- fread(cmd = paste("zcat", file.path(B, "METABRIC_expr.tsv.gz")))
s <- mapIds(org.Hs.eg.db, as.character(x$entrez), "SYMBOL", "ENTREZID")
x <- x[!is.na(s)]; sym <- s[!is.na(s)]
X <- as.matrix(x[, -1]); rownames(X) <- sym
X <- X[order(-rowMeans(X, na.rm = TRUE)), , drop = FALSE]; X <- X[!duplicated(rownames(X)), , drop = FALSE]
cl <- fread(file.path(O, "BRCA_METABRIC_clin.tsv"))
X <- X[, cl$sample, drop = FALSE]
fwrite(data.table(symbol = rownames(X), X), file.path(O, "BRCA_METABRIC_expr.tsv.gz"), sep = "\t")
cat("METABRIC:", nrow(X), "symbols,", ncol(X), "patients,", sum(cl$event), "events\n")
