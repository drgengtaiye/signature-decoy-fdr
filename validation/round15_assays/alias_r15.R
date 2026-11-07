# Aliases of the assay genes (org.Hs.eg.db), so that a gene a cohort file carries under an older symbol still counts as measured; per cohort, the first alias present is used and only when the current symbol is absent.
suppressMessages(library(org.Hs.eg.db))
F <- read.csv("assays_r15_frozen.csv", stringsAsFactors = FALSE)
g <- unique(F$gene_symbol)
eg <- mapIds(org.Hs.eg.db, g, "ENTREZID", "SYMBOL", multiVals = "first")
cat("symbols without Entrez:", g[is.na(eg)], "\n")
a <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg[!is.na(eg)]), columns = "ALIAS", keytype = "ENTREZID")
a$gene_symbol <- g[match(a$ENTREZID, eg)]
a <- unique(a[a$ALIAS != a$gene_symbol & grepl("^[A-Za-z0-9-]+$", a$ALIAS), c("gene_symbol", "ALIAS")])
amb <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(a$ALIAS), columns = "ENTREZID", keytype = "ALIAS")
n <- tapply(amb$ENTREZID, amb$ALIAS, function(x) length(unique(x)))
a <- a[n[a$ALIAS] == 1, ]                                # an alias shared with another gene is not used
write.csv(a, "assays_alias_r15.csv", row.names = FALSE); cat(nrow(a), "aliases of", length(unique(a$gene_symbol)), "genes\n")
