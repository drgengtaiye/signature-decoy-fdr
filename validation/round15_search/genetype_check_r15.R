# Round 15 amendment 1, rule 2: NCBI gene type (org.Hs.eg.db) of every printed symbol of the newly eligible models.
suppressMessages(library(org.Hs.eg.db))
g <- readLines("tmp_genes_r15.txt")
q <- toupper(gsub("\\s+", "", g))
sym <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, "SYMBOL"), columns = c("SYMBOL", "GENETYPE"), keytype = "SYMBOL")
al <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, "ALIAS"), columns = c("ALIAS", "SYMBOL", "GENETYPE"), keytype = "ALIAS")
out <- do.call(rbind, lapply(seq_along(g), function(i) {
  s <- sym[toupper(sym$SYMBOL) == q[i], ]
  if (nrow(s)) return(data.frame(printed = g[i], how = "symbol", symbol = s$SYMBOL[1], genetype = s$GENETYPE[1]))
  a <- unique(al[toupper(al$ALIAS) == q[i], c("SYMBOL", "GENETYPE")])
  if (nrow(a) == 1) return(data.frame(printed = g[i], how = "alias", symbol = a$SYMBOL, genetype = a$GENETYPE))
  if (nrow(a) > 1) return(data.frame(printed = g[i], how = "ambiguous alias", symbol = paste(a$SYMBOL, collapse = "|"),
                                     genetype = paste(unique(a$GENETYPE), collapse = "|")))
  data.frame(printed = g[i], how = "not found", symbol = NA, genetype = NA)
}))
write.table(out, "tmp_biotype_r15.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(table(out$how)); print(table(out$genetype, useNA = "ifany"))
