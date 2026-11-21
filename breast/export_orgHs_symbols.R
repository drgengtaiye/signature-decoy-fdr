## Export official-symbol and alias -> Entrez Gene ID tables from the installed org.Hs.eg.db
## (no downloaded data are read).  Usage: Rscript export_orgHs_symbols.R OUT_DIR
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "/path/to/analysis_root/breast/annotation"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(org.Hs.eg.db))
sym <- toTable(org.Hs.egSYMBOL)          # gene_id, symbol
ali <- toTable(org.Hs.egALIAS2EG)        # gene_id, alias_symbol
write.table(sym, file.path(out, "orgHs_symbol2entrez.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ali, file.path(out, "orgHs_alias2entrez.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db"))),
             capture.output(print(org.Hs.eg_dbInfo()))), file.path(out, "orgHs_version.txt"))
cat("symbols", nrow(sym), "aliases", nrow(ali), "\n")
