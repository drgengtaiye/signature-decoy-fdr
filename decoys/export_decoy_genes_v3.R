# Gene lists of the 2,912 target-decoy modules, for Supplementary Data 2.
#   Rscript export_decoy_genes_v3.R  -> calib_decoys_v3_genes.csv
x <- readRDS("/path/to/analysis_root/icb/calib_decoys_v3.rds")$decoys
out <- data.frame(decoy = vapply(x, `[[`, "", "name"),
                  genes_ensembl = vapply(x, function(d) paste(d$genes, collapse = ";"), ""))
write.csv(out, "/path/to/analysis_root/icb/calib_decoys_v3_genes.csv", row.names = FALSE)
cat(nrow(out), "decoys written\n")
