# Gene lists of the 2,912 target-decoy modules, for Additional file 4.
#   Rscript export_decoy_genes_v3.R  -> calib_decoys_v3_genes.csv
x <- readRDS("/groups/adv2105_gp/yichen/Yi/com/gastric/icb/calib_decoys_v3.rds")$decoys
out <- data.frame(decoy = vapply(x, `[[`, "", "name"),
                  genes_ensembl = vapply(x, function(d) paste(d$genes, collapse = ";"), ""))
write.csv(out, "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/calib_decoys_v3_genes.csv", row.names = FALSE)
cat(nrow(out), "decoys written\n")
