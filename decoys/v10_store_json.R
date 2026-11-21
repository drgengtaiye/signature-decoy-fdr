#!/usr/bin/env Rscript
# Convert the score-store metadata (v10_store.R) to JSON for v10_calib_sim.py, and add the decoy table.
#   Rscript v10_store_json.R  -> v10_store/<cohort>_meta.json, v10_store/decoy_info.csv, v10_store/targets.json
suppressMessages(library(jsonlite))
D <- "/path/to/analysis_root/icb"; ST <- file.path(D, "v10_store")
for (f in list.files(ST, pattern="_meta\\.rds$", full.names=TRUE)) {
  m <- readRDS(f)
  m$tgt_genes <- lapply(m$tgt_genes, as.character); m$dcy_genes <- lapply(m$dcy_genes, as.character)
  write_json(m, sub("\\.rds$", ".json", f), auto_unbox=TRUE, digits=NA)
}
DV <- readRDS(file.path(D, "calib_decoys_v3.rds"))
I <- DV$info; I$genes <- vapply(DV$decoys[I$decoy], function(d) paste(d$genes, collapse=";"), "")
write.csv(I, file.path(ST, "decoy_info.csv"), row.names=FALSE)
write_json(lapply(DV$targets, function(t) list(sig=t$sig, k=t$k, engine=t$engine, rt=as.list(t$rt))),
           file.path(ST, "targets.json"), auto_unbox=TRUE, digits=NA)
cat("done\n")
