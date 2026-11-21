#!/usr/bin/env Rscript
# Secondary readings of the applications from their saved calls (no new draws): the estimate with only the decoys
# whose internal correlation lies within 0.05 of their target's (closest-matched, as for the response analysis), and
# the calls by family.  Rscript v11_post.R tcga|breast
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
OUT <- file.path(D, "v11_pkg"); app <- commandArgs(trailingOnly = TRUE)[1]
info <- read.csv(file.path(OUT, sprintf("%s_decoy_info.csv", app)), stringsAsFactors = FALSE)
close_ids <- info$decoy[abs(info$r_achieved_med - info$r_target_med) <= 0.05 & info$class == "coex"]
cat(sprintf("%s: %d of %d co-expressed decoys within 0.05 of their target's r\n", app, length(close_ids),
            sum(info$class == "coex")))
files <- if (app == "tcga") c(OS = "tcga_OS_calls.csv", DSS = "tcga_DSS_calls.csv", PFI = "tcga_PFI_calls.csv") else
  c(pc1 = "breast_pc1_calls.csv", mean = "breast_mean_calls.csv")
out <- list()
for (nm in names(files)) {
  f <- file.path(OUT, files[[nm]]); if (!file.exists(f)) next
  calls <- read.csv(f, stringsAsFactors = FALSE)
  for (st in unique(calls$stat)) {
    filt <- calls$kind == "target" | calls$unit %in% close_ids
    r <- decoy_fdr(calls, st, "coex", decoy_filter = filt)
    out[[length(out) + 1]] <- cbind(set = nm, r, decoys = "closest")
  }
}
R <- do.call(rbind, out); print(R)
write.csv(R, file.path(OUT, sprintf("%s_fdr_close.csv", app)), row.names = FALSE)
if (app == "tcga") {
  fam <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
  calls <- read.csv(file.path(OUT, "tcga_OS_calls.csv"), stringsAsFactors = FALSE)
  t <- calls[calls$kind == "target", ]; t$family <- fam$family[match(t$target, fam$signature)]
  print(with(t, table(stat, family, call)))
  write.csv(t, file.path(OUT, "tcga_OS_target_calls.csv"), row.names = FALSE)
}
