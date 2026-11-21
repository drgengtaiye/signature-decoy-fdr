# Data for the survival examples of Figures 6 and 7: a called signature and its median co-expressed decoy (by folded C in
# that cohort) in one cohort, with the patients' survival, from the set-1 fits of the package runs (applications/v11_*.R).
#   Rscript figprep_km.R   -> icb/v13_fig/km_<tag>.csv
suppressMessages(library(decoyfdr, lib.loc = "/path/to/analysis_root/icb/v12_pkg/Rlib6"))
D <- "/path/to/analysis_root/icb"
ex <- list(list(tag = "tcga", fit = "tcga_OS_fit.rds", units = "tcga_OS_units.csv.gz", cohort = "KIRP", target = "CIN70_Carter"),
           list(tag = "breast", fit = "breast_pc1_fit.rds", units = "breast_pc1_units.csv.gz", cohort = "METABRIC",
                target = "SOTIRIOU-GGI"))
for (e in ex) {
  u <- read.csv(gzfile(file.path(D, "v11_pkg", e$units)), stringsAsFactors = FALSE)
  u <- u[u$cohort == e$cohort & u$target == e$target & u$kind == "decoy" & u$class == "coex", ]
  u <- u[order(u$fold_obs), ]; dmed <- u$unit[ceiling(nrow(u) / 2)]
  fit <- readRDS(file.path(D, "v11_pkg", e$fit))
  co <- fit$cohorts[[match(e$cohort, vapply(fit$cohorts, `[[`, "", "name"))]]
  sc <- fit$scores[[e$cohort]]
  out <- data.frame(sample = co$samples, time = co$time, event = co$event, target = sc$targets[e$target, ],
                    decoy = sc$decoys[dmed, ])
  out$decoy_name <- dmed; out$target_name <- e$target; out$cohort <- e$cohort
  write.csv(out, file.path(D, "v13_fig", sprintf("km_%s.csv", e$tag)), row.names = FALSE)
  cat(e$tag, nrow(out), "patients; decoy", dmed, "C target", cindex_rows(rbind(out$target), out$time, out$event),
      "C decoy", cindex_rows(rbind(out$decoy), out$time, out$event), "\n")
  rm(fit); invisible(gc())
}
