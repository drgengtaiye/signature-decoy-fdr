# Round 17 (validation/PLAN_round17_decoy_size.md): target and decoy sizes per cohort test, and the estimate with mismatched tests dropped.
#   DFDR_LIB=<frozen install> Rscript v22_size_check.R
suppressMessages(library(decoyfdr, lib.loc = Sys.getenv("DFDR_LIB")))
G <- "/path/to/analysis_root/icb"
P <- file.path(G, "v12_pkg"); O <- file.path(G, "v22_size"); dir.create(O, showWarnings = FALSE)
LIBS <- list(response = c("response_units.csv.gz", "response_fdr.csv"), tcga_OS = c("tcga_OS_units.csv.gz", "tcga_OS_fdr.csv"),
             tcga_DSS = c("tcga_DSS_units.csv.gz", "tcga_DSS_fdr.csv"), tcga_PFI = c("tcga_PFI_units.csv.gz", "tcga_PFI_fdr.csv"),
             breast_pc1 = c("breast_pc1_units.csv.gz", "breast_pc1_fdr.csv"))
S1 <- list(); S2 <- list()
for (nm in names(LIBS)) {
  U <- read.csv(gzfile(file.path(P, LIBS[[nm]][1])), stringsAsFactors = FALSE)
  pub <- read.csv(file.path(P, LIBS[[nm]][2]), stringsAsFactors = FALSE)
  pub <- pub[pub$boot == "targets", ]
  Tt <- U[U$kind == "target", c("cohort", "target", "k")]
  Dd <- U[U$kind == "decoy" & U$class == "coex", c("cohort", "target", "k")]
  kd <- aggregate(k ~ cohort + target, Dd, median); names(kd)[3] <- "kd"
  sp <- aggregate(k ~ cohort + target, Dd, function(x) diff(range(x)))                               # decoys of a test can differ where a cohort lost a gene
  M <- merge(Tt, kd); M$ratio <- M$k / M$kd
  S1[[nm]] <- data.frame(library = nm, tests = nrow(M), differ = sum(M$k != M$kd), target_smaller = sum(M$k < M$kd), target_larger = sum(M$k > M$kd),
                         over_10pct = sum(abs(M$k - M$kd) / M$kd > 0.10), ratio_min = min(M$ratio), ratio_q05 = quantile(M$ratio, .05), ratio_q95 = quantile(M$ratio, .95),
                         ratio_max = max(M$ratio), targets_any = length(unique(M$target[M$k != M$kd])),
                         tests_with_mixed_decoy_sizes = sum(sp$k > 0), max_decoy_spread = max(sp$k))
  write.csv(M, gzfile(file.path(O, sprintf("sizes_%s.csv.gz", nm))), row.names = FALSE)
  excl <- list(none = NULL, exact = M[M$k != M$kd, c("target", "cohort")], within10 = M[abs(M$k - M$kd) / M$kd > 0.10, c("target", "cohort")])
  stats_ <- intersect(c("fold", "dir"), unique(pub$stat))
  for (rule in names(excl)) {
    C <- call_units(U, stats = stats_, exclude = excl[[rule]])
    for (st in stats_) for (cl in c("coex", "rand")) {
      f <- decoy_fdr(C, st, cl)
      S2[[length(S2) + 1]] <- data.frame(library = nm, rule = rule, tests_dropped = if (is.null(excl[[rule]])) 0L else nrow(excl[[rule]]), f)
    }
    if (rule == "none") {                                                         # the control
      for (st in stats_) for (cl in c("coex", "rand")) {
        a <- S2[[length(S2) - (length(stats_) * 2) + (match(st, stats_) - 1) * 2 + match(cl, c("coex", "rand"))]]
        b <- pub[pub$stat == st & pub$decoys == cl, ]
        ok <- nrow(b) == 1 && abs(a$fdr - b$fdr) <= 0.1 && a$called == b$called      # amendment 1
        cat(nm, st, cl, "control", if (ok) "OK" else "MISMATCH", round(a$fdr, 3), if (nrow(b)) round(b$fdr, 3), "\n")
        if (!ok && nrow(b)) stop("control failed: ", nm, " ", st, " ", cl)
      }
    }
  }
  cat(nm, "done\n")
}
write.csv(do.call(rbind, S1), file.path(O, "size_counts.csv"), row.names = FALSE)
write.csv(do.call(rbind, S2), file.path(O, "size_sensitivity.csv"), row.names = FALSE)
