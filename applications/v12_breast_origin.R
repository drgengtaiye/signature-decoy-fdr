#!/usr/bin/env Rscript
# Analysis 5 of validation/PLAN_round7_estimand_and_bias.md: the breast application with each signature's origin cohorts
# left out. The provenance matrix (validation/venet_provenance.csv, frozen with its SHA-256 before this runs) codes each
# of the 45 measured Venet sets against the 13 cohorts: D (outcome used to pick genes or fit weights, cut-offs or a
# model), C (the cohort's outcome prompted or shaped the signature), V (only tested), blank. D and C tests are dropped
# for the target and its decoys; the conservative variant of the plan also drops every borderline entry the reader
# listed (validation/venet_provenance_borderline.csv: tiers 1, 2, 3 and inherited), and an intermediate variant only
# tiers 1 and inherited (the reader's own suggestion, added after the plan). Reads the ten-set fits of v12_rerun.R
# breast; no new draws.
#   DFDR_LIB=<frozen 0.2.0 install> Rscript v12_breast_origin.R
source("/path/to/analysis_root/signature-decoy-fdr/applications/common.R")
VAL <- "/path/to/analysis_root/signature-decoy-fdr/validation"
OUT <- file.path(D, "v12_pkg")
PV <- read.csv(file.path(VAL, "venet_provenance.csv"), check.names = FALSE, stringsAsFactors = FALSE)
BL <- read.csv(file.path(VAL, "venet_provenance_borderline.csv"), stringsAsFactors = FALSE)
COH <- c("NKI", "Loi_GSE6532", "TCGA_BRCA", "GSE7390", "GSE11121", "GSE1456", "GSE4922", "GSE20685", "GSE21653",
         "GSE42568", "GSE17705", "GSE96058", "METABRIC")
long <- do.call(rbind, lapply(COH, function(co) data.frame(target = PV$signature, cohort = co, code = trimws(PV[[co]]))))
ex_dc <- long[long$code %in% c("D", "C"), c("target", "cohort")]
ex_t1 <- unique(rbind(ex_dc, BL[BL$tier %in% c("tier 1", "inherited"), c("target", "cohort")]))
ex_cons <- unique(rbind(ex_dc, BL[, c("target", "cohort")]))
cat(nrow(ex_dc), "D/C tests (", length(unique(ex_dc$target)), "signatures );", nrow(ex_t1), "with tier 1 and inherited;",
    nrow(ex_cons), "with every borderline entry\n")
print(ex_dc)
CLb <- read.csv(file.path(OUT, "breast_clusters.csv"), stringsAsFactors = FALSE); cl <- setNames(CLb$cluster, CLb$signature)
rows <- list()
for (sc in c("pc1", "mean")) {
  U <- read.csv(gzfile(file.path(OUT, if (sc == "pc1") "breast_pc1_units.csv.gz" else "breast_mean_units.csv.gz")),
                stringsAsFactors = FALSE)
  for (rule in c("none", "D and C", "D, C, tier 1 and inherited", "D, C and every borderline")) {
    ex <- switch(rule, "none" = NULL, "D and C" = ex_dc, "D, C, tier 1 and inherited" = ex_t1, "D, C and every borderline" = ex_cons)
    ca <- call_units(U, stats = "fold", exclude = ex)
    for (b in c("targets", "clusters")) {
      f <- decoy_fdr(ca, "fold", "coex", clusters = if (b == "clusters") cl else NULL)
      rows[[length(rows) + 1]] <- cbind(score = sc, rule = rule, boot = b, f)
    }
    if (rule == "D and C") write.csv(ca, gzfile(file.path(OUT, sprintf("breast_%s_origin_calls.csv.gz", sc))), row.names = FALSE)
  }
}
R <- do.call(rbind, rows)
print(R[, c("score", "rule", "boot", "called", "fdr", "fdr_lo", "fdr_hi", "fdr_set_min", "fdr_set_max", "pi0", "fdr_pi0",
            "fdr_pi0_lo", "fdr_pi0_hi")])
write.csv(R, file.path(OUT, "breast_origin_fdr.csv"), row.names = FALSE)
write.csv(ex_dc, file.path(OUT, "breast_origin_excluded_tests.csv"), row.names = FALSE)
