# Round 15 (validation/PLAN_round15_literature_cohorts_controls.md, AMENDMENT_1): the shared pieces of round 13
# (v19_common.R) on the frozen table of all eligible models, outputs in icb/v20_round15/.
source("/path/to/analysis_root/signature-decoy-fdr/applications/v19_common.R")
OUT <- file.path(G, "icb", "v20_round15"); dir.create(OUT, showWarnings = FALSE)
MODELS <- file.path(G, "signature-decoy-fdr", "validation", "published_models_r15.csv")   # frozen by amendment 1
