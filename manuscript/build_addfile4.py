"""Additional file 4: target-decoy tables.

    python build_addfile4.py <outdir>  -> Additional_file_4_target_decoy.xlsx
"""
import os, sys, pandas as pd
I = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
SHEETS = [
 ("decoys", "The 2,912 decoy modules: target, pool (A, genes in no signature; B, also outside the 90-signature TME "
            "library), structure (coex, co-expressed; rand, uniform random genes), size, engine, declared direction "
            "(1 sensitivity, -1 resistance), internal r of target and decoy (median over cohorts), correlation of "
            "their cohort profiles, and Ensembl gene identifiers",
  pd.read_csv(I + "calib_decoys_v3_info.csv").merge(pd.read_csv(I + "calib_decoys_v3_genes.csv"), on="decoy")),
 ("tests", "Every signature and decoy under each null (N0, N1, N2 from the original pool; N1nat, N2nat from the "
           "natural pool) and statistic (fold, dir): median delta AUC over 12 cohorts, Wilcoxon P, call (P < 0.05, "
           "positive median), significantly negative, share of cohorts above the null 95th percentile",
  pd.read_csv(I + "calib_v3_units.csv")),
 ("fdr", "Target-decoy FDR for each null, statistic and decoy set, with 95% cluster-bootstrap intervals (6,000 replicates)",
  pd.read_csv(I + "calib_v3_fdr.csv")),
 ("fdr_paired", "Paired bootstrap differences between analyses (same resamples of targets)",
  pd.read_csv(I + "calib_v3_fdr_diff.csv")),
 ("fdr_matching", "FDR with only decoys within 0.05 of their target's internal r, and for targets with internal r <= 0.5",
  pd.read_csv(I + "calib_v3_fdr_matched.csv")),
 ("fdr_paired_close", "Paired differences with closest-matched decoys, overall and by target co-expression",
  pd.read_csv(I + "calib_v3_fdr_diff_close.csv")),
 ("decoy_r_curve", "Decoy call and negative rates by the decoy's internal co-expression",
  pd.read_csv(I + "calib_v3_rcurve.csv")),
 ("per_signature", "Each signature's median delta AUC against its 32 co-expressed decoys",
  pd.read_csv(I + "calib_v3_per_signature.csv")),
 ("own_bank", "208 decoys read against their target's null bank and against a bank of their own",
  pd.read_csv(I + "calib_v3_ownbank_compare.csv")),
]
out = os.path.join(sys.argv[1], "Additional_file_4_target_decoy.xlsx")
with pd.ExcelWriter(out) as xw:
    pd.DataFrame([(n, d, len(t)) for n, d, t in SHEETS], columns=["sheet", "content", "rows"]).to_excel(xw, sheet_name="key", index=False)
    for n, _, t in SHEETS: t.to_excel(xw, sheet_name=n, index=False)
print(out, os.path.getsize(out) // 1024, "kB")
