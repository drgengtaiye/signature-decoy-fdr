"""Per-signature dissociation table, both sides against a size- and engine-matched null.

    python build_dissociation_v3.py        -> icb/v3_dissociation.csv

Inputs
  native_response_matched_v2.csv   AUC per signature x checkpoint cohort (null_matched_v2.R)
  native_prognosis_matched.csv     C-index per signature x TCGA type x endpoint
                                   (null_matched_prog.R); the dissociation uses overall survival
Definitions (identical to native_dissociation.csv, verified to machine precision)
  d_resp / d_prog   median over cohorts / cancer types of (observed - null median)
  f_resp / f_prog   fraction of cohorts / cancer types above the null 95th percentile
  p_resp            two-sided Wilcoxon signed-rank test of the per-cohort gains
  shift             d_resp - d_prog
Why v3
  v1 (native_dissociation.csv) used a fixed 48-gene mean-z null on both sides; v2 matched
  the response side only. Both lost three signatures (WNT-beta-catenin_Takeuchi, CTL_Jiang,
  Tcell-inflamed_Ayers) to an R parsing error in sig_methods.csv. v3 matches both sides and
  restores the three.
ICB_VanDenEnde is excluded: CD8A, IFNG, CXCL9, GZMB and PTPRC are absent from its matrix.
"""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
R = pd.read_csv(f"{D}/native_response_matched_v2.csv")
R = R[R.cohort != "ICB_VanDenEnde"].copy()
P = pd.read_csv(f"{D}/native_prognosis_matched.csv")
P = P[P.endpoint == "OS"].copy()
meth = pd.read_csv(f"{D}/sig_methods_clean.csv").set_index("signature").method

# families as used throughout the manuscript; the three restored signatures are assigned
# by the same rule - T-cell/cytotoxic/IFN programmes are inflamed (both are listed in
# positive.R's INFLAMED set), a tumour-intrinsic pathway signature is 'other'
fam = pd.read_csv(f"{D}/native_dissociation.csv").set_index("signature").family.to_dict()
fam.update({"Tcell-inflamed_Ayers": "Inflamed/IFN", "CTL_Jiang": "Inflamed/IFN",
            "WNT-beta-catenin_Takeuchi": "Other"})

R["dl"] = R.AUC - R.rand_med
R["clear"] = R.AUC > R.rand_p95
P["dl"] = P.C - P.rand_med
P["clear"] = P.C > P.rand_p95

r = R.groupby("signature").agg(d_resp=("dl", "median"), n_coh=("dl", "size"),
                               f_resp=("clear", "mean"), n_genes=("n_genes", "median"))
r["p_resp"] = R.groupby("signature").dl.apply(lambda s: wilcoxon(s).pvalue)
p = P.groupby("signature").agg(d_prog=("dl", "median"), n_cancer=("dl", "size"),
                               f_prog=("clear", "mean"))
X = r.join(p, how="inner")
X["shift"] = X.d_resp - X.d_prog
X["method"] = [meth.get(s, "") for s in X.index]
X["family"] = [fam.get(s, "Other") for s in X.index]
missing = [s for s in X.index if s not in fam]
X = X.reset_index()[["signature", "d_prog", "f_prog", "n_cancer", "d_resp", "n_coh", "f_resp",
                     "p_resp", "shift", "n_genes", "method", "family"]]
X.to_csv(f"{D}/v3_dissociation.csv", index=False)

print(f"v3_dissociation: {len(X)} signatures; families "
      f"{X.family.value_counts().to_dict()}; unassigned (-> Other): {missing}")
old = pd.read_csv(f"{D}/v2_dissociation.csv").set_index("signature")
j = X.set_index("signature").join(old, rsuffix="_v2", how="inner")
from scipy.stats import spearmanr
for c in ["d_resp", "d_prog", "shift"]:
    print(f"  {c:7s} v3 vs v2 over {len(j)} shared: rho {spearmanr(j[c], j[c + '_v2'])[0]:+.3f}")
