"""The dissociation on a scale that is common to both sides.

    python standardised_effects.py      -> icb/v4_standardised_effects.csv

The shift is currently (AUC - null median) minus (C-index - null median). Both are
concordance measures on a 0.5 baseline, but that does not make the two gains
commensurable: one is a binary discrimination in a checkpoint cohort, the other a
censored survival ranking in a TCGA cohort, and the two have different null spreads.
Subtracting one from the other quietly assumes they do.

Here each side is divided by the standard deviation of its own matched null draw:

    z = (observed - null median) / null SD

which is the effect in units of what random gene sets of the same size, scored by the same
engine, achieve in that same cohort. The dissociation is then recomputed on the z scale and
compared with the raw version.
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]

R = pd.read_csv(f"{D}/native_response_matched_v2.csv")
R = R[R.cohort != "ICB_VanDenEnde"].copy()
P = pd.read_csv(f"{D}/native_prognosis_matched.csv")
P = P[P.endpoint == "OS"].copy()
for T in (R, P):
    assert "rand_sd" in T.columns, "rerun the null scripts: rand_sd column missing"

R["raw"] = R.AUC - R.rand_med
R["z"] = R.raw / R.rand_sd
P["raw"] = P.C - P.rand_med
P["z"] = P.raw / P.rand_sd
print(f"response  : {R.signature.nunique()} signatures x {R.cohort.nunique()} cohorts; "
      f"null SD median {R.rand_sd.median():.4f}")
print(f"prognosis : {P.signature.nunique()} signatures x {P.cancer.nunique()} types;   "
      f"null SD median {P.rand_sd.median():.4f}")
print(f"the two null spreads differ by a factor of {R.rand_sd.median()/P.rand_sd.median():.2f}, "
      "which is what the raw subtraction ignores\n")

fam = pd.read_csv(f"{D}/v3_dissociation.csv").set_index("signature").family
X = pd.DataFrame({
    "z_resp": R.groupby("signature").z.median(),
    "z_prog": P.groupby("signature").z.median(),
    "raw_resp": R.groupby("signature").raw.median(),
    "raw_prog": P.groupby("signature").raw.median()}).dropna()
X["z_shift"] = X.z_resp - X.z_prog
X["raw_shift"] = X.raw_resp - X.raw_prog
X["family"] = [fam.get(s, "Other") for s in X.index]
X.index.name = "signature"
X.reset_index().to_csv(f"{D}/v4_standardised_effects.csv", index=False)

print(f"{'':16s} {'n':>3s}  {'z shift':>18s}  {'raw shift':>18s}")
for f in FAMORD:
    g = X[X.family == f]
    pz = wilcoxon(g.z_shift).pvalue
    pr = wilcoxon(g.raw_shift).pvalue
    print(f"{f:16s} {len(g):3d}  {g.z_shift.median():+7.3f} (P={pz:7.1e})  "
          f"{g.raw_shift.median():+7.4f} (P={pr:7.1e})")

i, c = X[X.family == "Inflamed/IFN"], X[X.family == "Confounder"]
pz = mannwhitneyu(i.z_shift, c.z_shift).pvalue
pr = mannwhitneyu(i.raw_shift, c.raw_shift).pvalue
u = mannwhitneyu(i.z_shift, c.z_shift).statistic
print(f"\ninflamed vs confounder shift:  z scale P = {pz:.2e} (rank-biserial "
      f"r = {2*u/(len(i)*len(c))-1:+.2f});  raw scale P = {pr:.2e}")
pz = wilcoxon(X.z_prog, X.z_resp).pvalue
pr = wilcoxon(X.raw_prog, X.raw_resp).pvalue
print(f"paired prognosis vs response:  z scale P = {pz:.2e};  raw scale P = {pr:.2e}")
print(f"\nmedian z: prognosis {X.z_prog.median():+.3f} SD of its null, "
      f"response {X.z_resp.median():+.3f} SD of its null")
print(f"inflamed family: prognosis {i.z_prog.median():+.3f}, response {i.z_resp.median():+.3f}")
print(f"confounder family: prognosis {c.z_prog.median():+.3f}, response {c.z_resp.median():+.3f}")
print("\nsaved v4_standardised_effects.csv")
