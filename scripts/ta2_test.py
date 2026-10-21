"""Test A of validation/PLAN_three_tests.md: does a cohort's inflamed spread predict its mean effect?

    python ta2_test.py   -> icb/ta_results.csv
"""
import functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
S = pd.read_csv(f"{D}/ta_cohort_spread.csv")
dev = pd.read_csv(f"{D}/native_response_matched_v2.csv"); dev["eff"] = dev.AUC - dev.rand_med
val = pd.read_csv(f"{V}/validation_scored.csv"); val["eff"] = val.AUC - val.rand_med
E = pd.concat([dev.groupby("cohort").eff.mean(), val.groupby("cohort").eff.mean()])
FAM = pd.read_csv(f"{D}/v3_dissociation.csv").set_index("signature").family
dev["fam"] = FAM.reindex(dev.signature).to_numpy(); val["fam"] = FAM.reindex(val.signature).to_numpy()
EI = pd.concat([dev[dev.fam == "Inflamed/IFN"].groupby("cohort").eff.mean(),
                val[val.fam == "Inflamed/IFN"].groupby("cohort").eff.mean()])
S["eff"] = E.reindex(S.cohort).to_numpy()
S["eff_inflamed"] = EI.reindex(S.cohort).to_numpy()
S = S.dropna(subset=["spread", "eff"])
print(f"cohorts {len(S)} ({(S.set=='development').sum()} development, {(S.set=='validation').sum()} validation)")


def one_sided(x, y, lab):
    r = spearmanr(x, y); p = r.pvalue / 2 if r.statistic > 0 else 1 - r.pvalue / 2
    print(f"  {lab:46s} n={len(x):3d}  rho={r.statistic:+.3f}  one-sided P={p:.4g}")
    return r.statistic, p


print("\nPRIMARY")
rho, p = one_sided(S.spread, S.eff, "inflamed spread against mean effect")
print(f"  -> {'POSITIVE' if p < 0.05 else 'negative'} by the pre-specified rule")

print("\nSECONDARY")
RENAL = ["ICB_Braun", "ICB_Miao1"]
one_sided(S.spread[~S.cohort.isin(RENAL)], S.eff[~S.cohort.isin(RENAL)], "renal cohorts removed")
resid = lambda a, b: a - np.polyval(np.polyfit(b, a, 1), b)
ln = np.log(S.n.to_numpy(float))
one_sided(resid(S.spread.to_numpy(float), ln), resid(S.eff.to_numpy(float), ln), "partial, given log cohort size")
d2 = S.dropna(subset=["eff_inflamed"])
one_sided(d2.spread, d2.eff_inflamed, "against the inflamed family's own effect")
print("\nby set")
for s_ in ("development", "validation"):
    d = S[S.set == s_]; one_sided(d.spread, d.eff, f"{s_} only")
S.to_csv(f"{D}/ta_results.csv", index=False)
print("\n" + S[["cohort", "set", "n", "spread", "eff"]].sort_values("spread").round(3).to_string(index=False))
