"""Shared-bank vs own-bank reading of the same 208 decoys (4 pool-A coex per target).

    python calib_v3_ownbank_compare.py -> calib_v3_ownbank_compare.csv
Same observed AUC on both sides; only the null bank differs: the target's 1,000-draw bank, or a
300-draw bank built for the decoy itself by the published procedure (calib_ownbank_v3.R).
"""
import glob, numpy as np, pandas as pd, multiprocessing as mp
from scipy.stats import wilcoxon
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
O = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(D + "calib_v3/ICB_*_own.csv"))])
S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(D + "calib_v3/ICB_*_v3.csv"))])
ON = sorted(glob.glob(D + "calib_v3/ICB_*_ownnat.csv")); SN = sorted(glob.glob(D + "calib_v3/ICB_*_nat.csv"))
if len(ON) == 12 and len(SN) == 12:          # natural-pool N1/N2: own banks vs the targets' natural banks
    O = pd.concat([O, pd.concat([pd.read_csv(f) for f in ON])])
    S = pd.concat([S, pd.concat([pd.read_csv(f) for f in SN])])
S = S[S.unit.isin(O.unit.unique())]
print(f"own-bank cohorts {O.cohort.nunique()} | decoys {O.unit.nunique()}")
def _one(a):
    (u, nl, src), g = a; out = []
    for st in ("fold", "dir"):
        d = (g[st + "_obs"] - g[st + "_med"]).to_numpy()
        p = wilcoxon(d).pvalue if len(d) >= 6 and np.any(d != 0) else 1.0
        out.append(dict(unit=u, null=nl, src=src, stat=st, dirn=g.dirn.iloc[0], med=float(np.median(d)),
                        call=(p < .05) and np.median(d) > 0, clear95=float(np.mean(g[st+"_obs"] > g[st+"_p95"]))))
    return out
A = pd.concat([O.assign(src="own"), S.assign(src="shared")])
with mp.get_context("fork").Pool(32) as pool:
    R = pd.DataFrame([r for rr in pool.map(_one, list(A.groupby(["unit", "null", "src"])), chunksize=20) for r in rr])
R.to_csv(D + "calib_v3_ownbank_compare.csv", index=False)
print(f"\n{'':10s} {'shared bank':>24s} {'own bank':>24s}   per-test p95 (shared/own)   agreement of median delta")
for nl in [n for n in ("N0", "N1", "N2", "N1nat", "N2nat") if n in set(R["null"])]:
    for st in ("fold", "dir"):
        q = R[(R["null"] == nl) & (R.stat == st)].pivot_table(index="unit", columns="src", values=["call", "med", "clear95"])
        dn = R.drop_duplicates("unit").set_index("unit").dirn.reindex(q.index)
        def rate(src):
            c = q[("call", src)]
            return f"{100*c.mean():4.1f}% (+1 {100*c[dn==1].mean():4.1f} / -1 {100*c[dn==-1].mean():4.1f})"
        r = np.corrcoef(q[("med", "shared")], q[("med", "own")])[0, 1]
        print(f"{nl:5s} {st:4s}  {rate('shared'):>24s} {rate('own'):>24s}   "
              f"{100*q[('clear95','shared')].mean():4.1f}% / {100*q[('clear95','own')].mean():4.1f}%            r = {r:.2f}")
