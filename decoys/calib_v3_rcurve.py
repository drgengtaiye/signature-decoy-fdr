"""Decoy call rate against the decoy's own internal co-expression (pool A co-expressed decoys).

    python calib_v3_rcurve.py  -> calib_v3_rcurve.csv
"""
import numpy as np, pandas as pd
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
U = pd.read_csv(D + "calib_v3_units.csv")
I = pd.read_csv(D + "calib_decoys_v3_info.csv").set_index("decoy")
U = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex")].copy()
U["r"] = U.unit.map(I.r_achieved_med)
U["rbin"] = pd.cut(U.r, [0, .2, .3, .4, .5, .6, 1])
rows = []
for (nl, st, b), g in U.groupby(["null", "stat", "rbin"], observed=True):
    rows.append(dict(null=nl, stat=st, rbin=str(b), decoys=len(g), call_rate=100 * g.call.mean(),
                     neg_rate=100 * g.neg.mean(), med_delta=g.med_delta.median()))
R = pd.DataFrame(rows); R.to_csv(D + "calib_v3_rcurve.csv", index=False)
print(R[R["null"].isin(["N0", "N2nat"])].pivot_table(index="rbin", columns=["null", "stat"],
      values=["call_rate", "neg_rate"]).round(1).to_string())
print(R[R["null"].eq("N0")].groupby("rbin").decoys.first().to_string())
