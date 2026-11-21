"""Data for Figure 1b and l-o: two signatures read against their null banks and their decoys in every response cohort
(first decoy set, N0, observed labels), from the score store of decoys/v10_store.R.

    python figprep_example.py   -> icb/v13_fig/example_units.csv
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/decoys")
from v10_calib_lib import Library, D

L = Library()
TARGETS = ["Chemokine_Messina", "CIN70_Carter"]
rows = []
for c in L.C:
    y = c.y[:, None].astype(np.float32)
    for T in TARGETS:
        if T not in c.targets:
            continue
        ti = c.targets.index(T)
        a_t = c.auc(c.RT[[ti]], y)[0, 0]
        a_b = c.auc(c.RB[ti * c.ndraw:(ti + 1) * c.ndraw], y)[:, 0]
        di = [i for i, (d, t) in enumerate(zip(c.decoys, c.target_of)) if t == T]
        info = L.info.loc[[c.decoys[i] for i in di]]
        keep = [i for i, (p, s) in zip(di, zip(info.pool, info.structure)) if p == "A" and s == "coex"]
        a_d = c.auc(c.RD[keep], y)[:, 0]
        for kind, vals, names in (("target", [a_t], [T]), ("bank", a_b, [f"bank{i}" for i in range(len(a_b))]),
                                  ("decoy", a_d, [c.decoys[i] for i in keep])):
            for v, n in zip(vals, names):
                rows.append(dict(cohort=c.name, target=T, kind=kind, unit=n, auc=float(v), fold=float(max(v, 1 - v))))
out = pd.DataFrame(rows)
out.to_csv(f"{D}/v13_fig/example_units.csv", index=False)
print(out.groupby(["target", "kind"]).size())
