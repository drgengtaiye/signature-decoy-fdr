"""Data for Figure 1h: ROC curves of Chemokine (Messina) in the Gide cohort, of its 32 co-expressed decoys and of its
1,000-draw null bank (pointwise 5th, 50th and 95th percentiles), from the score store (first decoy set, observed labels).
    python figprep_roc.py   -> icb/v13_fig/roc_example.csv
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/decoys")
from v10_calib_lib import Cohort, D

T, CO = "Chemokine_Messina", "ICB_Gide"
c = Cohort(CO)
y = c.y.astype(bool)
grid = np.linspace(0, 1, 101)


def roc(r):
    """TPR on the FPR grid for a score whose ranks are r (higher rank, higher score); response = positive class."""
    o = np.argsort(-r, kind="stable")
    yy = y[o]
    tpr = np.r_[0, np.cumsum(yy) / yy.sum()]; fpr = np.r_[0, np.cumsum(~yy) / (~yy).sum()]
    return np.interp(grid, fpr, tpr)


ti = c.targets.index(T)
rows = {"fpr": grid, "target": roc(c.RT[ti])}
# the folded statistic reads a score in its better direction: flip the decoys and bank draws whose AUC is below 0.5
def folded(r):
    v = roc(r); auc = np.trapezoid(v, grid)
    return v if auc >= 0.5 else roc(-r)


rows["target"] = folded(c.RT[ti])
info = pd.read_csv(f"{D}/v10_store/decoy_info.csv").set_index("decoy")
di = [i for i, t in enumerate(c.target_of) if t == T and info.loc[c.decoys[i], "pool"] == "A" and info.loc[c.decoys[i], "structure"] == "coex"]
D_ = np.array([folded(c.RD[i]) for i in di])
B_ = np.array([folded(c.RB[ti * c.ndraw + k]) for k in range(c.ndraw)])
for j in range(len(di)):
    rows[f"decoy{j + 1}"] = D_[j]
rows["bank_p05"], rows["bank_p50"], rows["bank_p95"] = np.percentile(B_, [5, 50, 95], axis=0)
out = pd.DataFrame(rows)
out.to_csv(f"{D}/v13_fig/roc_example.csv", index=False)
print(len(di), "decoys; AUC target", round(np.trapezoid(out.target, grid), 3))
