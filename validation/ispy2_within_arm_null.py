"""I-SPY2, arm by arm, against the same matched random-gene-set null.

The companion of pancan/v8_within_arm_null.py for the trial whose endpoint is pathological
complete response: within each arm, every signature's folded AUC is compared with random gene
sets of its own size drawn from that matrix and scored the same way (mean of gene-wise
z-scores), in the same patients.

    python ispy2_within_arm_null.py   -> validation/ispy2_within_arm_null.csv
"""
import warnings, functools, sys
warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
LADDER = [4, 7, 12, 20, 32, 50, 80, 130, 260, 440]
NRAND = 200
rng = np.random.default_rng(8)


def auc(y, v):
    y = np.asarray(y); v = np.asarray(v, float)
    ok = np.isfinite(v); y, v = y[ok], v[ok]
    if len(np.unique(y)) < 2:
        return np.nan
    r = pd.Series(v).rank().values
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    a = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    return max(a, 1 - a) if False else a          # direction kept: higher score -> more pCR


SIGS = L.signatures()
rows = []
for name, arms in (("ISPY2_pembro", {1: "checkpoint", 0: "comparator"}),
                   ("ISPY2_durva", {1: "checkpoint", 0: "comparator"})):
    M = pd.read_csv(f"{V}/cohorts/{name}.expr.tsv.gz", sep="\t", index_col=0)
    C = pd.read_csv(f"{V}/cohorts/{name}.clin.tsv", sep="\t")
    M.columns = [str(c) for c in M.columns]; C["sample"] = C["sample"].astype(str)
    M = M[C["sample"]]
    M = M[np.isfinite(M).all(axis=1) & (M.std(axis=1) > 0)]
    Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
    pool = list(Z.index)
    print(f"\n=== {name}: {len(C)} patients, pCR {int(C.y.sum())}, pool {len(pool)} genes")
    for a_, arm_name in arms.items():
        m = (C.arm.values == a_)
        y = C.y.values[m]
        if len(y) < 20 or y.sum() < 5:
            print(f"  {arm_name}: skipped (n={len(y)}, pCR={int(y.sum())})"); continue
        null = {}
        for k in LADDER:
            if k > len(pool):
                continue
            v = [auc(y, Z.loc[rng.choice(pool, k, replace=False)].values[:, m].mean(0)) for _ in range(NRAND)]
            null[k] = np.array([x for x in v if np.isfinite(x)])
        for nm, gl in SIGS.items():
            g = [x for x in gl if x in Z.index]
            if len(g) < 3:
                continue
            fam = "Immune" if nm in L.IMM else "Control" if nm.startswith("CTRL_") else "Other"
            if fam == "Other":
                continue
            k = min(LADDER, key=lambda x: abs(x - len(g)))
            a = auc(y, Z.loc[g].values[:, m].mean(0)); v = null[k]
            rows.append(dict(trial=name, arm=arm_name, n=int(m.sum()), pcr=int(y.sum()), signature=nm,
                             family=fam, n_genes=len(g), null_size=k, AUC=a,
                             null_med=float(np.median(v)), null_p95=float(np.quantile(v, .95)),
                             null_sd=float(v.std(ddof=1)), z=(a - float(np.median(v))) / float(v.std(ddof=1)),
                             clears=bool(a > np.quantile(v, .95))))
R = pd.DataFrame(rows); R.to_csv(f"{V}/ispy2_within_arm_null.csv", index=False)
print("\n=== clearing the null's 95th percentile (higher score -> more pCR)")
for trial in R.trial.unique():
    for fam in ("Immune", "Control"):
        for arm_name in ("checkpoint", "comparator"):
            h = R[(R.trial == trial) & (R.family == fam) & (R.arm == arm_name)]
            if len(h):
                print(f"  {trial:13s} {fam:8s} {arm_name:11s} n={h.n.iloc[0]:3d} pCR={h.pcr.iloc[0]:3d}  "
                      f"{int(h.clears.sum()):2d} of {len(h):2d}  median AUC {h.AUC.median():.3f}  median z {h.z.median():+.2f}")
