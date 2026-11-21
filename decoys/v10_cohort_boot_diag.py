"""Why the cohort bootstrap of analysis 5 sits above the point estimate: inside resamples of the 12 cohorts with
replacement, the call rates of targets and decoys (N0, folded and directed), against the observed ones.
Descriptive check of v10_fdr_intervals.py (d).   python v10_cohort_boot_diag.py -> icb/v10_cohort_boot_diag.csv"""
import glob, numpy as np, pandas as pd
from v10_calib_lib import calls
D = "/path/to/analysis_root/icb"
rng = np.random.default_rng(5)
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{D}/calib_v3/ICB_*_v3.csv"))])
X = X[X["null"] == "N0"]; COH = sorted(X.cohort.unique())
DIRN = X.drop_duplicates("unit").set_index("unit").dirn
rows = []
for st in ("fold", "dir"):
    x = X[(X.kind == "target") | ((X.pool == "A") & (X.structure == "coex"))].copy()
    if st == "dir":
        x = x[(x.kind == "target") | (x.unit.map(DIRN) == x.target.map(DIRN))]
    x["delta"] = x[f"{st}_obs"] - x[f"{st}_med"]
    M = x.pivot_table(index="cohort", columns="unit", values="delta").reindex(COH)
    kind = x.drop_duplicates("unit").set_index("unit").kind.reindex(M.columns).to_numpy()
    M = M.to_numpy(float)
    c0 = calls(np.nan_to_num(M)[:, :, None], np.isfinite(M))[:, 0]
    for lab, idxs in (("observed", [np.arange(12)]),
                      ("with replacement", [rng.integers(0, 12, 12) for _ in range(300)]),
                      ("9 of 12 without replacement", [rng.choice(12, 9, replace=False) for _ in range(300)])):
        tc, dc, nd = [], [], []
        for ix in idxs:
            Mb = M[ix]; c = calls(np.nan_to_num(Mb)[:, :, None], np.isfinite(Mb))[:, 0]
            tc.append(c[kind == "target"].sum()); dc.append(c[kind == "decoy"].mean()); nd.append(len(set(ix)))
        rows.append(dict(stat=st, resampling=lab, mean_distinct_cohorts=np.mean(nd), targets_called=np.mean(tc),
                         decoy_call_rate=100 * np.mean(dc), fdr_ratio_of_means=100 * np.mean(dc) * 52 / np.mean(tc)))
R = pd.DataFrame(rows); R.to_csv(f"{D}/v10_cohort_boot_diag.csv", index=False); print(R.round(2).to_string())
