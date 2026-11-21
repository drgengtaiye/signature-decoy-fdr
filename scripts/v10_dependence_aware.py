"""Assemble v5_dependence_aware.csv (the dependence-aware signature-level statistics) from v5_cluster_stats.csv (v5_clusters.py) and
v5_perm_family.csv (v5_perm_family.py).

    python v10_dependence_aware.py [input_dir]   -> icb/v5_dependence_aware.csv

The table was first assembled by hand from those two files; this script writes the same table and, run on
the inputs of the original benchmark (icb/_pre_v10), must reproduce the file made from them.
"""
import sys, numpy as np, pandas as pd

D = "/path/to/analysis_root/icb"
IN = sys.argv[1] if len(sys.argv) > 1 else D
C = pd.read_csv(f"{IN}/v5_cluster_stats.csv")
F = pd.read_csv(f"{IN}/v5_perm_family.csv")
rows = []
for key in ("cluster_fixed", "cluster_dynamic"):
    for _, r in C[C.clustering == key].iterrows():
        rows.append(dict(statistic=r.statistic, method=f"cluster bootstrap ({key}, {r.n_clusters} clusters)",
                         effect=r.observed_effect, naive_P=r.observed_P, lo=r.cluster_boot_lo,
                         hi=r.cluster_boot_hi, P=r.cluster_boot_P_two_sided))
        rows.append(dict(statistic=r.statistic, method=f"one per cluster ({key})",
                         effect=r.onepercluster_median_effect, naive_P=r.observed_P, lo=np.nan, hi=np.nan,
                         P=r.onepercluster_median_P))
for _, r in F.iterrows():
    rows.append(dict(statistic=f"family contrast: {r.statistic}",
                     method=f"within-cohort response-label permutation ({r.n_perm:,})",
                     effect=r.observed, naive_P=np.nan, lo=r.null_p05, hi=r.null_p95, P=r.perm_P))
out = pd.DataFrame(rows)
if IN != D:
    ref = pd.read_csv(f"{IN}/v5_dependence_aware.csv")
    assert list(ref.statistic) == list(out.statistic) and list(ref.method) == list(out.method)
    num = ["effect", "naive_P", "lo", "hi", "P"]
    assert np.allclose(ref[num].to_numpy(float), out[num].to_numpy(float), rtol=1e-12, atol=0, equal_nan=True)
    print(f"reproduces {IN}/v5_dependence_aware.csv")
else:
    out.to_csv(f"{D}/v5_dependence_aware.csv", index=False)
    print("wrote v5_dependence_aware.csv")
