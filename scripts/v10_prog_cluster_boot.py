"""The prognostic gain over the matched null, tested by resampling whole correlation clusters (Methods, "Statistical
tests"; 20,000 resamples, two-sided). The script that produced the published values (P = 0.006, 5 x 10^-4 and 0.046
for OS, DSS and PFI) was not kept; this re-implementation gives 0.006, 0.0009 and 0.050 on the same table
(icb/_pre_v10/native_prognosis_matched.csv) and is the one reported.

    python v10_prog_cluster_boot.py [table]   -> printed
"""
import sys, numpy as np, pandas as pd
D = "/path/to/analysis_root/icb"
P = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else f"{D}/native_prognosis_matched.csv")
CL = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature").cluster_fixed
for ep in ("OS", "DSS", "PFI"):
    Q = P[P.endpoint == ep]; m = (Q.C - Q.rand_med).groupby(Q.signature).median()
    cl = m.index.map(CL); groups = [m[cl == c].to_numpy() for c in pd.unique(cl)]
    rng = np.random.default_rng(0)
    dr = np.array([np.median(np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))]))
                   for _ in range(20000)])
    print(f"{ep}: median gain {m.median():+.4f}; cluster bootstrap two-sided P = {2 * min((dr <= 0).mean(), (dr >= 0).mean()):.4f}")
