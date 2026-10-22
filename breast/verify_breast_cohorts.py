#!/usr/bin/env python3
"""Format checks on the prepared breast cohorts (no expression-outcome computation of any kind).

    python3 -I verify_breast_cohorts.py [/groups/adv2105_gp/yichen/Yi/com/gastric/breast/cohorts]
"""
import sys

import numpy as np
import pandas as pd

d = sys.argv[1] if len(sys.argv) > 1 else "/groups/adv2105_gp/yichen/Yi/com/gastric/breast/cohorts"
man = pd.read_csv(f"{d}/manifest.csv")
bad = 0
for _, r in man[man.included].iterrows():
    cid = r["cohort"]
    cl = pd.read_csv(f"{d}/{cid}_clin.tsv", sep="\t")
    x = pd.read_csv(f"{d}/{cid}_expr.tsv.gz", sep="\t", index_col=0)
    v = x.values
    checks = {
        "clin columns": list(cl.columns) == ["sample", "time_months", "event", "endpoint"],
        "first column 'entrez'": x.index.name == "entrez",
        "entrez integer, unique": x.index.dtype.kind == "i" and x.index.is_unique,
        "samples == clin samples": list(x.columns) == list(cl["sample"]),
        "samples unique": cl["sample"].is_unique,
        "time > 0, finite": bool((cl.time_months > 0).all() and np.isfinite(cl.time_months).all()),
        "event in {0,1}": set(cl.event.unique()) <= {0, 1},
        "one endpoint": cl.endpoint.nunique() == 1 and cl.endpoint.iloc[0] == r["endpoint"],
        "n matches manifest": len(cl) == r["n_samples"] and int(cl.event.sum()) == r["n_events"],
        "genes match manifest": len(x) == r["n_genes"],
        "log scale (max < 50)": np.nanmax(v) < 50,
        "n>=100, events>=25, genes>=10000": len(cl) >= 100 and cl.event.sum() >= 25 and len(x) >= 10000,
    }
    fin = np.isfinite(v).all(axis=1)
    ok_genes = int((fin & (np.nanstd(v, axis=1) > 0)).sum())
    failed = [k for k, ok in checks.items() if not ok]
    bad += len(failed)
    print(f"{cid:12s} n={len(cl):5d} events={int(cl.event.sum()):4d} {cl.endpoint.iloc[0]:4s} genes={len(x):6d} "
          f"finite&nonconstant={ok_genes:6d} range=[{np.nanmin(v):.2f},{np.nanmax(v):.2f}] NA={int(np.isnan(v).sum())} "
          f"median_t={cl.time_months.median():.1f}mo " + ("OK" if not failed else "FAILED: " + ", ".join(failed)))
print("all checks passed" if bad == 0 else f"{bad} failed checks")
sys.exit(1 if bad else 0)
