"""Two checks on the within-arm comparison of v8_within_arm_null.py.

1. Power. The checkpoint arms are larger than the comparator arms (174 vs 89 in IMmotion150,
   65 vs 27 in PRINCE), so the contrast could be a difference in power rather than in
   information. The checkpoint arm is therefore subsampled to the comparator arm's size and
   its event count, with the null rebuilt inside each subsample, and the clearance recomputed.

2. Direction. Clearance was counted in the benefit direction only. This also counts the
   comparator arms in the opposite direction, to say whether the signatures carry no
   information there or information pointing the other way.

    python v8b_arm_power_check.py   -> pancan/v8b_arm_power_check.csv
"""
import warnings, functools, sys
warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L
from v8_within_arm_null import matrices, coef, LADDER, NRAND

P = L.P
NSUB = 40
rng = np.random.default_rng(81)
SIGS = L.signatures()
rows = []
for trial in ("IMmotion150", "PRINCE"):
    Z, arm, t, e = matrices(trial)
    pool = list(Z.index)
    sigs = {nm: [x for x in gl if x in Z.index] for nm, gl in SIGS.items()}
    sigs = {nm: g for nm, g in sigs.items() if len(g) >= 3 and (nm in L.IMM or nm.startswith("CTRL_"))}
    n_cmp = int((arm == 0).sum())
    idx_ck = np.where((arm == 1).values)[0]
    print(f"\n=== {trial}: checkpoint {len(idx_ck)} -> subsampled to {n_cmp}, {NSUB} draws")
    counts = []
    for b in range(NSUB):
        take = rng.choice(idx_ck, n_cmp, replace=False)
        ta, ea = t.values[take], e.values[take]
        null = {}
        for k in LADDER:
            v = [coef(Z.loc[rng.choice(pool, k, replace=False)].values[:, take].mean(0), ta, ea) for _ in range(NRAND)]
            null[k] = np.array([x for x in v if np.isfinite(x)])
        c = 0
        for nm, g in sigs.items():
            if nm not in L.IMM:
                continue
            k = min(LADDER, key=lambda x: abs(x - len(g)))
            c += int(coef(Z.loc[g].values[:, take].mean(0), ta, ea) < np.quantile(null[k], .05))
        counts.append(c)
        rows.append(dict(trial=trial, draw=b, n=n_cmp, events=int(ea.sum()), immune_clearing=c, of=17))
    print(f"  immune signatures clearing the null in the subsampled checkpoint arm: "
          f"median {int(np.median(counts))} of 17 (range {min(counts)}-{max(counts)}); full arm gives "
          f"{9 if trial == 'IMmotion150' else 6} of 17, the comparator arm 0 of 17")

# direction in the comparator arms, on the full data
W = pd.read_csv(f"{P}/v8_within_arm_null.csv")
print("\n=== comparator arms, both directions (immune signatures)")
for trial in ("IMmotion150", "PRINCE"):
    g = W[(W.trial == trial) & (W.arm == "comparator") & (W.family == "Immune")]
    worse = int((g.coef > g.null_med + 1.645 * g.null_sd).sum())
    print(f"  {trial}: benefit direction {int(g.clears.sum())} of {len(g)}; opposite direction "
          f"{worse} of {len(g)}; median z {g.z.median():+.2f}")
pd.DataFrame(rows).to_csv(f"{P}/v8b_arm_power_check.csv", index=False)
