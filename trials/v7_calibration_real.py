"""Type-I error in the regime these two trials actually occupy - and on PRINCE's design.

v7_calibration_grid.py swept beta_score up to 0.70 and beta_arm down to -1.20. The fitted
reduced models say the real trials are nowhere near that: per signature,

    IMmotion150   median beta_score -0.104,  median beta_arm -0.035 (max |0.058|)
    PRINCE        median beta_score -0.143,  median beta_arm +0.065 (max |0.154|)

Every non-zero cell of that grid is therefore larger than anything these data show, which
makes it the wrong evidence for the claim "the nulls used in this paper are valid". This
script calibrates at the observed values instead, and - more importantly - runs the whole
exercise on PRINCE's own design (n=92, 83 events), where the per-signature exact
bootstrap already showed the Wald test is anti-conservative by a factor of 1.58.

More calibration replicates than the grid, because the question here is whether a size is
0.05 or 0.07, which 400 replicates cannot resolve.

Writes pancan/v7_calibration_real.csv
"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
_a = sys.argv; sys.argv = [_a[0]]
from v7_calibration import (p_armperm, p_restricted, p_scoreperm, p_parametric,
                            p_freedman, p_wild, METHODS)
sys.argv = _a
from v7_lib import P, build, signatures, cox_interaction
from v7_calibration_grid import make_dgp

NCAL  = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
NNULL = int(sys.argv[2]) if len(sys.argv) > 2 else 400
NJOBS = 40

# (trial design, beta_score, beta_arm, label)
CELLS = [
    ("IMmotion150", -0.104, -0.035, "IMmotion150 observed median"),
    ("IMmotion150", -0.300, -0.058, "IMmotion150, strongest arm effect seen"),
    ("PRINCE",      -0.143, +0.065, "PRINCE observed median"),
    ("PRINCE",      -0.300, +0.154, "PRINCE, strongest arm effect seen"),
]


def one(seed, sc, arm, t, e, sim, strata, nnull):
    rng = np.random.default_rng(seed)
    tn, en = sim(rng)
    obs, _ = cox_interaction(sc, arm, tn, en)
    if not np.isfinite(obs):
        return None
    return dict(armperm=p_armperm(rng, sc, arm, tn, en, obs, nnull),
                restricted=p_restricted(rng, sc, arm, tn, en, obs, nnull, strata),
                scoreperm=p_scoreperm(rng, sc, arm, tn, en, obs, nnull),
                parametric=p_parametric(rng, sc, arm, tn, en, obs, nnull),
                freedman=p_freedman(rng, sc, arm, tn, en, obs, nnull),
                wild=p_wild(rng, sc, arm, tn, en, obs, nnull))


def main():
    SIGS = signatures()
    D = {tr: build(tr, SIGS) for tr in ("IMmotion150", "PRINCE")}
    print(f"{NCAL} simulated trials x {NNULL} null replicates per cell; "
          f"interaction exactly zero everywhere\n")
    rows = []
    for trial, b_s, b_a, label in CELLS:
        d = D[trial]
        sc, arm, t, e = d["S"][0], d["arm"], d["t"], d["e"]
        strata = pd.qcut(sc, 5, labels=False, duplicates="drop")
        sim = make_dgp(arm, t, e, sc, b_s, b_a)
        t0 = time.time()
        out = Parallel(n_jobs=NJOBS)(
            delayed(one)(s, sc, arm, t, e, sim, strata, NNULL)
            for s in np.random.SeedSequence(abs(hash(label)) % 2**31).spawn(NCAL))
        R = pd.DataFrame([r for r in out if r])
        line = f"  {label:42s} n={d['n']:3d} ev={e.sum():3d} |"
        for m in METHODS:
            v = R[m].dropna()
            size = float((v <= .05).mean())
            se = np.sqrt(max(size, 1e-9) * (1 - size) / max(len(v), 1))
            rows.append(dict(trial=trial, label=label, b_score=b_s, b_arm=b_a, method=m,
                             n=len(v), size=size, se=se, ci_lo=max(size - 1.96 * se, 0),
                             ci_hi=size + 1.96 * se, inflated=bool(size - 1.96 * se > .05)))
            line += f" {m[:4]} {size:.3f}{'*' if size - 1.96*se > .05 else ' '}|"
        print(line + f"  [{time.time()-t0:.0f}s]")
    T = pd.DataFrame(rows)
    T.to_csv(f"{P}/v7_calibration_real.csv", index=False)
    print("\n=== size in the regime the data actually occupy (* = above 5%) ===")
    print(T.pivot_table(index="label", columns="method", values="size").round(3).to_string())
    print("\nper method, across the four realistic cells:")
    for m in METHODS:
        s = T[T.method == m]
        print(f"  {m:12s} sizes " + " ".join(f"{x:.3f}" for x in s["size"]) +
              f"   inflated {int(s.inflated.sum())}/{len(s)}")
    print(f"\nwrote {P}/v7_calibration_real.csv")


if __name__ == "__main__":
    main()
