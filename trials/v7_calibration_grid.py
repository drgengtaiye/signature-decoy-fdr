"""Where does each null actually break, as a function of the main effects it destroys?

v7_calibration.py showed that at beta_score=0.35, beta_arm=-0.40 only the
score-within-arm permutation is anti-conservative, while the two arm-label schemes held
their size. That is not a contradiction of the reviewer's objection - it is the
quantitative version of it. Breaking the treatment main effect can only distort the null
in proportion to how large that main effect is, and in these two trials the fitted arm
coefficient is almost zero (-0.003 in IMmotion150, +0.040 in PRINCE). Breaking the SCORE
main effect is a different matter, because the score main effect is not small.

So the honest statement is not "three nulls are invalid" but "each null is exact only
under a condition, and those conditions differ". This script maps the size of every
scheme over a grid of true main effects, with the interaction held at exactly zero, so
the manuscript can say which condition each null needs and whether the data satisfy it.

Writes pancan/v7_calibration_grid.csv
"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
_argv = sys.argv
sys.argv = [_argv[0]]                       # v7_calibration reads sys.argv at import
from v7_calibration import (p_armperm, p_restricted, p_scoreperm, p_parametric,
                            p_freedman, p_wild, METHODS)
sys.argv = _argv
from v7_lib import P, build, cox_interaction

NCAL  = int(sys.argv[1]) if len(sys.argv) > 1 else 400
NNULL = int(sys.argv[2]) if len(sys.argv) > 2 else 300
NJOBS = 48

B_SCORE = [0.0, 0.35, 0.70]
B_ARM   = [0.0, -0.20, -0.40, -0.80, -1.20]


def make_dgp(arm, t, e, sc, b_s, b_a):
    """Breslow-baseline generator with the given main effects and NO interaction."""
    order = np.argsort(t, kind="mergesort")
    uniq, first_idx, inv = np.unique(t[order], return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=e[order], minlength=len(uniq))
    w = np.exp(b_s * sc + b_a * arm)
    S0 = np.cumsum(w[order][::-1])[::-1][first_idx]
    H0 = np.cumsum(d / S0)
    keep = d > 0
    h_t, h_H = uniq[keep], H0[keep]

    def sim(rng):
        target = -np.log(rng.random(len(w))) / np.maximum(w, 1e-12)
        j = np.searchsorted(h_H, target, side="left")
        tstar = np.where(j < len(h_t), h_t[np.minimum(j, len(h_t) - 1)], np.inf)
        return np.minimum(tstar, t), (tstar <= t).astype(int)

    return sim


def one(seed, sc, arm, t, e, sim, strata, nnull):
    rng = np.random.default_rng(seed)
    tn, en = sim(rng)
    obs, _ = cox_interaction(sc, arm, tn, en)
    if not np.isfinite(obs):
        return None
    return dict(obs=obs,
                armperm=p_armperm(rng, sc, arm, tn, en, obs, nnull),
                restricted=p_restricted(rng, sc, arm, tn, en, obs, nnull, strata),
                scoreperm=p_scoreperm(rng, sc, arm, tn, en, obs, nnull),
                parametric=p_parametric(rng, sc, arm, tn, en, obs, nnull),
                freedman=p_freedman(rng, sc, arm, tn, en, obs, nnull),
                wild=p_wild(rng, sc, arm, tn, en, obs, nnull))


def main():
    D = build("IMmotion150")
    sc, arm, t, e = D["S"][0], D["arm"], D["t"], D["e"]
    strata = pd.qcut(sc, 5, labels=False, duplicates="drop")
    print(f"IMmotion150 n={D['n']} events={e.sum()}; focal signature '{D['names'][0]}'")
    print(f"{NCAL} simulated trials x {NNULL} null replicates per cell; "
          f"interaction is exactly zero everywhere\n")

    rows = []
    for b_s in B_SCORE:
        for b_a in B_ARM:
            sim = make_dgp(arm, t, e, sc, b_s, b_a)
            t0 = time.time()
            out = Parallel(n_jobs=NJOBS)(
                delayed(one)(s, sc, arm, t, e, sim, strata, NNULL)
                for s in np.random.SeedSequence(int(abs(b_s * 100) + abs(b_a * 1000))).spawn(NCAL))
            R = pd.DataFrame([r for r in out if r])
            line = f"  b_score={b_s:+.2f} b_arm={b_a:+.2f} |"
            for m in METHODS:
                v = R[m].dropna()
                size = float((v <= .05).mean())
                se = np.sqrt(max(size, 1e-9) * (1 - size) / max(len(v), 1))
                rows.append(dict(b_score=b_s, b_arm=b_a, method=m, n=len(v), size=size,
                                 se=se, ci_lo=max(size - 1.96 * se, 0), ci_hi=size + 1.96 * se,
                                 inflated=bool(size - 1.96 * se > 0.05)))
                line += f" {m[:4]} {size:.3f}{'*' if size - 1.96*se > .05 else ' '}|"
            print(line + f"  [{time.time()-t0:.0f}s]")

    T = pd.DataFrame(rows)
    T.to_csv(f"{P}/v7_calibration_grid.csv", index=False)
    print("\n=== size by method (rows: true main effects; * = significantly above 5%) ===")
    print(T.pivot_table(index=["b_score", "b_arm"], columns="method",
                        values="size").round(3).to_string())
    print("\ncells where each scheme is significantly anti-conservative:")
    for m in METHODS:
        bad = T[(T.method == m) & T.inflated]
        if len(bad):
            print(f"  {m:12s} " + ", ".join(
                f"(bs={r.b_score:+.2f},ba={r.b_arm:+.2f}) {r.size:.3f}" for r in bad.itertuples()))
        else:
            print(f"  {m:12s} none")
    print(f"\nwrote {P}/v7_calibration_grid.csv")


if __name__ == "__main__":
    main()
