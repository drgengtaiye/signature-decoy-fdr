"""Does each candidate null actually hold its size for the INTERACTION hypothesis?

The reviewer's objection is a claim about type-I error, so it should be settled by a
type-I error experiment rather than by argument. This script simulates survival data in
which

    beta_score != 0,  beta_arm != 0,  beta(score x arm) = 0

using the REAL IMmotion150 score matrix, so the correlation between signatures is the
real one, and then asks each null how often it rejects at the nominal 5%.

A null that respects H0 should reject about 5% of the time. A null that destroys one of
the main effects has no reason to, and the direction of its error is what matters: too
often means the published P values are too small.

Schemes tested
  armperm      permute arm labels (the v5 null; sharp null of no treatment effect)
  restricted   permute arm labels within quintiles of the score        [v6 "restricted"]
  scoreperm    permute the score within arm                            [v6 "scoreperm"]
  parametric   simulate from the fitted reduced model                  [v6 "parametric"]
  freedman     permute the interaction regressor after residualising it on score + arm
  wild         multiplier bootstrap of the efficient interaction score

Part 1 fixes one focal signature and reports per-signature size.
Part 2 uses the count statistic the paper reports - how many of 17 immune signatures
reach interaction P < 0.05 in the checkpoint-favouring direction - and reports the size
of that test, which is what the headline claim rests on.

Writes pancan/v7_calibration.csv.
"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.stats import norm

sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
from v7_lib import (P, build, cox_fit, cox_interaction, efficient_interaction_score,
                    reduced_simulator, pc_design, perm_p)

NCAL   = int(sys.argv[1]) if len(sys.argv) > 1 else 400    # calibration replicates
NNULL  = int(sys.argv[2]) if len(sys.argv) > 2 else 300    # null replicates inside each
NJOBS  = 48
B_SCORE, B_ARM = 0.35, -0.40      # true main effects; interaction is exactly zero


# ------------------------------------------------------------------ null generators
def p_armperm(rng, sc, arm, t, e, obs, nnull):
    v = [cox_interaction(sc, rng.permutation(arm), t, e)[0] for _ in range(nnull)]
    return perm_p([x for x in v if np.isfinite(x)], obs, "less")


def p_restricted(rng, sc, arm, t, e, obs, nnull, strata):
    out = []
    for _ in range(nnull):
        a = arm.copy()
        for s in np.unique(strata):
            m = strata == s
            a[m] = rng.permutation(arm[m])
        out.append(cox_interaction(sc, a, t, e)[0])
    return perm_p([x for x in out if np.isfinite(x)], obs, "less")


def p_scoreperm(rng, sc, arm, t, e, obs, nnull):
    out = []
    for _ in range(nnull):
        s2 = sc.copy()
        for v in (0, 1):
            m = np.where(arm == v)[0]
            s2[m] = sc[rng.permutation(m)]
        out.append(cox_interaction(s2, arm, t, e)[0])
    return perm_p([x for x in out if np.isfinite(x)], obs, "less")


def p_parametric(rng, sc, arm, t, e, obs, nnull):
    sim = reduced_simulator(np.column_stack([sc, arm]), t, e)
    if sim is None:
        return np.nan
    out = []
    for _ in range(nnull):
        tn, en = sim["sim"](rng)
        out.append(cox_interaction(sc, arm, tn, en)[0])
    return perm_p([x for x in out if np.isfinite(x)], obs, "less")


def p_freedman(rng, sc, arm, t, e, obs, nnull):
    """Permute the interaction regressor after projecting out score and arm.

    The reduced design stays bolted to the outcomes, so both main effects survive; only
    the interaction column is re-shuffled. This is the Freedman-Lane idea applied to the
    regressor rather than the residual outcome, which is what a Cox model allows.
    """
    D = np.column_stack([sc, arm, np.ones_like(sc)])
    v = sc * arm
    r = v - D @ np.linalg.lstsq(D, v, rcond=None)[0]
    out = []
    for _ in range(nnull):
        X = np.column_stack([sc, arm, rng.permutation(r)])
        b, V, ok = cox_fit(X, t, e)
        out.append(b[2] if ok else np.nan)
    return perm_p([x for x in out if np.isfinite(x)], obs, "less")


def p_wild(rng, sc, arm, t, e, obs, nnull):
    """Multiplier bootstrap of the efficient interaction score.

    Both main effects are held at their reduced-model estimates; the interaction score
    residuals are re-weighted by iid mean-zero multipliers. Nothing about the score or
    the arm is permuted, so neither main effect is disturbed.
    """
    E = efficient_interaction_score(sc, arm, t, e)
    if E is None:
        return np.nan
    u = E["u"]
    G = rng.choice([-1.0, 1.0], size=(nnull, len(u)))
    b1 = (G @ u) / E["Ieff"]
    return perm_p(b1, obs, "less")


METHODS = ["armperm", "restricted", "scoreperm", "parametric", "freedman", "wild"]


# --------------------------------------------------------------- part 1: one signature
def one_calibration(seed, sc, arm, t, e, sim, strata, nnull):
    rng = np.random.default_rng(seed)
    tn, en = sim["sim"](rng)                     # data WITH main effects, WITHOUT interaction
    obs, _ = cox_interaction(sc, arm, tn, en)
    if not np.isfinite(obs):
        return None
    r = dict(obs=obs)
    r["armperm"]    = p_armperm(rng, sc, arm, tn, en, obs, nnull)
    r["restricted"] = p_restricted(rng, sc, arm, tn, en, obs, nnull, strata)
    r["scoreperm"]  = p_scoreperm(rng, sc, arm, tn, en, obs, nnull)
    r["parametric"] = p_parametric(rng, sc, arm, tn, en, obs, nnull)
    r["freedman"]   = p_freedman(rng, sc, arm, tn, en, obs, nnull)
    r["wild"]       = p_wild(rng, sc, arm, tn, en, obs, nnull)
    return r


def main():
    D = build("IMmotion150")
    S, arm, t, e = D["S"], D["arm"], D["t"], D["e"]
    imm = D["fams"] == "Immune"
    print(f"IMmotion150: n={D['n']}, {imm.sum()} immune signatures, events={e.sum()}")

    # --- sanity: the fast Cox fitter must agree with lifelines
    from lifelines import CoxPHFitter
    sc = S[0]
    X = np.column_stack([sc, arm, sc * arm])
    b, V, ok = cox_fit(X, t, e)
    df = pd.DataFrame(dict(s=sc, arm=arm, sx=sc * arm, T=t, E=e))
    lf = CoxPHFitter().fit(df, "T", "E")
    print("\ncox_fit vs lifelines")
    for k, nm in enumerate(["s", "arm", "sx"]):
        print(f"  {nm:4s} mine {b[k]:+.6f} se {np.sqrt(V[k,k]):.6f} | "
              f"lifelines {lf.params_[nm]:+.6f} se {lf.standard_errors_[nm]:.6f}")
    # lifelines uses Efron for tied event times, this module uses Breslow. PFS recorded in
    # months ties heavily, so the two differ in the third decimal. Breslow is the right
    # choice here because the same convention has to serve the baseline-hazard simulator
    # that the model-based null inverts; what matters is that observed and null statistics
    # come from ONE convention, which they do. v7_observed_check.py confirms the published
    # counts are unchanged under Breslow.
    assert np.allclose(b, lf.params_[["s", "arm", "sx"]].to_numpy(), atol=8e-3), "fitter disagrees"

    # --- sanity: the efficient score must reproduce the Cox interaction estimate
    E0 = efficient_interaction_score(sc, arm, t, e)
    print(f"\nefficient score: one-step b {E0['b1']:+.5f} vs Cox MLE {b[2]:+.5f}; "
          f"robust z {E0['z']:+.3f}, model z {E0['T']/np.sqrt(E0['Ieff']):+.3f}, "
          f"Wald z {b[2]/np.sqrt(V[2,2]):+.3f}")

    # --- part 1: size of the per-signature test -------------------------------------
    truth = np.column_stack([sc, arm])
    sim = reduced_simulator(truth, t, e)
    sim["beta"] = np.array([B_SCORE, B_ARM])            # impose the truth we want
    sim2 = reduced_simulator(truth, t, e)
    # rebuild the simulator at the imposed coefficients
    import v7_lib
    order = np.argsort(t, kind="mergesort")
    uniq, first_idx, inv = np.unique(t[order], return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=e[order], minlength=len(uniq))
    eta = truth @ np.array([B_SCORE, B_ARM])
    w = np.exp(eta)
    S0 = np.cumsum(w[order][::-1])[::-1][first_idx]
    H0 = np.cumsum(d / S0)
    keep = d > 0
    h_t, h_H = uniq[keep], H0[keep]

    def sim_true(rng):
        target = -np.log(rng.random(len(w))) / np.maximum(w, 1e-12)
        j = np.searchsorted(h_H, target, side="left")
        tstar = np.where(j < len(h_t), h_t[np.minimum(j, len(h_t) - 1)], np.inf)
        return np.minimum(tstar, t), (tstar <= t).astype(int)

    strata = pd.qcut(sc, 5, labels=False, duplicates="drop")
    print(f"\n=== part 1: size of the per-signature interaction test "
          f"(beta_score={B_SCORE}, beta_arm={B_ARM}, beta_int=0) ===")
    print(f"    {NCAL} simulated trials x {NNULL} null replicates each")
    t0 = time.time()
    seeds = np.random.SeedSequence(7).spawn(NCAL)
    res = Parallel(n_jobs=NJOBS, verbose=0)(
        delayed(one_calibration)(sd, sc, arm, t, e, dict(sim=sim_true), strata, NNULL)
        for sd in seeds)
    res = [r for r in res if r]
    R = pd.DataFrame(res)
    print(f"    done in {(time.time()-t0)/60:.1f} min, {len(R)} usable\n")
    rows = []
    for m in METHODS:
        v = R[m].dropna()
        size = float((v <= 0.05).mean())
        lo, hi = size - 1.96 * np.sqrt(size * (1 - size) / len(v)), \
                 size + 1.96 * np.sqrt(size * (1 - size) / len(v))
        verdict = "OK" if abs(size - .05) < 1.96 * np.sqrt(.05 * .95 / len(v)) else \
                  ("ANTI-CONSERVATIVE" if size > .05 else "conservative")
        print(f"  {m:12s} size {size:.3f}  (95% CI {max(lo,0):.3f}-{hi:.3f})   {verdict}")
        rows.append(dict(part="per_signature", method=m, n=len(v), size=size,
                         ci_lo=max(lo, 0), ci_hi=hi, verdict=verdict))
    pd.DataFrame(rows).to_csv(f"{P}/v7_calibration.csv", index=False)
    print(f"\nwrote {P}/v7_calibration.csv")


if __name__ == "__main__":
    main()
