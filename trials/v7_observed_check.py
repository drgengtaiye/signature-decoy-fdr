"""Do the published observed statistics survive the switch to a Breslow fitter?

v4/v5/v6 computed the observed numbers with lifelines, which uses Efron for tied event
times. The v7 nulls need Breslow, because the model-based null simulates from a Breslow
baseline and the fitter that scores the simulated data has to match the one that
generated it. Before trusting any v7 P value, the observed side has to be shown to be the
same under both conventions - otherwise a shifted P could be a tie-handling artefact
rather than a change of null.

Also checks the efficient-score machinery that the multiplier bootstrap depends on:
  * the nuisance score residuals vanish at the reduced-model fit (they must, by definition)
  * the one-step interaction estimate tracks the Cox MLE
  * the multiplier null is centred on zero
"""
import sys, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
from v7_lib import (build, signatures, cox_fit, cox_interaction, cox_arm_slope,
                    efficient_interaction_score, score_residuals)
from lifelines import CoxPHFitter

SIGS = signatures()
allb, alle = [], []

for trial in ("IMmotion150", "PRINCE"):
    D = build(trial, SIGS)
    S, arm, t, e, fams, names = D["S"], D["arm"], D["t"], D["e"], D["fams"], D["names"]
    ev = t[e == 1]
    print(f"\n=== {trial}  n={D['n']}  events={e.sum()}  "
          f"distinct event times={len(np.unique(ev))}  "
          f"max tie multiplicity={int(pd.Series(ev).value_counts().max())}")

    hit_b = hit_e = 0
    diffs, ints_b, ints_e, b1s = [], [], [], []
    for i, fam in enumerate(fams):
        sc = S[i]
        bb, pb = cox_interaction(sc, arm, t, e)                       # Breslow (v7_lib)
        df = pd.DataFrame(dict(s=sc, arm=arm, sx=sc * arm, T=t, E=e))
        lf = CoxPHFitter().fit(df, "T", "E")                          # Efron (lifelines)
        be, pe = float(lf.params_["sx"]), float(lf.summary.loc["sx", "p"])
        if fam == "Immune":
            hit_b += int(np.isfinite(pb) and pb < .05 and bb < 0)
            hit_e += int(np.isfinite(pe) and pe < .05 and be < 0)
            li = cox_arm_slope(sc, arm, t, e, 1)
            lc = cox_arm_slope(sc, arm, t, e, 0)
            diffs.append(li - lc); ints_b.append(bb); ints_e.append(be)
            E = efficient_interaction_score(sc, arm, t, e)
            b1s.append(E["b1"] if E else np.nan)

    n_imm = int((fams == "Immune").sum())
    print(f"  immune hits  Breslow {hit_b}/{n_imm}   Efron(published) {hit_e}/{n_imm}")
    diffs, ints_b, ints_e, b1s = map(np.array, (diffs, ints_b, ints_e, b1s))
    print(f"  median(b_ICI - b_ctl)      {np.median(diffs):+.4f}")
    print(f"  median(b_int) Breslow      {np.median(ints_b):+.4f}   Efron {np.median(ints_e):+.4f}")
    print(f"  median(one-step b_int)     {np.nanmedian(b1s):+.4f}   "
          f"corr(one-step, MLE) {np.corrcoef(b1s[np.isfinite(b1s)], ints_b[np.isfinite(b1s)])[0,1]:.4f}")
    allb += list(ints_b); alle += list(ints_e)
    globals().setdefault("pool_diff", []).extend(list(diffs))
    globals().setdefault("pool_b1", []).extend(list(b1s))

print(f"\nPOOLED 34 tests: median(b_ICI-b_ctl) {np.median(pool_diff):+.4f} | "
      f"median(b_int) Breslow {np.median(allb):+.4f} Efron {np.median(alle):+.4f} | "
      f"median(one-step) {np.nanmedian(pool_b1):+.4f}")

# --- efficient-score sanity on one signature -----------------------------------------
D = build("IMmotion150", SIGS)
sc, arm, t, e = D["S"][0], D["arm"], D["t"], D["e"]
br, Vr, _ = cox_fit(np.column_stack([sc, arm]), t, e)
U = score_residuals(np.column_stack([sc, arm, sc * arm]), t, e, np.array([br[0], br[1], 0.0]))
E = efficient_interaction_score(sc, arm, t, e)
bfull, Vfull, _ = cox_fit(np.column_stack([sc, arm, sc * arm]), t, e)
rng = np.random.default_rng(0)
G = rng.choice([-1.0, 1.0], size=(20000, len(E["u"])))
null_b1 = (G @ E["u"]) / E["Ieff"]
print(f"\nefficient score checks (signature '{D['names'][0]}')")
print(f"  sum of nuisance score residuals  {np.abs(U[:, :2].sum(0)).max():.2e}  (should be ~0)")
print(f"  interaction score T = sum(u)     {E['T']:+.4f}")
print(f"  one-step b {E['b1']:+.5f}  vs Cox MLE {bfull[2]:+.5f}")
print(f"  robust se {np.sqrt((E['u']**2).sum())/E['Ieff']:.5f}  model se {E['se1']:.5f}  "
      f"Wald se {np.sqrt(Vfull[2,2]):.5f}")
print(f"  multiplier null: mean {null_b1.mean():+.5f}  sd {null_b1.std():.5f}  "
      f"(should be centred on 0)")
