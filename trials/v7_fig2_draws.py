"""Replicate-level draws from the valid (joint model-based) null, for Figure 2.

The summary CSVs keep only quantiles and P values; the figure needs the draws. This
re-runs the joint model-based bootstrap with the SAME SEEDS as
v7_interaction_null_valid.py, so every P value here reproduces that run exactly rather
than differing by Monte-Carlo noise, and writes out:

  v7_null_draws_hit_immune.csv   trial, K, replicate, hit_immune
  v7_null_draws_contrast.csv     trial, K, replicate, median_log_diff, wilcoxon_p
  v7_fig2_summary.csv            observed, null median and P for every statistic,
                                 per trial AND pooled, including the Wilcoxon calibration

Three statistics, all computed on the identical replicate:
  hit_immune       how many of the 17 immune signatures reach interaction P<0.05, beta<0
  median_log_diff  median over that trial's 17 immune signatures of (logHR_ICI - logHR_ctl)
  wilcoxon_p       paired signed-rank P over the same 17 pairs - kept so the legend's
                   claim that the Wilcoxon is badly calibrated can be restated under a
                   null that actually tests the interaction hypothesis
"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.stats import wilcoxon

sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
_a = sys.argv; sys.argv = [_a[0]]
from v7_interaction_null_valid import pass_stats          # identical statistic code path
sys.argv = _a
from v7_lib import (P, build, signatures, cox_interaction, cox_arm_slope,
                    reduced_simulator, pc_design, perm_p)

NREP  = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
NJOBS = 48
KS = [3, 4, 9, 15, 20]
TRIALS = ("IMmotion150", "PRINCE")


def wilx(d):
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    if d.size < 5 or np.allclose(d, 0):
        return np.nan
    try:
        return float(wilcoxon(d).pvalue)
    except Exception:
        return np.nan


def main():
    SIGS = signatures()
    D = {tr: build(tr, SIGS) for tr in TRIALS}
    obs = {}
    for tr in TRIALS:
        d = D[tr]
        imm = d["fams"] == "Immune"
        h, dif = pass_stats(d["S"], imm, d["arm"], d["t"], d["e"])
        obs[tr] = dict(hit=h, med=float(np.median(dif)), wp=wilx(dif), dif=dif)
        print(f"{tr}: observed hits {h}/{imm.sum()}  median log-HR diff {obs[tr]['med']:+.4f}  "
              f"paired Wilcoxon P {obs[tr]['wp']:.3g}")
    obs_pool_dif = np.concatenate([obs[tr]["dif"] for tr in TRIALS])
    obs_pool_med = float(np.median(obs_pool_dif))
    obs_pool_wp = wilx(obs_pool_dif)
    print(f"pooled: median {obs_pool_med:+.4f}  paired Wilcoxon P {obs_pool_wp:.3g}\n")

    hit_rows, con_rows, summ = [], [], []
    for K in KS:
        store = {}
        for tr in TRIALS:
            d = D[tr]
            S, arm, t, e = d["S"], d["arm"], d["t"], d["e"]
            imm = d["fams"] == "Immune"
            pcs, cum = pc_design(S, K)
            sim = reduced_simulator(np.column_stack([pcs, arm]), t, e)

            def one(seed):
                rng = np.random.default_rng(seed)
                tn, en = sim["sim"](rng)
                return pass_stats(S, imm, arm, tn, en)

            t0 = time.time()
            seeds = np.random.SeedSequence(90000 + K).spawn(NREP)   # SAME seeds as main run
            out = Parallel(n_jobs=NJOBS)(delayed(one)(s) for s in seeds)
            hits = np.array([o[0] for o in out], float)
            difs = np.vstack([o[1] for o in out])
            meds = np.median(difs, axis=1)
            wps = np.array([wilx(r) for r in difs])
            store[tr] = dict(hits=hits, difs=difs, meds=meds, wps=wps)

            for i in range(NREP):
                hit_rows.append(dict(trial=tr, K=K, replicate=i, hit_immune=int(hits[i])))
                con_rows.append(dict(trial=tr, K=K, replicate=i,
                                     median_log_diff=meds[i], wilcoxon_p=wps[i]))
            ev = int(e.sum()); npar = K + 1
            for stat, o, null, tail in [
                    ("hit_immune", obs[tr]["hit"], hits, "greater"),
                    ("median_log_diff", obs[tr]["med"], meds, "less"),
                    ("wilcoxon_p", obs[tr]["wp"], wps, "less")]:
                summ.append(dict(trial=tr, K=K, statistic=stat, observed=o,
                                 null_median=float(np.nanmedian(null)),
                                 null_p05=float(np.nanquantile(null, .05)),
                                 null_p95=float(np.nanquantile(null, .95)),
                                 n_rep=int(np.isfinite(null).sum()),
                                 perm_p=perm_p(null, o, tail),
                                 events_per_parameter=round(ev / npar, 1)))
            print(f"  {tr} K={K:2d} (var {cum[min(K,len(cum))-1]:.3f}, {ev/npar:.1f} ev/par): "
                  f"hits P={perm_p(hits,obs[tr]['hit'],'greater'):.4f} | "
                  f"contrast obs {obs[tr]['med']:+.4f} null med {np.median(meds):+.4f} "
                  f"P={perm_p(meds,obs[tr]['med'],'less'):.4f} | "
                  f"Wilcoxon obs {obs[tr]['wp']:.3g} null med {np.nanmedian(wps):.4g} "
                  f"P={perm_p(wps,obs[tr]['wp'],'less'):.4f}  [{time.time()-t0:.0f}s]")

        pooled_dif = np.hstack([store[tr]["difs"] for tr in TRIALS])
        pm = np.median(pooled_dif, axis=1)
        pw = np.array([wilx(r) for r in pooled_dif])
        for stat, o, null, tail in [("median_log_diff", obs_pool_med, pm, "less"),
                                    ("wilcoxon_p", obs_pool_wp, pw, "less")]:
            summ.append(dict(trial="pooled", K=K, statistic=stat, observed=o,
                             null_median=float(np.nanmedian(null)),
                             null_p05=float(np.nanquantile(null, .05)),
                             null_p95=float(np.nanquantile(null, .95)),
                             n_rep=int(np.isfinite(null).sum()),
                             perm_p=perm_p(null, o, tail), events_per_parameter=np.nan))
        print(f"  pooled K={K:2d}: contrast P={perm_p(pm,obs_pool_med,'less'):.4f} | "
              f"Wilcoxon obs {obs_pool_wp:.3g} null med {np.nanmedian(pw):.4g} "
              f"P={perm_p(pw,obs_pool_wp,'less'):.4f}\n")

    pd.DataFrame(hit_rows).to_csv(f"{P}/v7_null_draws_hit_immune.csv", index=False)
    pd.DataFrame(con_rows).to_csv(f"{P}/v7_null_draws_contrast.csv", index=False)
    pd.DataFrame(summ).to_csv(f"{P}/v7_fig2_summary.csv", index=False)
    print(f"wrote v7_null_draws_hit_immune.csv ({len(hit_rows)} rows), "
          f"v7_null_draws_contrast.csv, v7_fig2_summary.csv")


if __name__ == "__main__":
    main()
