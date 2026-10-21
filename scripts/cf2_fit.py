"""The closed-form matched null, fitted and validated leave-one-cohort-out.

    python cf2_fit.py   -> icb/cf_null_predictions.csv, icb/cf_null_model.json

Pre-specified in validation/PLAN_closed_form_null.md (SHA-256 edd47e92343ea3bc…, 2026-09-23 21:08 UTC).
"""
import json, functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import norm

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
S = pd.read_csv(f"{D}/cf_cohort_stats.csv").set_index("cohort")
N = pd.read_csv(f"{D}/native_response_matched_v2.csv")
R = N.drop_duplicates(["cohort", "engine", "null_size"])[
    ["cohort", "engine", "null_size", "rand_med", "rand_p95", "rand_mean", "rand_sd"]].copy()
R = R[R.cohort.isin(S.index)].reset_index(drop=True)
R["sigma1"] = S.sigma1.reindex(R.cohort).to_numpy()
R["rho"] = S.rho.reindex(R.cohort).to_numpy()
ENG = sorted(R.engine.unique())
print(f"{len(R)} rungs | {R.cohort.nunique()} cohorts | engines {ENG}")

# the unfolded AUC spread implied by the drawn folded null: for a half-normal at 0.5,
# sd(folded) = sigma * sqrt(1 - 2/pi), so sigma = rand_sd / sqrt(1 - 2/pi)
K = np.sqrt(1 - 2 / np.pi)
R["sigma_obs"] = R.rand_sd / K


def predict(p, d):
    """sigma(k) = c_e * (f * sigma1 + (sigma1 - f * sigma1) / sqrt(k_eff))"""
    f = p[0]
    c = {e: p[1 + i] for i, e in enumerate(ENG)}
    k = d.null_size.to_numpy(float)
    keff = k / (1 + (k - 1) * d.rho.to_numpy(float))
    s1 = d.sigma1.to_numpy(float)
    floor = f * s1
    ce = np.array([c[e] for e in d.engine])
    return ce * (floor + (s1 - floor) / np.sqrt(keff))


def fit(d):
    p0 = np.r_[0.3, np.ones(len(ENG))]
    r = least_squares(lambda p: np.log(predict(p, d)) - np.log(d.sigma_obs.to_numpy(float)), p0,
                      bounds=(np.r_[0.0, np.zeros(len(ENG))], np.r_[1.0, np.full(len(ENG), 5.0)]))
    return r.x


rows = []
for c in sorted(R.cohort.unique()):                      # leave one cohort out
    p = fit(R[R.cohort != c])
    d = R[R.cohort == c].copy()
    sig = predict(p, d)
    d["sd_hat"] = sig * K
    d["med_hat"] = 0.5 + sig * norm.ppf(0.75)            # folded median
    d["p95_hat"] = 0.5 + sig * norm.ppf(0.975)           # folded 95th percentile
    rows.append(d)
P = pd.concat(rows)
P["err_sd"] = (P.sd_hat - P.rand_sd) / P.rand_sd
P["err_p95"] = (P.p95_hat - P.rand_p95) / P.rand_p95
P["err_med"] = (P.med_hat - P.rand_med) / P.rand_med
print(f"\nleave-one-cohort-out relative error, median absolute:")
for col, lab in [("err_sd", "null s.d."), ("err_med", "null median"), ("err_p95", "null 95th percentile")]:
    print(f"  {lab:22s} {np.median(np.abs(P[col])):.3f}   (90th percentile {np.quantile(np.abs(P[col]), .9):.3f})")

# what it changes for the 52 signatures: does a signature still clear its null at the 95th percentile?
M = N.merge(P[["cohort", "engine", "null_size", "p95_hat", "med_hat", "sd_hat"]],
            on=["cohort", "engine", "null_size"], how="inner")
flip = (M.AUC > M.rand_p95) != (M.AUC > M.p95_hat)
print(f"\nsignature x cohort tests: {len(M)}, calls that change with the predicted null: {flip.sum()} ({100*flip.mean():.1f}%)")
ok = np.median(np.abs(P.err_sd)) < 0.10 and flip.mean() < 0.05
print(f"pre-specified criterion (median |error| < 10% and < 5% of calls change): {'MET' if ok else 'NOT MET'}")
p_all = fit(R)
json.dump(dict(engines=ENG, floor_fraction=float(p_all[0]),
               engine_constant={e: float(p_all[1+i]) for i, e in enumerate(ENG)},
               note="sigma(k) = c_e * (f*sigma1 + (1-f)*sigma1/sqrt(k_eff)); k_eff = k/(1+(k-1)*rho)"),
          open(f"{D}/cf_null_model.json", "w"), indent=1)
P.to_csv(f"{D}/cf_null_predictions.csv", index=False)
print("\nfitted on all cohorts:", json.load(open(f"{D}/cf_null_model.json")))
