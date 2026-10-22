"""Summaries of analyses 1, 2, 3 and 6 of validation/PLAN_round7_estimand_and_bias.md.

    python v12_truth_summary.py   -> icb/v12_truth/summary_response.csv, summary_tcga.csv, summary_null.csv

Per setting, ratios of totals over replicates: the realized FDR under the call-probability truth (sum V / sum R) and
its binomial variant, the plain and pi0-corrected estimates, the decomposition of the plain estimate's expectation
(null targets' decoy rates, true targets' decoy rates, the decoys' excess over null targets), and the coverage of the
target- and cluster-bootstrap intervals of both estimates, of (a) the setting's realized FDR and (b) the replicate's own
realized proportion. Response replicates also carry the expected-delta truth of round 6 and the gene-content truth.
"""
import glob, numpy as np, pandas as pd

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v12_truth"


def summarize(S, keys, extra=()):
    rows = []
    for k, x in S.groupby(keys, sort=False):
        R = x.R.sum(); y = x[x.R > 0]
        fdr = x.V_cp.sum() / R
        F0 = x.pi0 * x.E
        own = (y.V_cp / y.R)
        cov = lambda lo, hi, v: float(np.mean((y[lo] <= v + 1e-12) & (v <= y[hi] + 1e-12)))
        r = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
        r.update(reps=len(x), R=x.R.mean(), realized=100 * fdr, realized_bin=100 * x.V_bin.sum() / R,
                 estimate=100 * x.E.sum() / R, estimate_pi0=100 * F0.sum() / R, pi0=x.pi0.mean(),
                 n_true_cp=x.n_true_cp.mean(), exp_V=100 * x.exp_V.sum() / R,
                 part_null=100 * x.pi_null.sum() / R, part_true=100 * x.pi_true.sum() / R,
                 excess_null=100 * x.excess_null.sum() / R,
                 cov_t_fdr=cov("F_lo_t", "F_hi_t", fdr), cov_c_fdr=cov("F_lo_c", "F_hi_c", fdr),
                 cov0_t_fdr=cov("F0_lo_t", "F0_hi_t", fdr), cov0_c_fdr=cov("F0_lo_c", "F0_hi_c", fdr),
                 cov_t_own=float(np.mean((y.F_lo_t <= own + 1e-12) & (own <= y.F_hi_t + 1e-12))),
                 cov_c_own=float(np.mean((y.F_lo_c <= own + 1e-12) & (own <= y.F_hi_c + 1e-12))),
                 cov0_t_own=float(np.mean((y.F0_lo_t <= own + 1e-12) & (own <= y.F0_hi_t + 1e-12))),
                 cov0_c_own=float(np.mean((y.F0_lo_c <= own + 1e-12) & (own <= y.F0_hi_c + 1e-12))))
        for e in extra:
            r[e] = 100 * x[e].sum() / R
        rows.append(r)
    return pd.DataFrame(rows)


fs = sorted(glob.glob(f"{D}/signal_*.csv"))
if fs:
    S = pd.concat([pd.read_csv(f) for f in fs])
    O = summarize(S, ["scenario", "beta", "cohorts", "stat"], extra=("V_ed", "V_gc"))
    O.to_csv(f"{D}/summary_response.csv", index=False)
    pd.set_option("display.width", 250)
    print(O[["scenario", "beta", "cohorts", "stat", "reps", "R", "realized", "realized_bin", "V_ed", "V_gc", "estimate",
             "estimate_pi0", "pi0", "part_true", "excess_null", "cov_t_fdr", "cov_c_fdr", "cov0_c_fdr"]].round(2).to_string())
fs = sorted(glob.glob(f"{D}/tcga_planted_*.csv"))
if fs:
    T = pd.concat([pd.read_csv(f) for f in fs])
    O = summarize(T, ["kind", "beta"])
    O.to_csv(f"{D}/summary_tcga.csv", index=False)
    print(O[["kind", "beta", "reps", "R", "realized", "realized_bin", "estimate", "estimate_pi0", "pi0", "n_true_cp",
             "part_true", "excess_null", "cov_t_fdr", "cov_c_fdr", "cov0_t_fdr", "cov0_c_fdr"]].round(2).to_string())
P = pd.read_csv(f"{D}/null_pi0.csv")
rng = np.random.default_rng(713000)
rows = []
for st, x in P.groupby("stat"):
    E, E0, R = x.E.to_numpy(), x.E0.to_numpy(), x.R.to_numpy()
    bb = [rng.integers(0, len(x), len(x)) for _ in range(4000)]
    r0 = np.array([E0[i].mean() / R[i].mean() for i in bb]); r1 = np.array([E[i].mean() / R[i].mean() for i in bb])
    rows.append(dict(stat=st, perms=len(x), realized=R.mean(), plain=E.mean(), ratio=E.mean() / R.mean(),
                     ratio_lo=np.percentile(r1, 2.5), ratio_hi=np.percentile(r1, 97.5), pi0_mean=x.pi0.mean(),
                     corrected=E0.mean(), ratio_pi0=E0.mean() / R.mean(), ratio_pi0_lo=np.percentile(r0, 2.5),
                     ratio_pi0_hi=np.percentile(r0, 97.5)))
pd.DataFrame(rows).to_csv(f"{D}/summary_null.csv", index=False)
print(pd.DataFrame(rows).round(3).to_string())
