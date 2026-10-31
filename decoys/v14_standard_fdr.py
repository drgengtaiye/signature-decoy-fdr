"""Analyses 1-3 of validation/PLAN_round9_standard_fdr_comparison.md: the decoy estimate beside the nominal
(Benjamini-Hochberg / Storey) and permutation (SAM-type) estimates of the false discovery rate of the same calls.

    python v14_standard_fdr.py   -> icb/v14_standard/standard_real.csv, standard_globalnull.csv, standard_planted.csv,
                                    standard_planted_replicates.csv

Inputs are the stored outputs of rounds 7 (package runs and planted signal); nothing is recomputed from expression.
"""
import glob, os, numpy as np, pandas as pd

ICB = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
PKG, TRU, OUT = f"{ICB}/v12_pkg", f"{ICB}/v12_truth", f"{ICB}/v14_standard"
os.makedirs(OUT, exist_ok=True)
T = 0.025  # null rate of the call: two-sided P < 0.05 with a positive median

# analysis, statistic, calls file, decoy FDR file, global-null per-permutation file, decoy global-null summary
ANALYSES = [
    ("response", "fold", "response", f"{TRU}/null_pi0.csv", f"{TRU}/summary_null.csv"),
    ("response", "dir", "response", f"{TRU}/null_pi0.csv", f"{TRU}/summary_null.csv"),
    ("tcga_OS", "fold", "tcga_OS", f"{PKG}/tcga_OS_globalnull_per_perm.csv", f"{PKG}/tcga_OS_globalnull_summary.csv"),
    ("breast_pc1", "fold", "breast_pc1", f"{PKG}/breast_pc1_globalnull_per_perm.csv",
     f"{PKG}/breast_pc1_globalnull_summary.csv"),
    ("hallmark_OS", "fold", "hallmark_OS", f"{PKG}/hallmark_OS_globalnull_per_perm.csv",
     f"{PKG}/hallmark_OS_globalnull_summary.csv"),
]


def one_sided(c):
    """One-sided P of the call (positive direction) from the two-sided signed-rank P and the sign of the median."""
    return np.where(c.med_delta > 0, c.p / 2, 1 - c.p / 2)


def bh_adjust(p):
    m = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1)
    return out


def perm_R(path, stat):
    P = pd.read_csv(path)
    return P.loc[P.stat == stat, "R"].to_numpy(float)


rng = np.random.default_rng(914000)
real, gnull = [], []
for name, stat, stem, gpath, gsum in ANALYSES:
    C = pd.read_csv(f"{PKG}/{stem}_calls.csv.gz")
    C = C[(C.kind == "target") & (C.stat == stat)]
    p1 = one_sided(C)
    assert ((p1 < T) == C.call.to_numpy()).all(), name
    m, R = len(C), int((p1 < T).sum())
    pi0_st = min(1.0, (p1 > 0.5).sum() / (0.5 * m))
    q = bh_adjust(p1)
    Rp = perm_R(gpath, stat)
    F = pd.read_csv(f"{PKG}/{stem}_fdr.csv")
    F = F[(F.stat == stat) & (F["class"] == "coex") & (F.boot == "clusters")].iloc[0]
    real.append(dict(analysis=name, stat=stat, m=m, R=R, nominal=100 * m * T / R, pi0_storey=pi0_st,
                     storey=100 * pi0_st * m * T / R, perm_Rbar=Rp.mean(), perm_perms=len(Rp),
                     permutation=100 * Rp.mean() / R, permutation_storey=100 * pi0_st * Rp.mean() / R,
                     bh_05=int((q <= 0.05).sum()), bh_10=int((q <= 0.10).sum()),
                     bh_level_R=100 * m * np.sort(p1)[R - 1] / R,
                     decoy=F.fdr, decoy_lo=F.fdr_lo, decoy_hi=F.fdr_hi, decoy_pi0=F.pi0, decoy_pi0corr=F.fdr_pi0))
    bb = np.array([Rp[rng.integers(0, len(Rp), len(Rp))].mean() for _ in range(4000)])
    G = pd.read_csv(gsum); G = G[G.stat == stat].iloc[0]
    gnull.append(dict(analysis=name, stat=stat, perms=len(Rp), realized=Rp.mean(), nominal_expect=m * T,
                      nominal_ratio=m * T / Rp.mean(), nominal_lo=np.percentile(m * T / bb, 2.5),
                      nominal_hi=np.percentile(m * T / bb, 97.5), decoy_ratio=G.ratio, decoy_lo=G.ratio_lo,
                      decoy_hi=G.ratio_hi, decoy_pi0_ratio=G.get("ratio_pi0", np.nan),
                      decoy_pi0_lo=G.get("ratio_pi0_lo", np.nan), decoy_pi0_hi=G.get("ratio_pi0_hi", np.nan)))
real, gnull = pd.DataFrame(real), pd.DataFrame(gnull)
real.to_csv(f"{OUT}/standard_real.csv", index=False)
gnull.to_csv(f"{OUT}/standard_globalnull.csv", index=False)

# planted signal: ratios of totals over replicates, as v12_truth_summary.py
Rbar = {("response", s): perm_R(f"{TRU}/null_pi0.csv", s).mean() for s in ("fold", "dir")}
Rbar[("tcga", "fold")] = perm_R(f"{PKG}/tcga_OS_globalnull_per_perm.csv", "fold").mean()
S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{TRU}/signal_*.csv"))])
S["setting"] = S.scenario + ", beta " + S.beta.astype(str) + ", " + S.cohorts.astype(str) + " cohorts"
S["data"] = "response"
U = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{TRU}/tcga_planted_*.csv"))])
U["setting"] = "TCGA " + U.kind + ", beta " + U.beta.astype(str)
U["data"], U["stat"], U["cohorts"], U["scenario"] = "tcga", "fold", 20, "TCGA " + U.kind
U["V_gc"] = np.nan
rows, reps = [], []
for (data, setting, stat), x in pd.concat([S, U]).groupby(["data", "setting", "stat"], sort=False):
    m, n, sR = 52, len(x), x.R.sum()
    rb = Rbar[(data, stat)] if (data == "tcga" or x.cohorts.iloc[0] == 12) else np.nan
    rows.append(dict(data=data, setting=setting, scenario=x.scenario.iloc[0], beta=x.beta.iloc[0],
                     cohorts=x.cohorts.iloc[0], stat=stat, reps=n, R=x.R.mean(),
                     realized=100 * x.V_cp.sum() / sR, realized_gc=100 * x.V_gc.sum(min_count=1) / sR,
                     nominal=100 * n * m * T / sR, permutation=100 * n * rb / sR,
                     decoy=100 * x.E.sum() / sR, decoy_pi0=100 * (x.pi0 * x.E).sum() / sR))
    y = x[x.R > 0]
    reps.append(pd.DataFrame(dict(data=data, setting=setting, scenario=x.scenario.iloc[0], beta=x.beta.iloc[0],
                                  cohorts=x.cohorts.iloc[0], stat=stat, R=y.R, realized=100 * y.V_cp / y.R,
                                  nominal=100 * m * T / y.R, permutation=100 * rb / y.R, decoy=100 * y.E / y.R,
                                  decoy_pi0=100 * y.pi0 * y.E / y.R)))
P = pd.DataFrame(rows)
P.to_csv(f"{OUT}/standard_planted.csv", index=False)
Q = pd.concat(reps)
Q.to_csv(f"{OUT}/standard_planted_replicates.csv", index=False)

pd.set_option("display.width", 250)
print(real.round(2).to_string()); print(); print(gnull.round(3).to_string()); print()
print(P.drop(columns=["data"]).round(1).to_string()); print()
for e in ("nominal", "permutation", "decoy", "decoy_pi0"):
    a = P[P[e].notna()]
    g = a[a.realized_gc.notna()]
    print(f"{e:12s} below realized (call-probability truth): {(a[e] < a.realized).sum()} of {len(a)};"
          f" below realized (gene-content truth, response): {(g[e] < g.realized_gc).sum()} of {len(g)}")
z = Q[(Q.data == "response") & (Q.cohorts == 12) & (Q.beta == 1.0) & (Q.stat == "fold")]
print(f"\nreplicates with calls, 12 cohorts, folded, beta 1: {len(z)}")
for e in ("nominal", "permutation", "decoy", "decoy_pi0"):
    print(f"  {e:12s} at or above the realized proportion: {100 * np.mean(z[e] >= z.realized - 1e-9):.1f}%")
