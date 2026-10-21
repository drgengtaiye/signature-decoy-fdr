"""Cluster-aware inference for every statistic whose unit is a signature.

    python v5_clusters.py    -> icb/v5_signature_clusters.csv, icb/v5_cluster_stats.csv

Figure 1's count of significant signatures is already referred to a permutation null that
preserves the correlation between signatures. Every other signature-level statistic in the
paper — the inflamed-versus-confounder family contrast, the attribution correlations, the
triage held-out correlation — still treats the 52 signatures as independent observations.
They are not, so this script measures the dependence and then repeats those statistics in
two ways that respect it:

  one per cluster     draw one signature from each correlated cluster, recompute, repeat
  cluster bootstrap   resample whole clusters with replacement, not individual signatures

Dependence is measured on the scores themselves: the signature x signature correlation
within each cohort, averaged across cohorts on the Fisher z scale weighted by (n - 3).
Redundancy is |r| — two anti-correlated signatures are no more independent than two
correlated ones — so the clustering distance is 1 - |r|, average linkage.
"""
import warnings, functools, glob, os; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu, spearmanr

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
EXCLUDE = {"ICB_VanDenEnde"}
CUT = 0.50            # stated cut: merge while average |r| >= 0.50
NDRAW = 20000         # one-per-cluster draws
NBOOT = 10000         # cluster bootstrap resamples
rng = np.random.default_rng(0)

EFF = pd.read_csv(f"{D}/v4_standardised_effects.csv")
fam = EFF.set_index("signature").family.to_dict()

# ---------------------------------------------------------------- dependence structure
cohorts = sorted(os.path.basename(p)[:-11] for p in glob.glob(f"{D}/scores/*_sig.csv.gz"))
cohorts = [c for c in cohorts if c not in EXCLUDE]
sigs = sorted(EFF.signature)
S_IDX = {s: i for i, s in enumerate(sigs)}
nS = len(sigs)

zsum = np.zeros((nS, nS)); wsum = np.zeros((nS, nS))
for co in cohorts:
    S = pd.read_csv(f"{D}/scores/{co}_sig.csv.gz")
    cols = [c for c in S.columns if c in S_IDX]
    X = S[cols].to_numpy(float)
    r = np.corrcoef(X, rowvar=False)
    r = np.clip(np.nan_to_num(r, nan=0.0), -0.999999, 0.999999)
    z = np.arctanh(r)
    w = max(len(S) - 3, 1)
    idx = [S_IDX[c] for c in cols]
    zsum[np.ix_(idx, idx)] += z * w
    wsum[np.ix_(idx, idx)] += w
    print(f"  {co:20s} n={len(S):4d} signatures={len(cols)}")
R = np.tanh(np.divide(zsum, wsum, out=np.zeros_like(zsum), where=wsum > 0))
np.fill_diagonal(R, 1.0)

dist = 1.0 - np.abs(R)
np.fill_diagonal(dist, 0.0)
dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average")

cl_fixed = fcluster(Z, t=1.0 - CUT, criterion="distance")
# dynamic cut: the largest gap between consecutive merge heights, above a floor, so the
# tree is cut where it separates most cleanly rather than at a number chosen by hand
h = Z[:, 2]
elig = np.flatnonzero(h[1:] > 0.25) + 1
gaps = np.diff(h)
k = elig[np.argmax(gaps[elig - 1])]
t_dyn = (h[k - 1] + h[k]) / 2
cl_dyn = fcluster(Z, t=t_dyn, criterion="distance")

CL = pd.DataFrame(dict(signature=sigs, family=[fam[s] for s in sigs],
                       cluster_fixed=cl_fixed, cluster_dynamic=cl_dyn))
for name, lab in [("cluster_fixed", f"|r| >= {CUT}"), ("cluster_dynamic", f"dynamic cut at 1-|r| = {t_dyn:.2f}")]:
    sizes = CL[name].value_counts()
    print(f"\n{lab}: {CL[name].nunique()} clusters; sizes {sorted(sizes.tolist(), reverse=True)}")
    for c, g in CL.groupby(name):
        if len(g) > 1:
            print(f"   cluster {c:2d} ({len(g)}): " + ", ".join(g.signature))
mean_r_within = []
for s in sigs:
    c = CL.loc[CL.signature == s, "cluster_fixed"].iloc[0]
    mates = [S_IDX[m] for m in CL.loc[CL.cluster_fixed == c, "signature"] if m != s]
    mean_r_within.append(float(np.mean(np.abs(R[S_IDX[s], mates]))) if mates else np.nan)
CL["mean_abs_r_within_cluster"] = mean_r_within
CL.to_csv(f"{D}/v5_signature_clusters.csv", index=False)

# ------------------------------------------------------------------ statistics to repeat
A = pd.read_csv(f"{D}/singlecell_attribution.csv")
T = pd.read_csv(f"{D}/v3_sigtriage_final.csv")[["signature", "pred_loso", "d_resp"]]
LYM = ["z_Tcell", "z_NK", "z_Bcell", "z_Plasma"]
OFF = ["z_Fibro", "z_Endo", "z_Peri", "z_Epi"]
A = A.assign(lymphoid=A[LYM].max(1), offtarget=A[OFF].max(1))
d = (EFF.merge(A[["signature", "lymphoid", "offtarget"]], on="signature", how="left")
        .merge(T, on="signature", how="left")
        .merge(CL[["signature", "cluster_fixed", "cluster_dynamic"]], on="signature"))


def family_contrast(sub):
    a = sub.loc[sub.family == "Inflamed/IFN", "z_shift"].to_numpy()
    b = sub.loc[sub.family == "Confounder", "z_shift"].to_numpy()
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    # every statistic in this script returns (effect, P) in that order
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue, np.median(a) - np.median(b)


def rho(sub, x, y):
    s = sub[[x, y]].dropna()
    if len(s) < 6:
        return np.nan, np.nan
    r, p = spearmanr(s[x], s[y])
    return r, p


STATS = {
    "family contrast, inflamed vs confounder (z shift)": lambda s: family_contrast(s)[:2],  # (rank-biserial, P)
    "lymphoid attribution vs response (z)": lambda s: rho(s, "lymphoid", "z_resp"),
    "lymphoid attribution vs prognosis (z)": lambda s: rho(s, "lymphoid", "z_prog"),
    "off-target attribution vs shift (z)": lambda s: rho(s, "offtarget", "z_shift"),
    "triage held-out prediction vs measured effect": lambda s: rho(s, "pred_loso", "d_resp"),
}

rows = []
for key in ["cluster_fixed", "cluster_dynamic"]:
    groups = {c: g.index.to_numpy() for c, g in d.groupby(key)}
    clusters = list(groups)
    # (a) one signature per cluster
    draws = {name: [] for name in STATS}
    for _ in range(NDRAW):
        pick = [rng.choice(groups[c]) for c in clusters]
        sub = d.loc[pick]
        for name, fn in STATS.items():
            draws[name].append(fn(sub))
    # (b) cluster bootstrap: resample clusters with replacement, keep all their members
    boot = {name: [] for name in STATS}
    for _ in range(NBOOT):
        take = rng.choice(len(clusters), size=len(clusters), replace=True)
        sub = d.loc[np.concatenate([groups[clusters[i]] for i in take])]
        for name, fn in STATS.items():
            boot[name].append(fn(sub)[0])          # effect, whatever the statistic
    for name, fn in STATS.items():
        obs_e, obs_p = fn(d)
        arr = np.array(draws[name], float)
        effs, pvals = arr[:, 0], arr[:, 1]
        bt = np.array(boot[name], float)
        bt = bt[np.isfinite(bt)]
        rows.append(dict(
            clustering=key, n_clusters=len(clusters), statistic=name,
            observed_effect=obs_e, observed_P=obs_p,
            onepercluster_median_effect=np.nanmedian(effs),
            onepercluster_median_P=np.nanmedian(pvals),
            onepercluster_frac_P_lt_05=float(np.nanmean(pvals < 0.05)),
            onepercluster_P_p95=float(np.nanquantile(pvals, 0.95)),
            cluster_boot_lo=float(np.quantile(bt, 0.025)),
            cluster_boot_hi=float(np.quantile(bt, 0.975)),
            cluster_boot_P_two_sided=float(2 * min((bt <= 0).mean(), (bt >= 0).mean())),
            n_signatures=int(d[["z_shift"]].notna().sum().iloc[0])))
        print(f"[{key}] {name}\n    observed effect {obs_e:+.3f} (P {obs_p:.2g}) | "
              f"one-per-cluster median P {np.nanmedian(pvals):.3g}, {100*np.nanmean(pvals<0.05):.0f}% < 0.05 | "
              f"cluster bootstrap 95% CI {np.quantile(bt,0.025):+.3f} to {np.quantile(bt,0.975):+.3f}, "
              f"P {2*min((bt<=0).mean(),(bt>=0).mean()):.3g}")

pd.DataFrame(rows).to_csv(f"{D}/v5_cluster_stats.csv", index=False)
print("\nsaved v5_signature_clusters.csv and v5_cluster_stats.csv")
