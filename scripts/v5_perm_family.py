"""Permutation null for the standardised family contrast, preserving signature correlation.

    NPERM=2000 python v5_perm_family.py    -> icb/v5_perm_family.csv

The inflamed-versus-confounder contrast on the standardised shift is currently a
Mann-Whitney test over 20 and 8 signatures treated as independent observations. They are
not independent, so the test is repeated against a null built the same way as the one
behind Figure 1: response labels are shuffled within cohort while the score matrices stay
fixed, so every correlation between signatures is preserved exactly, and the whole
pipeline is recomputed on each permutation --- per-cohort AUC for each signature and each
matched random set, the null median and standard deviation of the signature's own rung,
the standardised response effect, the shift, and finally the family statistic.

Only the response side is exchangeable: the prognostic effects come from different
patients in TCGA, so z_prog is held at its observed value. The null therefore asks what
the family contrast looks like when the response side carries no signal but retains its
correlation structure --- which is the question the Mann-Whitney P was implicitly claiming
to answer.
"""
import warnings, functools, glob, os; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
B = int(os.environ.get("NPERM", "2000"))
CHUNK = 250
EXCLUDE = {"ICB_VanDenEnde"}
rng = np.random.default_rng(0)

EFF = pd.read_csv(f"{D}/v4_standardised_effects.csv")
fam = EFF.set_index("signature").family.to_dict()
zprog = EFF.set_index("signature").z_prog.to_dict()
sigs = sorted(EFF.signature)
S_IDX = {s: i for i, s in enumerate(sigs)}
nS = len(sigs)
FAMV = np.array([fam[s] for s in sigs])
ZPROG = np.array([zprog[s] for s in sigs])

cohorts = sorted(os.path.basename(p)[:-10] for p in glob.glob(f"{D}/scores/*_ranks.npz"))
cohorts = [c for c in cohorts if c not in EXCLUDE]
print(f"cohorts: {len(cohorts)}   permutations: {B}")


def folded_auc(ranks, Y, n1, n0):
    a = (ranks.T @ Y - n1 * (n1 + 1) / 2) / (n1 * n0)
    return np.maximum(a, 1 - a)


Zresp = np.full((len(cohorts), nS, B + 1), np.nan)
for ci, co in enumerate(cohorts):
    z = np.load(f"{D}/scores/{co}_ranks.npz", allow_pickle=True)
    y, ranks_s, ranks_n = z["y"], z["ranks_s"], z["ranks_n"]
    sig_names, rung_key, sig_rung = list(z["sig_names"]), z["rung_key"], z["sig_rung"]
    n = len(y); n1, n0 = int(y.sum()), int((1 - y).sum())
    Y = np.empty((n, B + 1)); Y[:, 0] = y
    for b in range(1, B + 1):
        Y[:, b] = rng.permutation(y)
    auc_s = np.empty((len(sig_names), B + 1))
    auc_n = np.empty((ranks_n.shape[1], B + 1))
    for lo in range(0, B + 1, CHUNK):
        hi = min(lo + CHUNK, B + 1)
        auc_s[:, lo:hi] = folded_auc(ranks_s, Y[:, lo:hi], n1, n0)
        auc_n[:, lo:hi] = folded_auc(ranks_n, Y[:, lo:hi], n1, n0)
    med, sd = {}, {}
    for rk in np.unique(rung_key):
        m = rung_key == rk
        med[rk] = np.median(auc_n[m], axis=0)
        sd[rk] = auc_n[m].std(axis=0, ddof=1)
    for si, s in enumerate(sig_names):
        if s not in S_IDX:
            continue
        rk = sig_rung[si]
        Zresp[ci, S_IDX[s]] = (auc_s[si] - med[rk]) / np.where(sd[rk] > 0, sd[rk], np.nan)
    print(f"  {co:20s} done")

z_resp = np.nanmedian(Zresp, axis=0)              # nS x (B+1), median across cohorts
z_shift = z_resp - ZPROG[:, None]
inf_m, cnf_m = FAMV == "Inflamed/IFN", FAMV == "Confounder"

stat_p, stat_rb, stat_med = [], [], []
for b in range(B + 1):
    a, c = z_shift[inf_m, b], z_shift[cnf_m, b]
    ok_a, ok_c = a[np.isfinite(a)], c[np.isfinite(c)]
    u = mannwhitneyu(ok_a, ok_c, alternative="two-sided")
    stat_p.append(u.pvalue)
    stat_rb.append(2 * u.statistic / (len(ok_a) * len(ok_c)) - 1)
    stat_med.append(np.median(ok_a) - np.median(ok_c))
stat_p, stat_rb, stat_med = map(np.array, (stat_p, stat_rb, stat_med))

rows = []
for name, v, tail in [("Mann-Whitney P (smaller is more extreme)", stat_p, "lower"),
                      ("rank-biserial r", stat_rb, "upper"),
                      ("median z-shift difference", stat_med, "upper")]:
    obs, null = v[0], v[1:]
    p = ((np.sum(null <= obs) + 1) / (len(null) + 1) if tail == "lower"
         else (np.sum(null >= obs) + 1) / (len(null) + 1))
    rows.append(dict(statistic=name, observed=obs, null_median=float(np.median(null)),
                     null_p05=float(np.quantile(null, 0.05)),
                     null_p95=float(np.quantile(null, 0.95)),
                     null_extreme=float(null.min() if tail == "lower" else null.max()),
                     perm_P=p, n_perm=len(null), n_cohorts=len(cohorts),
                     n_inflamed=int(inf_m.sum()), n_confounder=int(cnf_m.sum())))
    print(f"{name:42s} observed {obs:8.4g}  null median {np.median(null):8.4g}  "
          f"perm P = {p:.4f}")

pd.DataFrame(rows).to_csv(f"{D}/v5_perm_family.csv", index=False)
print("saved v5_perm_family.csv")
