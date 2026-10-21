"""Dependence-aware inference for the response benchmark.

    python perm_inference.py            -> icb/v4_permutation_inference.csv

The headline "20 of 52 signatures reach P < 0.05, 2.6 expected, binomial P = 2.5e-13"
treats the 52 signatures as independent tests. They are not: several encode the same
programme, the model's whole-family holdout fails, and their scores are strongly
correlated within every cohort.

This replaces the binomial test with a permutation test that keeps that dependence intact.
Response labels are shuffled within cohort; the scores are untouched, so every correlation
among signatures — and between each signature and the matched random sets — is preserved
exactly. For each permutation the entire benchmark is recomputed: per-cohort AUC for each
signature and each random set, the null median and 95th percentile of the matching
size-and-engine rung, the across-cohort Wilcoxon signed-rank test per signature, and from
those the two headline statistics:

  count      how many signatures reach P < 0.05
  contrast   inflamed clearance % minus confounder clearance %

AUC is folded (max(a, 1-a)) exactly as the benchmark folds it, and because the scores are
fixed their ranks are fixed, so all AUCs for one permutation are a single matrix product.

Also reports the effective number of independent signatures (Li & Ji 2005) from the
eigenvalues of the signature-score correlation matrix in each cohort.
"""
import warnings, functools, glob, os, re; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.stats import rankdata, wilcoxon, norm

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
B = int(os.environ.get("NPERM", "2000"))
CHUNK = 250
EXCLUDE = {"ICB_VanDenEnde"}          # excluded from the headline analysis, as in v3
rng = np.random.default_rng(0)

fam = pd.read_csv(f"{D}/v3_dissociation.csv").set_index("signature").family.to_dict()


def folded_auc(ranks, Y, n1, n0):
    """AUC of every score column against every permuted label column, folded.
    ranks: n x S (fixed), Y: n x B of 0/1 labels."""
    a = (ranks.T @ Y - n1 * (n1 + 1) / 2) / (n1 * n0)
    return np.maximum(a, 1 - a)


@functools.lru_cache(maxsize=None)
def _signed_rank_cdf(n):
    """Exact null distribution of W+ for n pairs, by DP over subset sums of ranks 1..n."""
    cnt = np.zeros(n * (n + 1) // 2 + 1)
    cnt[0] = 1.0
    for r in range(1, n + 1):
        nxt = cnt.copy()
        nxt[r:] += cnt[:len(cnt) - r]
        cnt = nxt
    return np.cumsum(cnt) / 2.0 ** n


def signed_rank_p(Dg):
    """Two-sided exact Wilcoxon signed-rank P, vectorised over every axis but the first.

    Same statistic and same exact two-sided p as scipy.stats.wilcoxon for continuous data
    (no zeros, no ties), which is what the benchmark calls; validated against it below."""
    n = Dg.shape[0]
    R = rankdata(np.abs(Dg), axis=0)
    wp = np.where(Dg > 0, R, 0.0).sum(axis=0)
    T = np.minimum(wp, n * (n + 1) / 2.0 - wp)
    cdf = _signed_rank_cdf(n)
    return np.minimum(1.0, 2.0 * cdf[np.clip(np.rint(T).astype(int), 0, len(cdf) - 1)])


@functools.lru_cache(maxsize=None)
def _sign_matrix(n):
    """All 2^n sign assignments for n pairs, as a 2^n x n matrix of 0/1."""
    return ((np.arange(2 ** n)[:, None] >> np.arange(n)[None, :]) & 1).astype(float)


def signed_rank_p_tied(Dg):
    """Exact signed-rank P for columns with zeros or tied |differences|, vectorised.

    scipy (1.17) keeps the exact branch for n <= 50 even when zeros and ties are present:
    zeros are dropped, |differences| get average ranks, and the null is enumerated over
    sign flips of that observed rank vector — not over the ranks 1..n, which is why a
    fixed-rank distribution disagrees whenever two |differences| tie. Permuted AUCs are
    discrete, so this is the common case here rather than a rare one.

    P = 2 * P(W+ <= min(W+, W-)), capped at 1, which reproduces scipy on ties and zeros."""
    N, ncol = Dg.shape[0], Dg.shape[1] * Dg.shape[2]
    D2 = Dg.reshape(N, ncol)
    nz = D2 != 0
    R = np.where(nz, rankdata(np.where(nz, np.abs(D2), np.inf), axis=0), 0.0)
    out = np.ones(ncol)
    counts = nz.sum(axis=0)
    for n in np.unique(counts):
        if n == 0:
            continue
        sel = np.flatnonzero(counts == n)
        Rs = np.sort(R[:, sel], axis=0)[N - n:]          # the n non-zero ranks, ascending
        wp = np.where(D2[:, sel] > 0, R[:, sel], 0.0).sum(axis=0)
        T = np.minimum(wp, Rs.sum(axis=0) - wp)
        W = _sign_matrix(int(n)) @ Rs                    # 2^n x |sel| null values of W+
        out[sel] = np.minimum(1.0, 2.0 * (W <= T[None, :] + 1e-9).mean(axis=0))
    return out.reshape(Dg.shape[1], Dg.shape[2])


def li_ji(corr):
    lam = np.linalg.eigvalsh(corr)
    lam = np.clip(lam, 0, None)
    return float(np.sum((lam >= 1).astype(float) + (lam - np.floor(lam))))


cohorts = sorted(os.path.basename(p)[:-8] for p in glob.glob(f"{D}/scores/*_map.csv"))
print(f"cohorts with score matrices: {len(cohorts)}   permutations: {B}")

sig_union, per_cohort, meff_rows = [], [], []
for co in cohorts:
    cache = f"{D}/scores/{co}_ranks.npz"
    if os.path.exists(cache):       # ranks do not depend on the labels, so cache them
        z = np.load(cache, allow_pickle=True)
        y, ranks_s, ranks_n = z["y"], z["ranks_s"], z["ranks_n"]
        sig_names, rung_key, sig_rung = list(z["sig_names"]), z["rung_key"], z["sig_rung"]
        meff, mabs = float(z["meff"]), float(z["mean_abs_r"])
    else:
        S = pd.read_csv(f"{D}/scores/{co}_sig.csv.gz")
        N = pd.read_csv(f"{D}/scores/{co}_null.csv.gz")
        mp = pd.read_csv(f"{D}/scores/{co}_map.csv")
        y = S.response.to_numpy(int)
        sig_names = [c for c in S.columns if c not in ("sample", "response")]
        Smat = S[sig_names].to_numpy(float)
        null_names = [c for c in N.columns if c != "sample"]
        Nmat = N[null_names].to_numpy(float)
        rung_key = np.array([f"{m.group(1)}|{int(m.group(2))}" for m in
                             (re.match(r"([a-z]+)\|k(\d+)\|", c) for c in null_names)])
        sig_rung = pd.Series([f"{e}|{int(k)}" for e, k in zip(mp.engine, mp.null_size)],
                             index=mp.signature).reindex(sig_names).to_numpy()
        corr = np.corrcoef(Smat, rowvar=False)
        meff = li_ji(np.nan_to_num(corr, nan=0.0))
        mabs = float(np.nanmean(np.abs(corr[np.triu_indices_from(corr, 1)])))
        ranks_s, ranks_n = rankdata(Smat, axis=0), rankdata(Nmat, axis=0)
        np.savez_compressed(cache, y=y, ranks_s=ranks_s, ranks_n=ranks_n,
                            sig_names=np.array(sig_names, dtype=object),
                            rung_key=rung_key, sig_rung=sig_rung,
                            meff=meff, mean_abs_r=mabs)
    meff_rows.append(dict(cohort=co, n=len(y), n_signatures=len(sig_names),
                          meff_li_ji=meff, mean_abs_r=mabs))
    per_cohort.append(dict(cohort=co, y=y, ranks_s=ranks_s, ranks_n=ranks_n,
                           sig_names=sig_names, sig_rung=sig_rung, rung_key=rung_key))
    sig_union += [s for s in sig_names if s not in sig_union]
    print(f"  {co:20s} n={len(y):4d} sigs={len(sig_names):3d} "
          f"nullsets={ranks_n.shape[1]} Meff={meff:.1f}")

sig_union = sorted(sig_union)
S_IDX = {s: i for i, s in enumerate(sig_union)}
FAMV = np.array([fam.get(s, "Other") for s in sig_union])


def run(cohort_subset, label):
    """Observed statistics and their permutation null, over the given cohorts."""
    use = [c for c in per_cohort if c["cohort"] in cohort_subset]
    nC, nS = len(use), len(sig_union)
    # permutation 0 is the observed labelling
    delta = np.full((nC, nS, B + 1), np.nan)
    clear = np.zeros((nC, nS, B + 1), dtype=bool)
    valid = np.zeros((nC, nS), dtype=bool)
    for ci, c in enumerate(use):
        y, n = c["y"], len(c["y"])
        n1, n0 = int(y.sum()), int((1 - y).sum())
        Y = np.empty((n, B + 1))
        Y[:, 0] = y
        for b in range(1, B + 1):
            Y[:, b] = rng.permutation(y)
        auc_s = np.empty((len(c["sig_names"]), B + 1))
        auc_n = np.empty((c["ranks_n"].shape[1], B + 1))
        for lo in range(0, B + 1, CHUNK):
            hi = min(lo + CHUNK, B + 1)
            auc_s[:, lo:hi] = folded_auc(c["ranks_s"], Y[:, lo:hi], n1, n0)
            auc_n[:, lo:hi] = folded_auc(c["ranks_n"], Y[:, lo:hi], n1, n0)
        med, p95 = {}, {}
        for rk in np.unique(c["rung_key"]):
            m = c["rung_key"] == rk
            med[rk] = np.median(auc_n[m], axis=0)
            p95[rk] = np.quantile(auc_n[m], 0.95, axis=0)
        for si, s in enumerate(c["sig_names"]):
            j = S_IDX[s]; rk = c["sig_rung"][si]
            delta[ci, j] = auc_s[si] - med[rk]
            clear[ci, j] = auc_s[si] > p95[rk]
            valid[ci, j] = True
        print(f"    [{label}] {c['cohort']} done")

    # signatures are grouped by how many cohorts they appear in, so each group is one
    # vectorised exact signed-rank test over all signatures and all permutations at once
    counts = np.zeros(B + 1, int)
    by_n = {}
    for j in range(nS):
        k = int(valid[:, j].sum())
        if k >= 6:
            by_n.setdefault(k, []).append(j)
    checked = approx_checked = False
    nfall = ntest = 0
    for nj, js in sorted(by_n.items()):
        Dg = np.stack([delta[valid[:, j], j, :] for j in js], axis=1)   # nj x |js| x B+1
        Pg = signed_rank_p(Dg)
        # a permuted AUC can land exactly on its rung's null median, and |differences| can
        # tie; scipy leaves the exact distribution in both cases, so those columns take the
        # same route it would take rather than the vectorised exact one
        A = np.sort(np.abs(Dg), axis=0)
        bad = (Dg == 0).any(axis=0) | (np.diff(A, axis=0) == 0).any(axis=0)
        if bad.any():
            Pg = np.where(bad, signed_rank_p_tied(Dg), Pg)
            if not approx_checked:
                idx = np.argwhere(bad)[:100]
                ref = np.array([wilcoxon(Dg[:, a, b_]).pvalue for a, b_ in idx])
                assert np.allclose(ref, [Pg[a, b_] for a, b_ in idx], rtol=1e-8, atol=1e-10), \
                    "vectorised tie/zero branch != scipy"
                print(f"    [{label}] tie/zero branch validated against scipy on "
                      f"{len(idx)} columns")
                approx_checked = True
        nfall += int(bad.sum()); ntest += Pg.size
        counts += (Pg < 0.05).sum(axis=0)
        if not checked:          # validate the fast test against scipy at the real labels
            ref = np.array([wilcoxon(delta[valid[:, j], j, 0]).pvalue for j in js])
            assert np.allclose(ref, Pg[:, 0], atol=1e-10), "fast signed-rank != scipy"
            print(f"    [{label}] signed-rank validated against scipy for {len(js)} "
                  f"signatures (n = {nj} cohorts)")
            checked = True
    print(f"    [{label}] exact path for {100*(1-nfall/max(ntest,1)):.1f}% of "
          f"{ntest:,} tests; {nfall:,} fell back to scipy")
    inf_m = (FAMV == "Inflamed/IFN")[None, :] & valid
    cnf_m = (FAMV == "Confounder")[None, :] & valid
    contrast = np.array([100 * clear[:, :, b][inf_m].mean() - 100 * clear[:, :, b][cnf_m].mean()
                         for b in range(B + 1)])
    return counts, contrast, valid


out = []; draws = []
for label, subset in [("12 cohorts (headline)", set(c["cohort"] for c in per_cohort) - EXCLUDE),
                      ("13 cohorts (all)", set(c["cohort"] for c in per_cohort))]:
    counts, contrast, valid = run(subset, label)
    for name, v in [("signatures with P<0.05", counts), ("inflamed-confounder clearance (pp)", contrast)]:
        obs, null = v[0], v[1:]
        draws.append(pd.DataFrame({"analysis": label, "statistic": name, "perm": np.arange(1, len(null) + 1), "value": null}))
        p = (np.sum(null >= obs) + 1) / (len(null) + 1)
        out.append(dict(analysis=label, statistic=name, observed=obs,
                        perm_null_median=float(np.median(null)),
                        perm_null_mean=float(np.mean(null)),
                        perm_null_p95=float(np.quantile(null, 0.95)),
                        perm_null_max=float(np.max(null)), perm_P=p, n_perm=len(null),
                        n_cohorts=int(valid.shape[0]), n_signatures=int(valid.any(0).sum())))
        print(f"  {label:22s} {name:34s} observed {obs:7.2f}  null median "
              f"{np.median(null):6.2f}  95th {np.quantile(null,.95):6.2f}  P = {p:.4f}")

R = pd.DataFrame(out)
M = pd.DataFrame(meff_rows)
M["analysis"] = "per-cohort effective tests (Li & Ji)"
M["statistic"] = "Meff"
R = pd.concat([R, M], ignore_index=True)
R.to_csv(f"{D}/v4_permutation_inference.csv", index=False)
pd.concat(draws, ignore_index=True).to_csv(f"{D}/v4_permutation_null_draws.csv.gz", index=False)   # the null itself, for plotting
print(f"\nMeff across cohorts: median {M.meff_li_ji.median():.1f} of {len(sig_union)} "
      f"signatures (mean |r| = {M.mean_abs_r.mean():.2f})")
print("saved v4_permutation_inference.csv")
