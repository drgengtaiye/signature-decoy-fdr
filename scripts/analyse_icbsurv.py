"""Survival under checkpoint blockade: the pre-specified analysis of null_matched_icbsurv.R.

    python analyse_icbsurv.py    -> icb/icbsurv_summary.csv, printed summary

Primary: median gain in C-index per family for overall survival, the share of tests clearing
the null's 95th percentile against the 5% expected, and the inflamed-versus-confounder contrast
on the z scale, tested by a one-sided bootstrap over the 13 signature correlation clusters.
The comparison of interest is with the same signatures' survival effect without checkpoint
blockade (TCGA) and their response effect, both already measured with the same instrument.
"""
import numpy as np, pandas as pd

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
S = pd.read_csv(f"{D}/native_icbsurv_matched.csv")
X = pd.read_csv(f"{D}/v3_dissociation.csv").set_index("signature")
CL = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature").cluster_fixed
S["family"] = S.signature.map(X.family)
S["d"] = S.C - S.rand_med
S["z"] = (S.C - S.rand_med) / S.rand_sd
S["clear"] = S.C > S.rand_p95
rng = np.random.default_rng(13)


def per_signature(ep):
    g = S[S.endpoint == ep]
    return g.groupby("signature").agg(d=("d", "median"), z=("z", "median"), clear=("clear", "mean"),
                                      n_cohorts=("cohort", "nunique"))


def cluster_boot(v, fam, n=20000):
    """One-sided bootstrap over correlation clusters: inflamed shift above confounder."""
    cl = v.index.map(CL)
    groups = {c: v[cl == c] for c in pd.unique(cl.dropna())}
    keys = list(groups)
    obs = v[fam == "Inflamed/IFN"].median() - v[fam == "Confounder"].median()
    draws = []
    for _ in range(n):
        pick = rng.choice(keys, len(keys), replace=True)
        vv = pd.concat([groups[k] for k in pick]); ff = fam.reindex(vv.index)
        a, b = vv[ff.values == "Inflamed/IFN"], vv[ff.values == "Confounder"]
        if len(a) and len(b):
            draws.append(a.median() - b.median())
    draws = np.array(draws)
    return obs, float(np.quantile(draws, .025)), float(np.quantile(draws, .975)), float((draws <= 0).mean())


print(f"cohorts: {S.cohort.nunique()} | tests: {len(S)} | endpoints: {sorted(S.endpoint.unique())}")
for ep in ["OS", "PFS"]:
    g = S[S.endpoint == ep]
    if not len(g):
        continue
    ps = per_signature(ep); fam = ps.index.map(X.family).to_series(index=ps.index)
    print(f"\n=== {ep}: {g.cohort.nunique()} cohorts, {g.n.groupby(g.cohort).max().sum()} patients, "
          f"{g.events.groupby(g.cohort).max().sum()} events, {len(g)} tests")
    print("  median gain in C-index by family:", ps.groupby(fam).d.median().round(4).to_dict())
    print("  median z by family:              ", ps.groupby(fam).z.median().round(2).to_dict())
    print("  tests clearing the null's 95th percentile:",
          {f: round(100 * g[g.family == f].clear.mean(), 1) for f in ["Inflamed/IFN", "Confounder", "Other"]}, "% (5% expected)")
    obs, lo, hi, p = cluster_boot(ps.z, fam)
    print(f"  inflamed − confounder (z): {obs:+.2f} (95% CI {lo:+.2f} to {hi:+.2f}), one-sided cluster bootstrap P = {p:.3f}")
    if ep == "OS":
        ps.assign(family=fam).to_csv(f"{D}/icbsurv_summary.csv")
print("\nFor comparison, same instrument elsewhere (null s.d.):")
V = pd.read_csv(f"{D}/v4_standardised_effects.csv")
print("  response, inflamed  %+.2f | survival without checkpoint blockade, inflamed %+.2f" %
      (V[V.family == "Inflamed/IFN"].z_resp.median(), V[V.family == "Inflamed/IFN"].z_prog.median()))
print("  response, confounder %+.2f | survival without checkpoint blockade, confounder %+.2f" %
      (V[V.family == "Confounder"].z_resp.median(), V[V.family == "Confounder"].z_prog.median()))
