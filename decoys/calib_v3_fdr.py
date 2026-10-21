"""Target-decoy false discovery rate for the signature library (calibration v3).

    python calib_v3_fdr.py  -> calib_v3_units.csv, calib_v3_fdr.csv, calib_v3_per_signature.csv

Each unit (one of the 52 signatures, or one of its decoys) is tested exactly as the paper tests a
signature: Wilcoxon signed-rank over the 12 cohorts on delta = observed AUC - median of the
target's null bank, called when P < 0.05. FDR for a library claim = (decoy call rate) x 52 /
(signatures called). Intervals: cluster bootstrap over the 52 targets, each carried with its
own decoys, so the dependence among one target's decoys is kept.
"""
import glob, numpy as np, pandas as pd
from scipy.stats import wilcoxon
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
RNG = np.random.default_rng(20261019)
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(D + "calib_v3/ICB_*_v3.csv"))])
NAT = sorted(glob.glob(D + "calib_v3/ICB_*_nat.csv"))          # natural-pool N1/N2 (null_natpool_v4.R)
if len(NAT) == 12: X = pd.concat([X, pd.concat([pd.read_csv(f) for f in NAT])])
NULLS = [n for n in ("N0", "N1", "N2", "N1nat", "N2nat") if n in set(X["null"])]
print(f"cohorts {X.cohort.nunique()} | targets {X[X.kind=='target'].unit.nunique()} | "
      f"decoys {X[X.kind=='decoy'].unit.nunique()}")

def _unit(args):
    (u, nl), g = args
    out = []
    for st in ("fold", "dir"):
        d = (g[st + "_obs"] - g[st + "_med"]).to_numpy()
        # scipy defaults, as in agg_calib_v1.py and the stored counts; a zero delta (discrete AUC
        # equal to the bank median) sends scipy to its slow exact-with-zeros path, hence the pool
        p = wilcoxon(d).pvalue if len(d) >= 6 and np.any(d != 0) else 1.0
        m = float(np.median(d))
        out.append(dict(unit=u, null=nl, stat=st, kind=g.kind.iloc[0], target=g.target.iloc[0],
                        pool=g.pool.iloc[0], structure=g.structure.iloc[0], ncoh=len(d),
                        med_delta=m, p=p, call=(p < .05) and m > 0, neg=(p < .05) and m < 0,
                        clear95=float(np.mean(g[st + "_obs"] > g[st + "_p95"]))))
    return out
import multiprocessing as mp
with mp.get_context("fork").Pool(32) as pool:
    rows = [r for rr in pool.map(_unit, list(X.groupby(["unit", "null"])), chunksize=40) for r in rr]
U = pd.DataFrame(rows); U.to_csv(D + "calib_v3_units.csv", index=False)

# ---- check: the targets re-read from the banks reproduce the stored counts
print("\nsignatures called, re-read from the banks (stored: 21 27 19 26 14 28):")
T = U[U.kind == "target"]
print("  " + "  ".join(f"{nl}-{st} {int(T[(T['null']==nl)&(T.stat==st)].call.sum())}"
                       for nl in NULLS for st in ("fold","dir")))

# ---- matching
R = X[X["null"] == "N0"]
tr = R[R.kind == "target"].groupby("cohort").r_int.median()
print("\ninternal r, median per cohort:   signatures | A coex | B coex | A rand")
for c in sorted(R.cohort.unique()):
    s = R[R.cohort == c]
    q = lambda p, st: s[(s.kind=="decoy")&(s.pool==p)&(s.structure==st)].r_int.median()
    print(f"  {c:18s} {tr[c]:.3f} | {q('A','coex'):.3f} | {q('B','coex'):.3f} | {q('A','rand'):.3f}")

# ---- FDR
DIRN = X.drop_duplicates("unit").set_index("unit").dirn
def fdr_block(nl, st, pool, struct, nb=6000):
    t = U[(U.kind=="target") & (U["null"]==nl) & (U.stat==st)].set_index("target")
    d = U[(U.kind=="decoy") & (U["null"]==nl) & (U.stat==st) & (U.pool==pool) & (U.structure==struct)]
    if st == "dir":                    # a decoy declared in its target's direction (target-decoy proper)
        d = d[d.unit.map(DIRN) == d.target.map(DIRN)]
    pi_t = d.groupby("target").call.mean().reindex(t.index)        # per-target decoy call rate
    tg = t.index.to_numpy(); call = t.call.to_numpy(float); pit = pi_t.to_numpy(float)
    def est(ix):
        Rn = call[ix].sum(); pi = np.nanmean(pit[ix]); E = pi * len(ix)
        return Rn, pi, E, (E / Rn if Rn else np.nan)
    Rn, pi, E, F = est(np.arange(len(tg)))
    bs = np.array([est(RNG.integers(0, len(tg), len(tg))) for _ in range(nb)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    return dict(null=nl, stat=st, pool=pool, decoys=struct, called=int(Rn), decoy_rate=100*pi,
                decoy_lo=100*lo[1], decoy_hi=100*hi[1], exp_false=E, exp_true=Rn-E,
                fdr=100*F, fdr_lo=100*lo[3], fdr_hi=100*hi[3], n_decoys=len(d))
F = pd.DataFrame([fdr_block(nl, st, p, s) for nl in NULLS for st in ("fold","dir")
                  for p, s in (("A","coex"), ("B","coex"), ("A","rand"))])
F.to_csv(D + "calib_v3_fdr.csv", index=False)
print("\ntarget-decoy FDR (cluster bootstrap over the 52 targets, 6,000 replicates; directed: decoys in their target's direction)")
print(f"{'null':4s} {'stat':4s} {'decoys':7s} {'called':>6s} | {'decoy call rate':>20s} | {'E[false]':>8s} {'E[true]':>7s} | {'FDR':>20s}")
for _, r in F.iterrows():
    print(f"{r['null']:4s} {r.stat:4s} {r.pool+' '+r.decoys:7s} {r.called:6d} | "
          f"{r.decoy_rate:5.1f}% ({r.decoy_lo:4.1f}-{r.decoy_hi:4.1f}) | {r.exp_false:8.1f} {r.exp_true:7.1f} | "
          f"{r.fdr:5.1f}% ({r.fdr_lo:4.1f}-{r.fdr_hi:4.1f})")

# ---- per signature: where it sits among its own 32 pool-A coex decoys
out = []
for nl, st in (("N0","fold"), ("N0","dir"), ("N2","dir")):
    t = U[(U.kind=="target") & (U["null"]==nl) & (U.stat==st)].set_index("target")
    d = U[(U.kind=="decoy") & (U["null"]==nl) & (U.stat==st) & (U.pool=="A") & (U.structure=="coex")]
    for s, r in t.iterrows():
        dd = d[d.target == s].med_delta.to_numpy()
        out.append(dict(signature=s, null=nl, stat=st, med_delta=r.med_delta, wilcoxon_p=r.p,
                        call=r.call, decoys=len(dd), decoys_at_or_above=int((dd >= r.med_delta).sum()),
                        p_decoy=(1 + (dd >= r.med_delta).sum()) / (1 + len(dd))))
P = pd.DataFrame(out); P.to_csv(D + "calib_v3_per_signature.csv", index=False)
print("\nsignatures beating every one of their 32 matched decoys (decoy-calibrated P = 1/33):")
for nl, st in (("N0","fold"), ("N0","dir"), ("N2","dir")):
    q = P[(P["null"]==nl) & (P.stat==st)]
    top = q[q.decoys_at_or_above == 0]
    print(f"  {nl}-{st}: {len(top)} of 52; called by Wilcoxon but beaten by >=2 decoys: "
          f"{int(((q.call) & (q.decoys_at_or_above >= 2)).sum())}")
