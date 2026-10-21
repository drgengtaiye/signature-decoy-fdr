"""Reviewer checks on the target-decoy FDR (calibration v3).

    python calib_v3_robust.py  -> calib_v3_fdr_diff.csv, calib_v3_fdr_matched.csv

1. Paired comparison: the same bootstrap draw of the 52 targets (each with its decoys) is used for
   every scoring/null combination, so the difference in FDR and in expected true discoveries
   between two combinations gets its own interval.
2. Matching sensitivity: decoys are under-matched for the most coherent signatures (internal r
   above 0.6). The FDR is recomputed (a) with only decoys whose median internal r is within 0.05
   of their target's and (b) without targets whose internal r exceeds 0.5.
"""
import glob, numpy as np, pandas as pd
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
U = pd.read_csv(D + "calib_v3_units.csv")
I = pd.read_csv(D + "calib_decoys_v3_info.csv").set_index("decoy")
X = pd.concat([pd.read_csv(f, usecols=["unit", "kind", "dirn"]) for f in sorted(glob.glob(D + "calib_v3/ICB_*_v3.csv"))])
DIRN = X.drop_duplicates("unit").set_index("unit").dirn
RT = I.groupby("target").r_target_med.first()
U = U[(U.kind == "target") | ((U.pool == "A") & (U.structure == "coex"))]
COMBOS = [("N0", "fold"), ("N0", "dir"), ("N1nat", "fold"), ("N1nat", "dir"), ("N2nat", "fold"), ("N2nat", "dir")]
TG = np.array(sorted(U[U.kind == "target"].target.unique()))

def arrays(nl, st, keep_decoy=None):
    t = U[(U.kind == "target") & (U["null"] == nl) & (U.stat == st)].set_index("target").reindex(TG)
    d = U[(U.kind == "decoy") & (U["null"] == nl) & (U.stat == st)]
    if st == "dir": d = d[d.unit.map(DIRN) == d.target.map(DIRN)]
    if keep_decoy is not None: d = d[d.unit.map(keep_decoy).fillna(False).astype(bool)]
    pit = d.groupby("target").call.mean().reindex(TG).to_numpy(float)
    return t.call.to_numpy(float), pit

def est(call, pit, ix):
    ok = ix[~np.isnan(pit[ix])]                     # targets left without a decoy drop out of both terms
    Rn = call[ok].sum(); E = pit[ok].mean() * len(ok)
    return np.array([Rn, E, Rn - E, E / Rn if Rn else np.nan])

RNG = np.random.default_rng(20261019)
NB = 6000
BOOT = [RNG.integers(0, len(TG), len(TG)) for _ in range(NB)]
ALL = np.arange(len(TG))

# ---- 1. paired differences
A = {c: arrays(*c) for c in COMBOS}
pt = {c: est(*A[c], ALL) for c in COMBOS}
bs = {c: np.array([est(*A[c], ix) for ix in BOOT]) for c in COMBOS}
rows = []
for ref, alt in ((("N0", "fold"), ("N0", "dir")), (("N0", "fold"), ("N2nat", "dir")),
                 (("N0", "fold"), ("N1nat", "dir")), (("N0", "fold"), ("N2nat", "fold")),
                 (("N0", "dir"), ("N2nat", "dir"))):
    dd = bs[alt] - bs[ref]
    for j, q in ((2, "exp_true"), (3, "fdr")):
        v = dd[:, j]; v = v[~np.isnan(v)]
        sc = 100 if q == "fdr" else 1
        rows.append(dict(ref="-".join(ref), alt="-".join(alt), quantity=q,
                         ref_value=sc * pt[ref][j], alt_value=sc * pt[alt][j],
                         diff=sc * (pt[alt][j] - pt[ref][j]),
                         lo=sc * np.percentile(v, 2.5), hi=sc * np.percentile(v, 97.5),
                         p_one_sided=float(np.mean(v <= 0) if q == "exp_true" else np.mean(v >= 0))))
Dd = pd.DataFrame(rows); Dd.to_csv(D + "calib_v3_fdr_diff.csv", index=False)
print("paired bootstrap (same 6,000 resamples of the 52 targets for both combinations)")
for _, r in Dd.iterrows():
    print(f"  {r.alt:10s} vs {r.ref:8s} {r.quantity:8s}: {r.alt_value:6.1f} vs {r.ref_value:6.1f}  "
          f"diff {r['diff']:+6.1f} ({r.lo:+6.1f} to {r.hi:+6.1f})  one-sided P {r.p_one_sided:.4f}")

# ---- 2. matching sensitivity
dev = (I.r_achieved_med - I.r_target_med)
close = (dev.abs() <= 0.05)
low_r = set(RT[RT <= 0.5].index)
print(f"\ndecoys within 0.05 of target r: {close[I.pool.eq('A') & I.structure.eq('coex')].mean():.1%}; "
      f"targets with r <= 0.5: {len(low_r)} of {len(RT)}")
rows = []
for nl, st in COMBOS:
    for label, keep, tmask in (("all", None, None), ("decoys within 0.05 r", close, None),
                               ("targets r <= 0.5", None, low_r)):
        call, pit = arrays(nl, st, keep)
        ix = ALL if tmask is None else np.array([i for i, t in enumerate(TG) if t in tmask])
        p = est(call, pit, ix)
        b = np.array([est(call, pit, ix[RNG.integers(0, len(ix), len(ix))]) for _ in range(NB)])
        lo, hi = np.nanpercentile(b[:, 3], [2.5, 97.5])
        ntg = int((~np.isnan(pit[ix])).sum())
        rows.append(dict(null=nl, stat=st, subset=label, targets=ntg, called=int(p[0]),
                         exp_false=p[1], exp_true=p[2], fdr=100 * p[3], fdr_lo=100 * lo, fdr_hi=100 * hi))
M = pd.DataFrame(rows); M.to_csv(D + "calib_v3_fdr_matched.csv", index=False)
print("\nFDR under stricter matching")
for _, r in M.iterrows():
    print(f"  {r['null']:6s} {r.stat:4s} {r.subset:22s} targets {r.targets:2d} called {r.called:2d} "
          f"E[false] {r.exp_false:4.1f} FDR {r.fdr:5.1f}% ({r.fdr_lo:4.1f}-{r.fdr_hi:4.1f})")

# ---- 3. paired differences with only close-matched decoys, overall and by target coherence
rows = []
B = {c: arrays(*c, keep_decoy=close) for c in COMBOS}
for label, tset in (("all targets", None), ("targets r <= 0.5", low_r), ("targets r > 0.5", set(TG) - low_r)):
    ix0 = ALL if tset is None else np.array([i for i, t in enumerate(TG) if t in tset])
    boots = [ix0[RNG.integers(0, len(ix0), len(ix0))] for _ in range(NB)]
    for ref, alt in ((("N0", "fold"), ("N0", "dir")), (("N0", "fold"), ("N2nat", "dir"))):
        pr, pa = est(*B[ref], ix0), est(*B[alt], ix0)
        dd = np.array([est(*B[alt], ix) - est(*B[ref], ix) for ix in boots])
        for j, q in ((2, "exp_true"), (3, "fdr")):
            v = dd[:, j]; v = v[~np.isnan(v)]; sc = 100 if q == "fdr" else 1
            rows.append(dict(subset=label, ref="-".join(ref), alt="-".join(alt), quantity=q,
                             ref_value=sc * pr[j], alt_value=sc * pa[j], diff=sc * (pa[j] - pr[j]),
                             lo=sc * np.percentile(v, 2.5), hi=sc * np.percentile(v, 97.5),
                             ref_called=int(pr[0]), alt_called=int(pa[0])))
C = pd.DataFrame(rows); C.to_csv(D + "calib_v3_fdr_diff_close.csv", index=False)
print("\npaired differences, decoys within 0.05 r of their target")
for _, r in C.iterrows():
    print(f"  {r.subset:16s} {r.alt:9s} vs {r.ref:7s} {r.quantity:8s}: {r.alt_value:6.1f} vs {r.ref_value:6.1f} "
          f"(called {r.alt_called} vs {r.ref_called})  diff {r['diff']:+6.1f} ({r.lo:+6.1f} to {r.hi:+6.1f})")
