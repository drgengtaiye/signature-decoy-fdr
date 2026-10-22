"""Analyses 5 and 6 of validation/PLAN_calibration_round5.md: dependence-aware intervals for the observed
target-decoy FDR, and the FDR with each signature's origin cohorts left out.

    python v10_fdr_intervals.py   -> icb/v10_fdr_intervals.csv, icb/v10_origin_out.csv

Input: the per-cohort reads of calib_score_v3.R and calib_read_natbank.R (calib_v3/ICB_*_v3.csv, *_nat.csv),
the same rows calib_v3_fdr.py tests. Calls: Wilcoxon signed-rank over cohorts, P < 0.05 with a positive median.
Intervals for each FDR:
  (a) targets resampled with their decoys (the published interval)
  (b) the 13 correlation clusters of Fig. 2d resampled, each with its targets and their decoys
  (c) decoys alone resampled within each target, targets fixed (Monte Carlo error of the decoys)
  (d) the 12 cohorts resampled, every call of every target and decoy recomputed
  (e) clusters and cohorts resampled together
plus (b) with clusters cut at 1 - |r| = 0.4, 0.5 and 0.6. Origin cohorts: signature_provenance.csv.
"""
import glob, functools, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from v10_calib_lib import calls
print = functools.partial(print, flush=True)

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
rng = np.random.default_rng(20261021)
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{D}/calib_v3/ICB_*_v3.csv"))
               + sorted(glob.glob(f"{D}/calib_v3/ICB_*_nat.csv"))])
COH = sorted(X.cohort.unique()); assert len(COH) == 12
TG = sorted(X[X.kind == "target"].unit.unique())
I = pd.read_csv(f"{D}/v10_store/decoy_info.csv").set_index("decoy")
close = ((I.r_achieved_med - I.r_target_med).abs() <= 0.05)
DIRN = X.drop_duplicates("unit").set_index("unit").dirn

# ---------------------------------------------------------------- clusters at three cuts (0.5 is Fig. 2d)
EFF = pd.read_csv(f"{D}/v4_standardised_effects.csv"); sigs = sorted(EFF.signature); S_IDX = {s: i for i, s in enumerate(sigs)}
zs = np.zeros((52, 52)); ws = np.zeros((52, 52))
for co in COH:
    S = pd.read_csv(f"{D}/scores/{co}_sig.csv.gz"); cols = [c for c in S.columns if c in S_IDX]
    r = np.clip(np.nan_to_num(np.corrcoef(S[cols].to_numpy(float), rowvar=False)), -0.999999, 0.999999)
    ix = [S_IDX[c] for c in cols]; w = max(len(S) - 3, 1)
    zs[np.ix_(ix, ix)] += np.arctanh(r) * w; ws[np.ix_(ix, ix)] += w
R = np.tanh(np.divide(zs, ws, out=np.zeros_like(zs), where=ws > 0)); np.fill_diagonal(R, 1)
dist = 1 - np.abs(R); np.fill_diagonal(dist, 0); dist = (dist + dist.T) / 2
LK = linkage(squareform(dist, checks=False), method="average")
CUTS = {c: pd.Series(fcluster(LK, t=1 - c, criterion="distance"), index=sigs).reindex(TG).to_numpy() for c in (0.4, 0.5, 0.6)}
ref = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature").cluster_fixed.reindex(TG).to_numpy()
assert (pd.crosstab(CUTS[0.5], ref).gt(0).sum(1) == 1).all() and len(set(ref)) == len(set(CUTS[0.5])), "Fig. 2d clusters"
print("clusters at cuts 0.4/0.5/0.6:", {c: len(set(v)) for c, v in CUTS.items()})


# ---------------------------------------------------------------- per-cohort deltas, units aligned
def block(nl, st, struct_pool=("A", "coex"), keep_close=False):
    x = X[X["null"] == nl]
    t = x[x.kind == "target"]
    d = x[(x.kind == "decoy") & (x.pool == struct_pool[0]) & (x.structure == struct_pool[1])]
    if st == "dir":
        d = d[d.unit.map(DIRN) == d.target.map(DIRN)]
    if keep_close:
        d = d[d.unit.map(close).fillna(False).astype(bool)]
    units = TG + sorted(d.unit.unique())
    tgt_of = np.array(TG + list(d.drop_duplicates("unit").set_index("unit").target.reindex(units[len(TG):])))
    xx = pd.concat([t, d]); xx["delta"] = xx[f"{st}_obs"] - xx[f"{st}_med"]
    M = xx.pivot_table(index="cohort", columns="unit", values="delta").reindex(index=COH, columns=units)
    return M.to_numpy(float), np.array(units), tgt_of


def fdr(call, kind_t, tgt_of, tsel=None):
    """call: units (x B); the first 52 units are the targets. tsel: multiset of target indices (resampling)."""
    tc = call[:52]; dc = call[52:]; tix = np.array([TG.index(t) for t in tgt_of[52:]])
    num = np.zeros((52,) + tc.shape[1:]); den = np.zeros(52)
    np.add.at(num, tix, dc.astype(float)); np.add.at(den, tix, 1.0)
    pit = num / np.where(den > 0, den, np.nan)[(...,) + (None,) * (tc.ndim - 1)]
    ix = np.arange(52) if tsel is None else tsel
    p = pit[ix]; ok = ~np.isnan(p)
    Rn = np.where(ok, tc[ix], 0).sum(0); E = np.nansum(p, axis=0)
    return Rn, E, np.where(Rn > 0, E / np.maximum(Rn, 1), np.nan)


def calls_of(M):                                   # M: cohorts x units (x B)
    if M.ndim == 2: M = M[:, :, None]
    v = np.isfinite(M)
    return calls(np.nan_to_num(M), v[:, :, 0] if np.all(v == v[:, :, :1]) else v[:, :, 0])


def pct(v):
    v = v[np.isfinite(v)]; return np.percentile(v, 2.5), np.percentile(v, 97.5)


rows = []
COMBOS = [("N0", "fold", False), ("N0", "dir", False), ("N0", "fold", True), ("N0", "dir", True), ("N2nat", "dir", False)]
NB, NBC = 6000, 2000
for nl, st, kc in COMBOS:
    M, units, tgt_of = block(nl, st, keep_close=kc)
    c0 = calls_of(M)[:, 0]
    R0, E0, F0 = fdr(c0, None, tgt_of)
    lab = f"{nl} {'folded' if st == 'fold' else 'directed'}{', closest-matched decoys' if kc else ''}"
    out = dict(analysis=lab, called=int(R0), exp_false=float(E0), fdr=100 * float(F0))
    # (a) targets
    Fa = np.array([fdr(c0, None, tgt_of, rng.integers(0, 52, 52))[2] for _ in range(NB)])
    # (b) clusters, three cuts
    Fb = {}
    for cut, cl in CUTS.items():
        u = np.unique(cl); mem = {g: np.flatnonzero(cl == g) for g in u}
        Fb[cut] = np.array([fdr(c0, None, tgt_of, np.concatenate([mem[g] for g in rng.choice(u, len(u))]))[2]
                            for _ in range(NB)])
    # (c) decoys within targets
    tix = np.array([TG.index(t) for t in tgt_of[52:]]); grp = {i: np.flatnonzero(tix == i) for i in range(52)}
    Fc = []
    for _ in range(NB):
        sel = np.concatenate([52 + rng.choice(g, len(g)) for g in grp.values() if len(g)])
        cc = np.concatenate([c0[:52], c0[sel]]); to = np.concatenate([tgt_of[:52], tgt_of[sel]])
        Fc.append(fdr(cc, None, to)[2])
    Fc = np.array(Fc)
    # (d) cohorts, (e) cohorts and clusters (cut 0.5)
    cl = CUTS[0.5]; u = np.unique(cl); mem = {g: np.flatnonzero(cl == g) for g in u}
    Fd, Fe = [], []
    for b0 in range(0, NBC, 100):
        idx = rng.integers(0, 12, (12, 100))
        Mb = M[idx]                                         # 12 x 100 x units -> cohorts x units x B
        Mb = np.transpose(Mb, (0, 2, 1))
        v = np.isfinite(Mb)
        cb = np.zeros((M.shape[1], 100), bool)
        for j in range(100):                               # validity differs by resample (CRMA_Shukla)
            cb[:, j] = calls(np.nan_to_num(Mb[:, :, j:j + 1]), v[:, :, j])[:, 0]
        for j in range(100):
            Fd.append(fdr(cb[:, j], None, tgt_of)[2])
            Fe.append(fdr(cb[:, j], None, tgt_of, np.concatenate([mem[g] for g in rng.choice(u, len(u))]))[2])
    Fd, Fe = np.array(Fd), np.array(Fe)
    for nm, arr in [("(a) targets", Fa), ("(b) clusters, cut 0.4", Fb[0.4]), ("(b) clusters, cut 0.5", Fb[0.5]),
                    ("(b) clusters, cut 0.6", Fb[0.6]), ("(c) decoys only", Fc), ("(d) cohorts", Fd),
                    ("(e) clusters x cohorts", Fe)]:
        lo, hi = pct(arr)
        rows.append(dict(out, interval=nm, lo=100 * lo, hi=100 * hi, n_resamples=int(np.isfinite(arr).sum())))
    print(f"{lab}: called {R0}, E {E0:.2f}, FDR {100 * F0:.1f}% | " +
          " | ".join(f"{r['interval']} {r['lo']:.1f}-{r['hi']:.1f}" for r in rows[-7:]))
pd.DataFrame(rows).to_csv(f"{D}/v10_fdr_intervals.csv", index=False)

# ---------------------------------------------------------------- analysis 6: origin cohorts left out
PV = pd.read_csv(f"{D}/signature_provenance.csv").set_index("signature")
CNAME = {c: c.replace("ICB_", "") for c in COH}
CONSERV = {"C-ECM_Chakravarthy": ["Hugo", "Riaz"], "STAT1_Care": ["Mariathasan"],
           "ImmunoScore_Hao": ["Van_Allen", "Hugo", "Riaz"], "PTEN-MITF_Cabrita": ["Gide", "Riaz"],
           "KDM5A_Wang": ["Liu", "Riaz"], "TLS_Cabrita": ["Van_Allen", "Gide", "Riaz"],
           "APM_Wang": ["Hugo", "Van_Allen"], "CYT_Davoli": ["Van_Allen"]}   # signature_provenance_notes.md


def drop_mask(rule):
    m = np.ones((12, 52), bool)
    for j, s in enumerate(TG):
        for i, c in enumerate(COH):
            code = str(PV.loc[s, CNAME[c]]) if pd.notna(PV.loc[s, CNAME[c]]) else ""
            if code in ("D", "C") or (rule == "D, C, V" and code == "V") or \
               (rule == "D, C + borderline" and CNAME[c] in CONSERV.get(s, [])):
                m[i, j] = False
    return m


orow = []
for rule in ("none", "D, C", "D, C + borderline", "D, C, V"):
    keep = drop_mask(rule) if rule != "none" else np.ones((12, 52), bool)
    for nl, st, kc in COMBOS[:2] + COMBOS[4:]:
        M, units, tgt_of = block(nl, st, keep_close=kc)
        tix = np.array([TG.index(t) for t in tgt_of])
        Mk = np.where(keep[:, tix], M, np.nan)
        v = np.isfinite(Mk)
        c = calls(np.nan_to_num(Mk)[:, :, None], v)[:, 0]
        Rn, E, F = fdr(c, None, tgt_of)
        tested = int((v[:, :52].sum(0) >= 6).sum())
        called = [TG[i] for i in np.flatnonzero(c[:52])]
        orow.append(dict(rule=rule, analysis=f"{nl} {st}", cohort_tests_removed=int((~keep).sum()),
                         signatures_affected=int((~keep).any(0).sum()), signatures_testable=tested,
                         called=int(Rn), exp_false=float(E), fdr=100 * float(F), called_list=";".join(called)))
        print(f"origin out [{rule}] {nl} {st}: removed {int((~keep).sum())} cohort tests over "
              f"{int((~keep).any(0).sum())} signatures; called {Rn}, E {E:.2f}, FDR {100 * F:.1f}%")
O = pd.DataFrame(orow); O.to_csv(f"{D}/v10_origin_out.csv", index=False)
base = O[O.rule == "none"].set_index("analysis").called_list
for _, r in O[O.rule != "none"].iterrows():
    a, b = set(base[r.analysis].split(";")), set(r.called_list.split(";")) - {""}
    print(f"  [{r.rule}] {r.analysis}: lost {sorted(a - b)}, gained {sorted(b - a)}")
