"""Shared engine for analyses 3-6 of validation/PLAN_calibration_round5.md: the target-decoy pipeline of
calib_score_v3.R + calib_v3_fdr.py, rerun on any relabelling of the 12 response cohorts.

The score store (v10_store.R, v10_store_json.R) holds, per cohort, the ranks across patients of the 52 targets,
the 2,912 decoys and the 52 x 1,000 N0 bank sets, plus gene-wise z-scores. Ranks do not depend on the labels, so
the AUC of every unit under a label matrix Y (patients x replicates) is one matrix product, and the whole pipeline
— bank medians per statistic, Wilcoxon signed-rank over cohorts, calls, decoy call rates, FDR — is vectorized
over replicates. At the observed labels it must reproduce calib_v3_units.csv.
"""
import glob, json, functools, numpy as np, pandas as pd
from scipy.stats import rankdata

D = "/path/to/analysis_root/icb"
ST = f"{D}/v10_store"
COHORTS = ["ICB_Braun", "ICB_Gide", "ICB_Hugo", "ICB_Jung", "ICB_Kim", "ICB_Limagne1", "ICB_Liu",
           "ICB_Mariathasan", "ICB_Miao1", "ICB_Puch", "ICB_Riaz", "ICB_Van_Allen"]


def _f32(path, nrow, ncol):
    a = np.fromfile(path, dtype="<f4")
    assert a.size == nrow * ncol, path
    return a.reshape(nrow, ncol)


class Cohort:
    def __init__(self, co):
        m = json.load(open(f"{ST}/{co}_meta.json"))
        self.name, self.n = co, int(m["n"])
        self.y = np.asarray(m["y"], int)
        self.targets, self.decoys = list(m["targets"]), list(m["decoys"])
        self.dirn_t = np.asarray(m["dirn_t"], int)
        # named vectors arrive from jsonlite as arrays in the order of m["decoys"]
        self.dirn_d = np.asarray(m["dirn_d"], int); assert len(self.dirn_d) == len(self.decoys)
        self.target_of = list(m["target_of"]); assert len(self.target_of) == len(self.decoys)
        self.genes = list(m["genes"])
        self.tgt_genes = {k: list(v) for k, v in m["tgt_genes"].items()}
        nT, nD = len(self.targets), len(self.decoys)
        self.RT = _f32(f"{ST}/{co}_tgt.f32", nT, self.n)
        self.RD = _f32(f"{ST}/{co}_dcy.f32", nD, self.n)
        self.RB = _f32(f"{ST}/{co}_n0.f32", nT * int(m["ndraw"]), self.n)
        self.ndraw = int(m["ndraw"])
        self.Z = _f32(f"{ST}/{co}_z.f32", len(self.genes), self.n)
        self.gidx = {g: i for i, g in enumerate(self.genes)}
        ti = {t: i for i, t in enumerate(self.targets)}
        self.d_tidx = np.array([ti[t] for t in self.target_of])

    def auc(self, R, Y):
        """Signed AUC of every row of the rank matrix R against every column of the label matrix Y."""
        n1 = float(Y[:, 0].sum()); n0 = self.n - n1
        assert np.all(Y.sum(0) == n1)
        # ranks are multiples of 1/2 and rank sums stay far below 2^23, so the float32 product is exact;
        # the division is done in float64, as R does it
        return ((R @ Y).astype(np.float64) - n1 * (n1 + 1) / 2) / (n1 * n0)

    def deltas(self, Y):
        """Folded and directed delta (observed AUC - bank median, per statistic) for targets and decoys."""
        Y = Y.astype(np.float32)
        aT, aD = self.auc(self.RT, Y), self.auc(self.RD, Y)
        aB = self.auc(self.RB, Y).reshape(len(self.targets), self.ndraw, Y.shape[1])
        mf = np.median(np.maximum(aB, 1 - aB), axis=1)            # folded bank median, per target
        mp = np.median(aB, axis=1)                                  # signed bank median; 1 - mp the other way
        tf = np.maximum(aT, 1 - aT) - mf
        td = np.where(self.dirn_t[:, None] == 1, aT - mp, (1 - aT) - (1 - mp))
        j = self.d_tidx
        df = np.maximum(aD, 1 - aD) - mf[j]
        dd = np.where(self.dirn_d[:, None] == 1, aD - mp[j], (1 - aD) - (1 - mp[j]))
        return dict(tf=tf, td=td, df=df, dd=dd)


# ---------------------------------------------------------------- exact Wilcoxon signed-rank, as scipy
# (signed_rank_p and signed_rank_p_tied are perm_inference.py's own, validated there against scipy 1.17)
@functools.lru_cache(maxsize=None)
def _signed_rank_cdf(n):
    cnt = np.zeros(n * (n + 1) // 2 + 1); cnt[0] = 1.0
    for r in range(1, n + 1):
        nxt = cnt.copy(); nxt[r:] += cnt[:len(cnt) - r]; cnt = nxt
    return np.cumsum(cnt) / 2.0 ** n


@functools.lru_cache(maxsize=None)
def _sign_matrix(n):
    return ((np.arange(2 ** n)[:, None] >> np.arange(n)[None, :]) & 1).astype(float)


def signed_rank_p(D2):
    """Two-sided exact P for each column of D2 (cohorts x columns); ties and zeros take scipy's route."""
    n = D2.shape[0]
    A = np.sort(np.abs(D2), axis=0)
    bad = (D2 == 0).any(0) | (np.diff(A, axis=0) == 0).any(0)
    R = rankdata(np.abs(D2), axis=0)
    wp = np.where(D2 > 0, R, 0.0).sum(0)
    T = np.minimum(wp, n * (n + 1) / 2.0 - wp)
    cdf = _signed_rank_cdf(n)
    out = np.minimum(1.0, 2.0 * cdf[np.clip(np.rint(T).astype(int), 0, len(cdf) - 1)])
    if bad.any():
        X = D2[:, bad]; nz = X != 0
        Rk = np.where(nz, rankdata(np.where(nz, np.abs(X), np.inf), axis=0), 0.0)
        o = np.ones(X.shape[1]); cnts = nz.sum(0)
        for k in np.unique(cnts):
            if k == 0: continue
            sel = np.flatnonzero(cnts == k)
            Rs = np.sort(Rk[:, sel], axis=0)[n - k:]
            w = np.where(X[:, sel] > 0, Rk[:, sel], 0.0).sum(0)
            Tt = np.minimum(w, Rs.sum(0) - w)
            W = _sign_matrix(int(k)) @ Rs
            o[sel] = np.minimum(1.0, 2.0 * (W <= Tt[None, :] + 1e-9).mean(0))
        out[bad] = o
    return out


def calls(delta, valid):
    """delta: cohorts x units x B; valid: cohorts x units. Call = P < 0.05 with positive median (>= 6 cohorts)."""
    C, U, B = delta.shape
    out = np.zeros((U, B), bool)
    nval = valid.sum(0)
    for k in np.unique(nval):
        if k < 6: continue
        us = np.flatnonzero(nval == k)
        X = np.stack([delta[valid[:, u], u, :] for u in us], axis=1)        # k x |us| x B
        X2 = X.reshape(k, -1)
        p = signed_rank_p(X2).reshape(len(us), B)
        allzero = (X == 0).all(0)
        p[allzero] = 1.0
        out[us] = (p < .05) & (np.median(X, axis=0) > 0)
    return out


# ---------------------------------------------------------------- library-level structure
class Library:
    """The 12 cohorts with targets and decoys aligned across them."""

    def __init__(self, cohorts=COHORTS):
        self.C = [Cohort(c) for c in cohorts]
        self.targets = sorted(set().union(*[c.targets for c in self.C]))
        self.decoys = sorted(set().union(*[c.decoys for c in self.C]))
        I = pd.read_csv(f"{ST}/decoy_info.csv").set_index("decoy").loc[self.decoys]
        self.info = I
        self.d_target = I.target.to_numpy(); self.d_pool = I.pool.to_numpy(); self.d_struct = I.structure.to_numpy()
        self.d_dirn = I.dirn.to_numpy(int)
        self.d_close = ((I.r_achieved_med - I.r_target_med).abs() <= 0.05).to_numpy()
        self.d_genes = [set(s.split(";")) for s in I.genes]
        tdir = {}
        for c in self.C:
            tdir.update(dict(zip(c.targets, c.dirn_t)))
        self.t_dirn = np.array([tdir[t] for t in self.targets])
        self.t_genes = {t: set().union(*[set(c.tgt_genes.get(t, [])) for c in self.C]) for t in self.targets}
        self.ti = {t: i for i, t in enumerate(self.targets)}
        self.di = {d: i for i, d in enumerate(self.decoys)}
        self.d_tidx = np.array([self.ti[t] for t in self.d_target])
        self.vt = np.zeros((len(self.C), len(self.targets)), bool)
        self.vd = np.zeros((len(self.C), len(self.decoys)), bool)
        for ci, c in enumerate(self.C):
            self.vt[ci, [self.ti[t] for t in c.targets]] = True
            self.vd[ci, [self.di[d] for d in c.decoys]] = True

    def run(self, Ys, cohort_mask=None):
        """Ys: list (per cohort) of label matrices n_c x B. Returns calls per statistic for targets and decoys."""
        T, Dd = self.deltas(Ys)
        vt, vd = self.vt.copy(), self.vd.copy()
        if cohort_mask is not None:                    # cohorts x targets (e.g. origin cohorts left out)
            vt &= cohort_mask
            vd &= cohort_mask[:, self.d_tidx]
        return {s: (calls(np.nan_to_num(T[s]), vt), calls(np.nan_to_num(Dd[s]), vd)) for s in ("fold", "dir")}

    def deltas(self, Ys):
        """Per-cohort deltas, cohorts x units x B, for targets and decoys (NaN where a unit is not scored)."""
        B = Ys[0].shape[1]
        T = {s: np.full((len(self.C), len(self.targets), B), np.nan) for s in ("fold", "dir")}
        Dd = {s: np.full((len(self.C), len(self.decoys), B), np.nan) for s in ("fold", "dir")}
        for ci, (c, Y) in enumerate(zip(self.C, Ys)):
            d = c.deltas(Y)
            it = [self.ti[t] for t in c.targets]; idd = [self.di[x] for x in c.decoys]
            T["fold"][ci, it], T["dir"][ci, it] = d["tf"], d["td"]
            Dd["fold"][ci, idd], Dd["dir"][ci, idd] = d["df"], d["dd"]
        return T, Dd

    def decoy_mask(self, stat, pool="A", struct="coex", close=False, exclude_genes=None):
        m = (self.d_pool == pool) & (self.d_struct == struct)
        if stat == "dir":
            m &= self.d_dirn == self.t_dirn[self.d_tidx]
        if close:
            m &= self.d_close
        if exclude_genes is not None:
            m &= np.array([not (g & exclude_genes) for g in self.d_genes])
        return m

    def pi(self, dcalls, mask):
        """Per-target decoy call rate (targets x B); NaN where a target has no eligible decoy."""
        B = dcalls.shape[1]
        num = np.zeros((len(self.targets), B)); den = np.zeros(len(self.targets))
        np.add.at(num, self.d_tidx[mask], dcalls[mask].astype(float))
        np.add.at(den, self.d_tidx[mask], 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return num / den[:, None]


def fdr_point(call, pit, ix=None):
    """calib_v3_robust.est: targets without an eligible decoy drop out of both terms. call, pit: targets (x B)."""
    if ix is not None:
        call, pit = call[ix], pit[ix]
    ok = ~np.isnan(pit)
    R = np.where(ok, call, 0).sum(0)
    E = np.where(ok, pit, 0).sum(0) / np.maximum(ok.sum(0), 1) * ok.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        F = np.where(R > 0, E / R, np.nan)
    return R, E, F
