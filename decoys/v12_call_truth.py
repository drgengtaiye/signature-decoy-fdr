"""Analyses 1, 3 and 6 of validation/PLAN_round7_estimand_and_bias.md on the response cohorts.

    python v12_call_truth.py signal <scenario index 0-7> [R] [M]   -> icb/v12_truth/signal_<si>.csv
    python v12_call_truth.py null [B]                             -> icb/v12_truth/null_pi0.csv

signal: each planted-signal replicate of round 5 is regenerated exactly as in v10_calib_sim.py signal (same seeds,
same order of draws; its calls and estimate checked against v10_calib_signal.csv). For its causal set and beta,
M further label sets are drawn from a separate stream (seed 710000 + 1000 x scenario + replicate) and the whole
call procedure is run on each, with the replicate's cohorts. A target is null when its call frequency over the M
sets is at most the mean call frequency of its eligible decoys (primary truth); the binomial variant, the
expected-delta truth of round 6 (with these M sets) and the gene-content truth of round 5 are recorded beside it.
For the replicate's own labels: the plain estimate, the pi0-corrected estimate (decoy-referenced Storey, lambda
0.5), and 95% intervals from 1,000 bootstrap resamples over targets and over the 13 clusters for both.

null: the 1,000 within-cohort label permutations of round 5 regenerated exactly (checked against
v10_calib_null.csv), with the pi0-corrected estimate added.
"""
import os, sys, functools, numpy as np, pandas as pd
from scipy.stats import binom
from v10_calib_lib import Library, signed_rank_p, D
print = functools.partial(print, flush=True)

OUTD = f"{D}/v12_truth"; os.makedirs(OUTD, exist_ok=True)
MODE = sys.argv[1]
L = Library()
nT = len(L.targets)
FAM = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature")
CL = FAM.cluster_fixed.reindex(L.targets).to_numpy()
FAMV = FAM.family.reindex(L.targets).to_numpy()
LAMBDA = 0.5


def pm(delta, valid):
    """As v10_calib_lib.calls, also returning P and the median delta (units x B); fewer than 6 cohorts: P 1, median 0."""
    C, U, B = delta.shape
    P = np.ones((U, B)); MED = np.zeros((U, B)); CALL = np.zeros((U, B), bool)
    nval = valid.sum(0)
    for k in np.unique(nval):
        if k < 6: continue
        us = np.flatnonzero(nval == k)
        X = np.stack([delta[valid[:, u], u, :] for u in us], axis=1)
        p = signed_rank_p(X.reshape(k, -1)).reshape(len(us), B)
        p[(X == 0).all(0)] = 1.0
        med = np.median(X, axis=0)
        P[us], MED[us], CALL[us] = p, med, (p < .05) & (med > 0)
    return P, MED, CALL


def run_pm(Ys, cohort_mask=None):
    T, Dd = L.deltas(Ys)
    vt, vd = L.vt.copy(), L.vd.copy()
    if cohort_mask is not None:
        vt &= cohort_mask; vd &= cohort_mask[:, L.d_tidx]
    return {s: (pm(np.nan_to_num(T[s]), vt), pm(np.nan_to_num(Dd[s]), vd)) for s in ("fold", "dir")}, T, Dd


def sstat(P, MED):
    return np.sign(MED) * -np.log10(np.maximum(P, 1e-300))


def qvals(st_t, st_d, m):
    """Each target's rank among its eligible decoys: (1 + #{s_d >= s_t}) / (1 + n). NaN without eligible decoys."""
    q = np.full(nT, np.nan)
    for t in range(nT):
        x = st_d[m & (L.d_tidx == t)]
        if len(x): q[t] = (1 + np.sum(x >= st_t[t] - 1e-12)) / (1 + len(x))
    return q


def estimates(call, pit, q):
    """Plain and pi0-corrected estimate (targets without eligible decoys drop out of every term)."""
    ok = ~np.isnan(pit)
    R = call[ok].sum(); E = pit[ok].sum(); n = ok.sum()
    pi0 = min(1.0, (q[ok] > LAMBDA).sum() / (n * (1 - LAMBDA)))
    return R, E, pi0


def boot(call, pit, q, groups, rng, nb=1000):
    """95% percentile intervals of F and F0 over nb resamples of targets (groups None) or of clusters."""
    ok = ~np.isnan(pit)
    v = np.stack([np.where(ok, call, 0), np.where(ok, pit, 0), ok, ok & (q > LAMBDA)], 1).astype(float)
    if groups is None:
        W = rng.multinomial(nT, np.full(nT, 1 / nT), size=nb) @ v
    else:
        u = np.unique(groups); S = np.stack([v[groups == g].sum(0) for g in u])
        W = rng.multinomial(len(u), np.full(len(u), 1 / len(u)), size=nb) @ S
    R, E, n, Z = W.T
    with np.errstate(all="ignore"):
        F = np.where(R > 0, E / R, np.nan); p0 = np.minimum(1, Z / (n * (1 - LAMBDA))); F0 = p0 * F
        return [np.nanpercentile(x, q_) for x in (F, F0) for q_ in (2.5, 97.5)]


if MODE == "signal":
    si = int(sys.argv[2]); NREP = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    M = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    inflamed = [t for t, f in zip(L.targets, FAMV) if f == "Inflamed/IFN"]
    others = [t for t, f in zip(L.targets, FAMV) if f != "Inflamed/IFN"]
    poolB = [d for d, p, s in zip(L.decoys, L.d_pool, L.d_struct) if p == "B" and s == "coex"]
    SCEN = [("P1 inflamed", 0.5, 12), ("P1 inflamed", 1.0, 12), ("P2 non-inflamed", 0.5, 12),
            ("P2 non-inflamed", 1.0, 12), ("P3 outside signatures", 0.5, 12), ("P3 outside signatures", 1.0, 12),
            ("P1 inflamed", 1.0, 6), ("P2 non-inflamed", 1.0, 6)]
    scen, beta, ncoh = SCEN[si]
    REF = pd.read_csv(f"{D}/v10_calib_signal.csv")
    REF = REF[(REF.truth == "10%") & (REF.scenario == scen) & (REF.beta == beta) & (REF.cohorts == ncoh)].set_index(["rep", "stat"])

    def draw(c, G, rng, ncol):
        gi = [c.gidx[g] for g in G if g in c.gidx]
        u = c.Z[gi].mean(0); u = (u - u.mean()) / u.std()
        Y = np.zeros((c.n, ncol), np.float32); k = int(c.y.sum())
        for j in range(ncol):
            key = beta * u + rng.gumbel(size=c.n)
            Y[np.argsort(-key)[:k], j] = 1
        return Y

    rng = np.random.default_rng(20261021 + si)
    rows = []
    for r in range(NREP):
        if scen.startswith("P3"):
            dstar = poolB[rng.integers(len(poolB))]; G = L.d_genes[L.di[dstar]]; true_gc = np.zeros(nT, bool); excl = G
        else:
            src = inflamed if scen.startswith("P1") else others
            pick = rng.choice(src, 3, replace=False)
            G = set().union(*[L.t_genes[t] for t in pick]); excl = None
            true_gc = np.array([len(L.t_genes[t] & G) / len(L.t_genes[t]) for t in L.targets]) >= 0.10
        keep = np.sort(rng.choice(12, ncoh, replace=False)) if ncoh < 12 else np.arange(12)
        Ys = [draw(c, G, rng, 1) for c in L.C]
        mask = np.zeros((12, nT), bool); mask[keep] = True
        obs, _, _ = run_pm(Ys, mask)
        rng2 = np.random.default_rng(710000 + 1000 * si + r)
        Yx = [draw(c, G, rng2, M) for c in L.C]
        tru, Tx, Dx = run_pm(Yx, mask)
        brng = np.random.default_rng(712000 + 1000 * si + r)
        for st in ("fold", "dir"):
            (Pt, Mt, tc), (Pd, Md, dc) = obs[st]
            tc, dc = tc[:, 0], dc[:, 0]
            m = L.decoy_mask(st, exclude_genes=excl)
            pit = L.pi(dc[:, None], m)[:, 0]
            ok = ~np.isnan(pit)
            q = qvals(sstat(Pt[:, 0], Mt[:, 0]), sstat(Pd[:, 0], Md[:, 0]), m)
            Rn, E, pi0 = estimates(tc, pit, q)
            ref = REF.loc[(r, st)]
            assert int(ref.R) == int(tc.sum()) and abs(float(ref["E coex A"]) - E) < 1e-9, (scen, r, st)
            # truths from the M further label sets
            ft = tru[st][0][2].mean(1)                                   # targets' call frequency
            pid = np.nanmean(L.pi(tru[st][1][2], m), axis=1)             # decoys' mean call frequency, per target
            kt = tru[st][0][2].sum(1)
            null_cp = ok & (ft <= pid)
            true_bin = ok & (binom.sf(kt - 1, M, np.clip(pid, 0, 1)) < 0.05)
            et = np.nanmedian(np.nanmean(Tx[st][keep], axis=2), axis=0)
            ed = np.nanmedian(np.nanmean(Dx[st][keep], axis=2), axis=0)
            true_ed = np.array([np.isfinite(et[t]) and (lambda v: len(v) > 0 and et[t] > np.quantile(v, 0.95))(
                ed[m & (L.d_tidx == t)][np.isfinite(ed[m & (L.d_tidx == t)])]) for t in range(nT)])
            lo_t = boot(tc, pit, q, None, brng); lo_c = boot(tc, pit, q, CL, brng)
            rows.append(dict(scenario=scen, beta=beta, cohorts=ncoh, rep=r, stat=st, R=int(Rn), E=E, pi0=pi0,
                             V_cp=int((tc & null_cp).sum()), V_bin=int((tc & ok & ~true_bin).sum()),
                             V_ed=int((tc & ok & ~true_ed).sum()), V_gc=int((tc & ok & ~true_gc).sum()),
                             n_true_cp=int((ok & ~null_cp).sum()), n_true_bin=int(true_bin.sum()),
                             n_true_ed=int((ok & true_ed).sum()), n_true_gc=int((ok & true_gc).sum()),
                             exp_V=float(ft[null_cp].sum()), pi_null=float(pid[null_cp].sum()),
                             pi_true=float(pid[ok & ~null_cp].sum()), excess_null=float((pid - ft)[null_cp].sum()),
                             F_lo_t=lo_t[0], F_hi_t=lo_t[1], F0_lo_t=lo_t[2], F0_hi_t=lo_t[3],
                             F_lo_c=lo_c[0], F_hi_c=lo_c[1], F0_lo_c=lo_c[2], F0_hi_c=lo_c[3]))
        if r == 0 or (r + 1) % 25 == 0:
            print(f"  {scen}, beta {beta}, {ncoh} cohorts: {r + 1} replicates")
    pd.DataFrame(rows).to_csv(f"{OUTD}/signal_{si}.csv", index=False)

elif MODE == "null":
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    REF = pd.read_csv(f"{D}/v10_calib_null.csv").set_index("perm")
    rng = np.random.default_rng(20261021)
    CH = 50; rows = []
    for lo in range(0, B, CH):
        b = min(CH, B - lo)
        Ys = [np.stack([rng.permutation(c.y) for _ in range(b)], 1).astype(np.float32) for c in L.C]
        res, _, _ = run_pm(Ys)
        for st in ("fold", "dir"):
            (Pt, Mt, tc), (Pd, Md, dc) = res[st]
            m = L.decoy_mask(st)
            pit = L.pi(dc, m)
            St, Sd = sstat(Pt, Mt), sstat(Pd, Md)
            for j in range(b):
                q = qvals(St[:, j], Sd[:, j], m)
                Rn, E, pi0 = estimates(tc[:, j], pit[:, j], q)
                ref = REF.loc[lo + j + 1]
                assert int(ref[f"{st}_coex A_Rok"]) == Rn and abs(ref[f"{st}_coex A_E"] - E) < 1e-9, (lo + j, st)
                rows.append(dict(perm=lo + j + 1, stat=st, R=int(Rn), E=E, pi0=pi0, E0=pi0 * E))
        print(f"  {lo + b} permutations")
    P = pd.DataFrame(rows); P.to_csv(f"{OUTD}/null_pi0.csv", index=False)
    brng = np.random.default_rng(713000)
    for st, x in P.groupby("stat"):
        rb = [x.E0.to_numpy()[i].mean() / x.R.to_numpy()[i].mean() for i in (brng.integers(0, len(x), len(x)) for _ in range(4000))]
        print(f"{st}: realized {x.R.mean():.3f}, plain {x.E.mean():.3f} (ratio {x.E.mean() / x.R.mean():.3f}), "
              f"pi0-corrected {x.E0.mean():.3f} (ratio {x.E0.mean() / x.R.mean():.3f}, "
              f"{np.percentile(rb, 2.5):.3f}-{np.percentile(rb, 97.5):.3f}); mean pi0 {x.pi0.mean():.3f}")
