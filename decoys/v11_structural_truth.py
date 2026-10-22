"""Analysis 4 of validation/PLAN_round6_applications.md, as amended (AMENDMENT_2_round6_structural_truth.md):
the planted-signal replicates of round 5 read against a structural truth.

Each replicate's own labels, calls and decoy estimate are regenerated exactly as in v10_calib_sim.py signal (same
seeds, same order of draws) and checked against v10_calib_signal.csv. For that replicate's causal set and beta,
M further label sets are drawn from a separate stream; a unit's expected delta is its mean over them in each
cohort, then the median over the replicate's cohorts. A target is structurally true when its expected delta
exceeds the 95th percentile of its eligible decoys' expected deltas.

    python v11_structural_truth.py [R] [M]   -> icb/v11_structural_truth.csv, v11_structural_truth_summary.csv
"""
import sys, functools, numpy as np, pandas as pd
from v10_calib_lib import Library, D
print = functools.partial(print, flush=True)

NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 200
M = int(sys.argv[2]) if len(sys.argv) > 2 else 20
L = Library()
nT = len(L.targets)
FAM = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature")
FAMV = FAM.family.reindex(L.targets).to_numpy()
inflamed = [t for t, f in zip(L.targets, FAMV) if f == "Inflamed/IFN"]
others = [t for t, f in zip(L.targets, FAMV) if f != "Inflamed/IFN"]
poolB = [d for d, p, s in zip(L.decoys, L.d_pool, L.d_struct) if p == "B" and s == "coex"]
SCEN = [("P1 inflamed", 0.5, 12), ("P1 inflamed", 1.0, 12), ("P2 non-inflamed", 0.5, 12),
        ("P2 non-inflamed", 1.0, 12), ("P3 outside signatures", 0.5, 12), ("P3 outside signatures", 1.0, 12),
        ("P1 inflamed", 1.0, 6), ("P2 non-inflamed", 1.0, 6)]
REF = pd.read_csv(f"{D}/v10_calib_signal.csv")
REF = REF[REF.truth == "10%"].set_index(["scenario", "beta", "cohorts", "rep", "stat"])


def draw(c, G, beta, rng, ncol):
    gi = [c.gidx[g] for g in G if g in c.gidx]
    u = c.Z[gi].mean(0); u = (u - u.mean()) / u.std()
    Y = np.zeros((c.n, ncol), np.float32); k = int(c.y.sum())
    for j in range(ncol):
        key = beta * u + rng.gumbel(size=c.n)
        Y[np.argsort(-key)[:k], j] = 1
    return Y


rows = []
for si, (scen, beta, ncoh) in enumerate(SCEN):
    rng = np.random.default_rng(20261021 + si)
    for r in range(NREP):
        if scen.startswith("P3"):
            dstar = poolB[rng.integers(len(poolB))]; G = L.d_genes[L.di[dstar]]; true = np.zeros(nT, bool); excl = G
        else:
            src = inflamed if scen.startswith("P1") else others
            pick = rng.choice(src, 3, replace=False)
            G = set().union(*[L.t_genes[t] for t in pick]); excl = None
            frac = np.array([len(L.t_genes[t] & G) / len(L.t_genes[t]) for t in L.targets])
            true = frac >= 0.10
        keep = np.sort(rng.choice(12, ncoh, replace=False)) if ncoh < 12 else np.arange(12)  # before the labels, as round 5 ran
        Ys = [draw(c, G, beta, rng, 1) for c in L.C]
        mask = np.zeros((12, nT), bool); mask[keep] = True
        res = L.run(Ys, cohort_mask=mask)
        rng2 = np.random.default_rng(20261025 + 1000 * si + r)
        Yx = [draw(c, G, beta, rng2, M) for c in L.C]
        T, Dd = L.deltas(Yx)
        for st in ("fold", "dir"):
            tc, dc = res[st][0][:, 0], res[st][1][:, 0]
            m = L.decoy_mask(st, exclude_genes=excl)
            pit = L.pi(dc[:, None], m)[:, 0]
            ok = ~np.isnan(pit)
            Rn = int((tc & ok).sum()); E = float(np.nanmean(pit) * ok.sum())
            et = np.nanmedian(np.nanmean(T[st][keep], axis=2), axis=0)          # targets
            ed = np.nanmedian(np.nanmean(Dd[st][keep], axis=2), axis=0)         # decoys
            strue = np.zeros(nT, bool)
            for ti in range(nT):
                v = ed[m & (L.d_tidx == ti)]
                v = v[np.isfinite(v)]
                strue[ti] = len(v) > 0 and np.isfinite(et[ti]) and et[ti] > np.quantile(v, 0.95)
            ref = REF.loc[(scen, beta, ncoh, r, st)]
            assert int(ref.R) == int(tc.sum()) and int(ref.V) == int((tc & ~true).sum()), (scen, beta, ncoh, r, st)
            assert abs(float(ref["E coex A"]) - E) < 1e-9, (scen, r, st, ref["E coex A"], E)
            rows.append(dict(scenario=scen, beta=beta, cohorts=ncoh, rep=r, stat=st, R=Rn, E=E,
                             V_gene=int((tc & ok & ~true).sum()), V_struct=int((tc & ok & ~strue).sum()),
                             n_true_gene=int(true.sum()), n_true_struct=int(strue.sum()),
                             agree=float(np.mean(true == strue))))
        if (r + 1) % 25 == 0:
            print(f"  {scen}, beta {beta}, {ncoh} cohorts: {r + 1} replicates")
S = pd.DataFrame(rows); S.to_csv(f"{D}/v11_structural_truth.csv", index=False)
g = S.groupby(["scenario", "beta", "cohorts", "stat"])
out = g.apply(lambda x: pd.Series(dict(R=x.R.mean(), estimate=100 * x.E.sum() / x.R.sum(),
                                       fdp_gene=100 * x.V_gene.sum() / x.R.sum(),
                                       fdp_struct=100 * x.V_struct.sum() / x.R.sum(),
                                       n_true_gene=x.n_true_gene.mean(), n_true_struct=x.n_true_struct.mean(),
                                       agree=x.agree.mean()))).reset_index()
out.to_csv(f"{D}/v11_structural_truth_summary.csv", index=False)
pd.set_option("display.width", 200)
print(out.round(2).to_string())
