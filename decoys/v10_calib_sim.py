"""Analyses 3 and 4 of validation/PLAN_calibration_round5.md: the target-decoy FDR estimator where the truth is
known.

    python v10_calib_sim.py repro            reproduce calib_v3_units.csv at the observed labels
    python v10_calib_sim.py null [B]         analysis 3: B within-cohort label permutations (all targets false)
    python v10_calib_sim.py signal [R]       analysis 4: planted-signal scenarios, R replicates each
        -> icb/v10_calib_null.csv, v10_calib_null_summary.csv, v10_calib_signal.csv, v10_calib_signal_summary.csv
"""
import sys, json, functools, numpy as np, pandas as pd
from v10_calib_lib import Library, fdr_point, D
print = functools.partial(print, flush=True)

MODE = sys.argv[1]
L = Library()
nT = len(L.targets)
FAM = pd.read_csv(f"{D}/v5_signature_clusters.csv").set_index("signature")
CL = FAM.cluster_fixed.reindex(L.targets).to_numpy()
FAMV = FAM.family.reindex(L.targets).to_numpy()
ESTIMATORS = {"coex A": dict(pool="A", struct="coex"), "random A": dict(pool="A", struct="rand"),
              "coex B": dict(pool="B", struct="coex"), "coex A, closest": dict(pool="A", struct="coex", close=True)}


def observed_Y():
    return [c.y[:, None].astype(np.float32) for c in L.C]


def summarize(res, B, exclude_genes=None):
    """Per replicate: target calls, and for each estimator the expected false calls and FDR."""
    out = {}
    for st in ("fold", "dir"):
        tc, dc = res[st]
        out[f"{st}|R"] = tc.sum(0)
        out[f"{st}|tcalls"] = tc
        for nm, kw in ESTIMATORS.items():
            m = L.decoy_mask(st, exclude_genes=exclude_genes, **kw)
            pit = L.pi(dc, m)
            R, E, F = fdr_point(tc.astype(float), pit)
            out[f"{st}|{nm}|E"], out[f"{st}|{nm}|F"], out[f"{st}|{nm}|pit"] = E, F, pit
            out[f"{st}|{nm}|Rok"] = R
            out[f"{st}|{nm}|rate"] = dc[m].mean(0)
    return out


def boot_ci(call, pit, groups, nb=2000, rng=None):
    """Percentile intervals for the FDR, resampling targets (groups=None) or clusters, per replicate column."""
    B = call.shape[1]; lo = np.full(B, np.nan); hi = np.full(B, np.nan)
    if groups is None:
        idx = [rng.integers(0, nT, nT) for _ in range(nb)]
    else:
        u = np.unique(groups); mem = {g: np.flatnonzero(groups == g) for g in u}
        idx = [np.concatenate([mem[g] for g in rng.choice(u, len(u))]) for _ in range(nb)]
    F = np.full((nb, B), np.nan)
    for i, ix in enumerate(idx):
        F[i] = fdr_point(call[ix], pit[ix])[2]
    with np.errstate(all="ignore"):
        lo, hi = np.nanpercentile(F, 2.5, axis=0), np.nanpercentile(F, 97.5, axis=0)
    return lo, hi


if MODE == "repro":
    res = L.run(observed_Y())
    U = pd.read_csv(f"{D}/calib_v3_units.csv"); U = U[U["null"] == "N0"]
    for st in ("fold", "dir"):
        tc, dc = res[st]
        ref_t = U[(U.kind == "target") & (U.stat == st)].set_index("unit").call.reindex(L.targets).to_numpy()
        ref_d = U[(U.kind == "decoy") & (U.stat == st)].set_index("unit").call.reindex(L.decoys).to_numpy()
        print(f"{st}: targets called {tc[:, 0].sum()} (stored {ref_t.sum()}), agreement {np.mean(tc[:, 0] == ref_t):.4f}; "
              f"decoys agreement {np.mean(dc[:, 0] == ref_d):.4f} ({(dc[:, 0] != ref_d).sum()} differ)")
        s = summarize(res, 1)
        for nm in ESTIMATORS:
            print(f"   {nm:16s} called {int(s[f'{st}|{nm}|Rok'][0])}  E {s[f'{st}|{nm}|E'][0]:.2f}  "
                  f"FDR {100 * s[f'{st}|{nm}|F'][0]:.1f}%  decoy rate {100 * s[f'{st}|{nm}|rate'][0]:.1f}%")

elif MODE == "null":
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    rng = np.random.default_rng(20261021)
    CH = 50; rows = []; percall = {st: np.zeros(nT) for st in ("fold", "dir")}
    pisum = {(st, nm): np.zeros(nT) for st in ("fold", "dir") for nm in ESTIMATORS}
    for lo in range(0, B, CH):
        b = min(CH, B - lo)
        Ys = [np.stack([rng.permutation(c.y) for _ in range(b)], 1).astype(np.float32) for c in L.C]
        s = summarize(L.run(Ys), b)
        for j in range(b):
            r = dict(perm=lo + j + 1)
            for st in ("fold", "dir"):
                r[f"{st}_R"] = int(s[f"{st}|R"][j])
                r[f"{st}_target_rate"] = s[f"{st}|tcalls"][:, j].mean()
                for nm in ESTIMATORS:
                    r[f"{st}_{nm}_E"] = s[f"{st}|{nm}|E"][j]; r[f"{st}_{nm}_rate"] = s[f"{st}|{nm}|rate"][j]
                    r[f"{st}_{nm}_Rok"] = s[f"{st}|{nm}|Rok"][j]
            rows.append(r)
        for st in ("fold", "dir"):
            percall[st] += s[f"{st}|tcalls"].sum(1)
            for nm in ESTIMATORS:
                pisum[(st, nm)] += np.nan_to_num(s[f"{st}|{nm}|pit"]).sum(1)
        print(f"  {lo + b} permutations")
    P = pd.DataFrame(rows); P.to_csv(f"{D}/v10_calib_null.csv", index=False)
    out = []
    brng = np.random.default_rng(7)
    for st in ("fold", "dir"):
        for nm in ESTIMATORS:
            E, R = P[f"{st}_{nm}_E"].to_numpy(), P[f"{st}_{nm}_Rok"].to_numpy()
            rho = E.mean() / R.mean()
            bb = [brng.integers(0, len(P), len(P)) for _ in range(4000)]
            rb = np.array([E[i].mean() / R[i].mean() for i in bb])
            out.append(dict(stat=st, estimator=nm, perms=len(P), mean_target_calls=R.mean(),
                            mean_estimated_false=E.mean(), rho=rho, rho_lo=np.percentile(rb, 2.5),
                            rho_hi=np.percentile(rb, 97.5), target_test_rate=100 * P[f"{st}_target_rate"].mean(),
                            decoy_test_rate=100 * P[f"{st}_{nm}_rate"].mean(),
                            per_target_r=float(np.corrcoef(percall[st] / len(P), pisum[(st, nm)] / len(P))[0, 1])))
    O = pd.DataFrame(out); O.to_csv(f"{D}/v10_calib_null_summary.csv", index=False)
    print(O.round(3).to_string())
    pt = pd.DataFrame({f"{st}_target_callrate": percall[st] / len(P) for st in ("fold", "dir")}, index=L.targets)
    for st in ("fold", "dir"):
        pt[f"{st}_decoy_callrate"] = pisum[(st, "coex A")] / len(P)
    pt.to_csv(f"{D}/v10_calib_null_per_target.csv")

elif MODE == "signal":
    NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    inflamed = [t for t, f in zip(L.targets, FAMV) if f == "Inflamed/IFN"]
    others = [t for t, f in zip(L.targets, FAMV) if f != "Inflamed/IFN"]
    poolB = [d for d, p, s in zip(L.decoys, L.d_pool, L.d_struct) if p == "B" and s == "coex"]
    SCEN = [("P1 inflamed", 0.5, 12), ("P1 inflamed", 1.0, 12), ("P2 non-inflamed", 0.5, 12),
            ("P2 non-inflamed", 1.0, 12), ("P3 outside signatures", 0.5, 12), ("P3 outside signatures", 1.0, 12),
            ("P1 inflamed", 1.0, 6), ("P2 non-inflamed", 1.0, 6)]
    rows = []
    for si, (scen, beta, ncoh) in enumerate(SCEN):
        rng = np.random.default_rng(20261021 + si)
        for r in range(NREP):
            if scen.startswith("P3"):
                dstar = poolB[rng.integers(len(poolB))]; G = L.d_genes[L.di[dstar]]; true = np.zeros(nT, bool)
                excl = G
            else:
                src = inflamed if scen.startswith("P1") else others
                pick = rng.choice(src, 3, replace=False)
                G = set().union(*[L.t_genes[t] for t in pick]); excl = None
                frac = np.array([len(L.t_genes[t] & G) / len(L.t_genes[t]) for t in L.targets])
                true = frac >= 0.10
                true1, true25 = frac > 0, frac >= 0.25
            if scen.startswith("P3"):
                true1 = true25 = true
            # the cohorts of a 6-cohort replicate are drawn before its labels, as in the run that produced
            # v10_calib_signal.csv (a later edit had moved this line below the label draws)
            keep = np.sort(rng.choice(12, ncoh, replace=False)) if ncoh < 12 else np.arange(12)
            Ys = []
            for c in L.C:
                gi = [c.gidx[g] for g in G if g in c.gidx]
                u = c.Z[gi].mean(0); u = (u - u.mean()) / u.std()
                key = beta * u + rng.gumbel(size=c.n)                 # sampling without replacement, weights exp(beta u)
                y = np.zeros(c.n, np.float32); y[np.argsort(-key)[:int(c.y.sum())]] = 1
                Ys.append(y[:, None])
            mask = np.zeros((12, nT), bool); mask[keep] = True
            res = L.run(Ys, cohort_mask=mask)
            s = summarize(res, 1, exclude_genes=excl)
            for st in ("fold", "dir"):
                tc = res[st][0][:, 0]
                for tdef, tv in (("10%", true), ("any gene", true1), ("25%", true25)):
                    Rn = tc.sum(); V = (tc & ~tv).sum()
                    row = dict(scenario=scen, beta=beta, cohorts=ncoh, rep=r, stat=st, truth=tdef, R=int(Rn),
                               V=int(V), FDP=V / Rn if Rn else 0.0, n_true=int(tv.sum()))
                    for nm in ESTIMATORS:
                        row[f"E {nm}"] = float(s[f"{st}|{nm}|E"][0]); row[f"R {nm}"] = int(s[f"{st}|{nm}|Rok"][0])
                        row[f"FDR {nm}"] = float(s[f"{st}|{nm}|F"][0]) if s[f"{st}|{nm}|Rok"][0] else 0.0
                        # realized FDP over the targets this estimator can see (closest-matched drops some)
                        okt = ~np.isnan(s[f"{st}|{nm}|pit"][:, 0])
                        Rk = (tc & okt).sum(); row[f"FDP {nm}"] = (tc & okt & ~tv).sum() / Rk if Rk else 0.0
                    if tdef == "10%":
                        m = L.decoy_mask(st, exclude_genes=excl)
                        pit = L.pi(res[st][1], m)
                        brng = np.random.default_rng(1000 * si + r)
                        a, b = boot_ci(tc[:, None].astype(float), pit, None, rng=brng)
                        c_, d_ = boot_ci(tc[:, None].astype(float), pit, CL, rng=brng)
                        fdp = row["FDP coex A"]
                        row.update(ci_target_lo=a[0], ci_target_hi=b[0], ci_cluster_lo=c_[0], ci_cluster_hi=d_[0],
                                   cover_target=bool(a[0] <= fdp <= b[0]) if Rn else np.nan,
                                   cover_cluster=bool(c_[0] <= fdp <= d_[0]) if Rn else np.nan)
                    rows.append(row)
            if (r + 1) % 50 == 0:
                print(f"  {scen}, beta {beta}, {ncoh} cohorts: {r + 1} replicates")
    S = pd.DataFrame(rows); S.to_csv(f"{D}/v10_calib_signal.csv", index=False)
    g = S.groupby(["scenario", "beta", "cohorts", "stat", "truth"])
    out = g.agg(R=("R", "mean"), V=("V", "mean"), FDP=("FDP", "mean"), n_true=("n_true", "mean"),
                **{f"E_{k}": (f"E {k}", "mean") for k in ESTIMATORS},
                **{f"FDR_{k}": (f"FDR {k}", "mean") for k in ESTIMATORS},
                **{f"FDPv_{k}": (f"FDP {k}", "mean") for k in ESTIMATORS},
                cover_target=("cover_target", "mean"), cover_cluster=("cover_cluster", "mean")).reset_index()
    out.to_csv(f"{D}/v10_calib_signal_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[out.truth == "10%"][["scenario", "beta", "cohorts", "stat", "R", "V", "FDP", "FDR_coex A",
                                    "FDR_random A", "FDR_coex A, closest", "cover_target", "cover_cluster"]].round(3).to_string())
