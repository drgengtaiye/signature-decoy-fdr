"""Test C of validation/PLAN_three_tests.md (SHA-256 e2edd3db3cbcb4e4…, 2026-09-23 21:33 UTC):
the treatment specificity one gene at a time.

    python tc1_gene_arms.py   -> icb/tc_gene_arms.csv

Every gene of the inflamed/interferon family, and of the confounder controls, is fitted inside each arm
of each randomized comparison against 200 random genes of that trial in the same patients: a Cox
coefficient per s.d. for progression-free survival in IMmotion150 and PRINCE, a folded AUC for
pathological complete response in the two I-SPY2 contrasts. The larger arm is subsampled to the
smaller arm's size (40 draws) and those matched counts are the primary ones.
"""
import sys, functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L

P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
NRAND, NSUB = 200, 40
rng = np.random.default_rng(23)

SIGS = L.signatures()                      # the trial library, as the within-arm analysis used it
INF = [s for s in SIGS if s in L.IMM]      # the 17 immune sets that enter the interaction test
CTL = [s for s in SIGS if s.startswith("CTRL_")]
GI = sorted({g for s in INF for g in SIGS[s]})
GC = sorted({g for s in CTL for g in SIGS[s]})
print(f"immune genes {len(GI)} from {len(INF)} signatures | control genes {len(GC)} from {len(CTL)}")


def auc(y, v):
    y = np.asarray(y); v = np.asarray(v, float)
    ok = np.isfinite(v); y, v = y[ok], v[ok]
    if len(np.unique(y)) < 2:
        return np.nan
    r = pd.Series(v).rank().values
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def surv_trial(trial):
    C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0) if trial == "IMmotion150" \
        else pd.read_csv(f"{P}/v2_prince_clinical_api.csv", index_col=0)
    M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0) if trial == "IMmotion150" \
        else pd.read_csv(f"{P}/ici/paad_iatlas_prince_2022.csv", index_col=0)
    arm = (C["ICI_RX"].astype(str) == "Atezolizumab").astype(int) if trial == "IMmotion150" \
        else C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False).astype(int)
    keep = [s for s in M.columns if s in C.index]
    C, M, arm = C.loc[keep], M[keep], arm.loc[keep]
    t, e = L._pfs(C)
    ok = t.notna() & (t > 0) & arm.notna()
    Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
    idx = arm.index[ok.reindex(arm.index).fillna(False)]
    return dict(kind="surv", Z=Z[idx], arm=arm[idx].to_numpy(int), t=t[idx].to_numpy(float), e=e[idx].to_numpy(int))


def pcr_trial(name):
    M = pd.read_csv(f"{V}/cohorts/{name}.expr.tsv.gz", sep="\t", index_col=0)
    cl = pd.read_csv(f"{V}/cohorts/{name}.clin.tsv", sep="\t")
    cl = cl[cl.y.notna() & cl.arm.notna()]
    cl["sample"] = cl["sample"].astype(str)          # the matrix carries the ids as text
    idx = [s for s in M.columns if s in set(cl["sample"])]
    cl = cl.set_index("sample").loc[idx]
    Z = ((M[idx].T - M[idx].mean(1)) / M[idx].std(1).replace(0, 1)).T
    return dict(kind="pcr", Z=Z, arm=cl.arm.to_numpy(int), y=cl.y.to_numpy(int))


def stat_one(d, v, sel):
    """benefit-direction statistic for one gene inside one arm (higher is better)."""
    if d["kind"] == "surv":
        b, _, ok = L.cox_fit(v[sel][:, None], d["t"][sel], d["e"][sel])
        return -b[0] if ok else np.nan     # a negative log hazard ratio is benefit
    a = auc(d["y"][sel], v[sel])
    return a - 0.5


def run_arm(d, sel, genes):
    Z = d["Z"]
    pool = [g for g in Z.index if np.isfinite(Z.loc[g].to_numpy(float)).all()]
    null = np.array([stat_one(d, Z.loc[g].to_numpy(float), sel) for g in rng.choice(pool, NRAND, replace=False)])
    cut = np.nanpercentile(null, 95)
    if not np.isfinite(cut):
        return {}, np.nan
    obs = {g: stat_one(d, Z.loc[g].to_numpy(float), sel) for g in genes if g in Z.index}
    return obs, cut


rows = []
TRIALS = [("IMmotion150", surv_trial("IMmotion150")), ("PRINCE", surv_trial("PRINCE")),
          ("ISPY2_pembro", pcr_trial("ISPY2_pembro")), ("ISPY2_durva", pcr_trial("ISPY2_durva"))]
for name, d in TRIALS:
    n1, n0 = int((d["arm"] == 1).sum()), int((d["arm"] == 0).sum())
    m = min(n1, n0)
    for label, val in (("checkpoint", 1), ("comparator", 0)):
        full = np.where(d["arm"] == val)[0]
        for setname, genes in (("immune", GI), ("control", GC)):
            obs, cut = run_arm(d, full, genes)
            if not obs:
                continue
            frac = float(np.nanmean([float(v > cut) if np.isfinite(v) else np.nan for v in obs.values()]))
            sub = []
            for _ in range(NSUB):
                s = rng.choice(full, m, replace=False)
                o2, c2 = run_arm(d, s, genes)
                sub.append(np.nanmean([float(v > c2) if np.isfinite(v) else np.nan for v in o2.values()]))
            rows.append(dict(trial=name, arm=label, n_arm=len(full), n_matched=m, gene_set=setname,
                             n_genes=len(obs), frac_clear=frac, frac_matched=float(np.mean(sub))))
            print(f"  {name:14s} {label:11s} {setname:8s} n={len(full):4d} genes={len(obs):4d} "
                  f"clear {frac:.3f} | size-matched {np.mean(sub):.3f}", flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{D}/tc_gene_arms.csv", index=False)
print("\nwrote tc_gene_arms.csv")
