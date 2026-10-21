"""Test C, the pre-specified permutation (validation/PLAN_three_tests.md, SHA-256 e2edd3db3cbcb4e4…).

    python tc2_perm.py   -> icb/tc_perm.csv

The arm label is permuted within each trial, 2,000 times. A permutation keeps each arm's size, so the
statistic is the difference in the fraction of genes clearing their own within-arm null without
subsampling; the size-matched fractions of tc1_gene_arms.py stay as the descriptive version. The null
cut is redrawn inside every permuted arm, from the same 200 random genes of that trial, so the
comparison is like for like.
"""
import sys, functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
import multiprocessing as mp
from multiprocessing import Pool
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L

STATE = {}          # filled in main(), inherited by the workers through fork

P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
NPERM, NRAND = 2000, 200
SIGS = L.signatures()
GI = sorted({g for s in SIGS if s in L.IMM for g in SIGS[s]})
GC = sorted({g for s in SIGS if s.startswith("CTRL_") for g in SIGS[s]})


def load(trial):
    if trial in ("IMmotion150", "PRINCE"):
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
        idx = arm.index[ok.reindex(arm.index).fillna(False)]
        Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T[idx]
        return dict(kind="surv", Z=Z, arm=arm[idx].to_numpy(int), t=t[idx].to_numpy(float), e=e[idx].to_numpy(int))
    M = pd.read_csv(f"{V}/cohorts/{trial}.expr.tsv.gz", sep="\t", index_col=0)
    cl = pd.read_csv(f"{V}/cohorts/{trial}.clin.tsv", sep="\t")
    cl = cl[cl.y.notna() & cl.arm.notna()].copy(); cl["sample"] = cl["sample"].astype(str)
    idx = [s for s in M.columns if s in set(cl["sample"])]
    cl = cl.set_index("sample").loc[idx]
    Z = ((M[idx].T - M[idx].mean(1)) / M[idx].std(1).replace(0, 1)).T
    return dict(kind="pcr", Z=Z, arm=cl.arm.to_numpy(int), y=cl.y.to_numpy(int))


def stats_matrix(d, X, sel):
    """benefit-direction statistic for every row of X inside the patient subset sel."""
    if d["kind"] == "pcr":
        y = d["y"][sel]; v = X[:, sel]
        r = np.apply_along_axis(lambda z: pd.Series(z).rank().to_numpy(), 1, v)
        n1, n0 = (y == 1).sum(), (y == 0).sum()
        if n1 == 0 or n0 == 0:
            return np.full(len(X), np.nan)
        return (r[:, y == 1].sum(1) - n1 * (n1 + 1) / 2) / (n1 * n0) - 0.5
    t, e = d["t"][sel], d["e"][sel]
    out = np.empty(len(X))
    for i in range(len(X)):
        b, _, ok = L.cox_fit(X[i, sel][:, None], t, e)
        out[i] = -b[0] if ok else np.nan
    return out


def one(args):
    trial, seed = args                       # seed None means the observed arms
    d = STATE[trial]
    arm = d["arm"] if seed is None else np.random.default_rng(seed).permutation(d["arm"])
    res = {}
    for setname in ("immune", "control"):
        X = d[setname]
        f = {}
        for label, val in (("checkpoint", 1), ("comparator", 0)):
            sel = np.where(arm == val)[0]
            null = stats_matrix(d, d["rand"], sel)
            cut = np.nanpercentile(null, 95)
            obs = stats_matrix(d, X, sel)
            f[label] = float(np.nanmean(obs > cut))
        res[setname] = f["checkpoint"] - f["comparator"]
    return res


TRIALS = ["IMmotion150", "PRINCE", "ISPY2_pembro", "ISPY2_durva"]


def main():
    for tr in TRIALS:
        d = load(tr)
        Z = d["Z"]
        pool = [g for g in Z.index if np.isfinite(Z.loc[g].to_numpy(float)).all()]
        rng = np.random.default_rng(sum(map(ord, tr)))
        d["rand"] = Z.loc[rng.choice(pool, NRAND, replace=False)].to_numpy(float)
        d["immune"] = Z.loc[[g for g in GI if g in Z.index]].to_numpy(float)
        d["control"] = Z.loc[[g for g in GC if g in Z.index]].to_numpy(float)
        d.pop("Z")
        STATE[tr] = d
        print(f"{tr}: {len(d['arm'])} patients, {int((d['arm']==1).sum())} on checkpoint, "
              f"{d['immune'].shape[0]} immune and {d['control'].shape[0]} control genes")

    rows = []
    for tr in TRIALS:
        obs = one((tr, None))
        with Pool(6) as pl:
            null = pl.map(one, [(tr, s) for s in range(NPERM)], chunksize=8)
        for setname in ("immune", "control"):
            nd = np.array([x[setname] for x in null])
            p_ = float((np.sum(nd >= obs[setname]) + 1) / (len(nd) + 1))
            rows.append(dict(trial=tr, gene_set=setname, observed=obs[setname],
                             null_mean=float(np.nanmean(nd)), p_one_sided=p_))
            print(f"  {tr:14s} {setname:8s} difference {obs[setname]:+.3f} "
                  f"(permutation mean {np.nanmean(nd):+.3f})  one-sided P = {p_:.4f}")
    R = pd.DataFrame(rows)
    R.to_csv(f"{D}/tc_perm.csv", index=False)
    imm = R[R.gene_set == "immune"]
    from scipy.stats import chi2
    chi = -2 * np.sum(np.log(imm.p_one_sided))
    print(f"\nPRIMARY across the four comparisons (Fisher combination of the immune sets): "
          f"chi2 = {chi:.1f}, P = {chi2.sf(chi, 2*len(imm)):.3g}")
    print(f"  mean difference {imm.observed.mean():+.3f}; positive in {int((imm.observed>0).sum())} of {len(imm)}")


if __name__ == "__main__":
    mp.set_start_method("fork", force=True)
    main()
