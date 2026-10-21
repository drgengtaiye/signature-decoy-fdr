"""Inference for the within-arm contrast (validation/PLAN_within_arm_inference.md,
SHA-256 4f2ac0fc4227852e…, deposited 2026-09-23 22:15 UTC).

    python wa1_perm.py   -> icb/wa_perm.csv

The statistic is the one the manuscript already reports: the difference between arms in the fraction of
the 17 immune signatures clearing their own within-arm null in the benefit direction. The arm label is
permuted 2,000 times within each trial and the whole procedure, null included, is repeated inside the
permuted arms.
"""
import sys, functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
import multiprocessing as mp
from multiprocessing import Pool
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L

P = L.P
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
LADDER = [4, 7, 12, 20, 32, 50, 80, 130, 260, 440]
NPERM, NRAND = 2000, 200
STATE = {}


def auc_rows(Y, X):
    """directional AUC of every row of X against the binary Y."""
    r = np.apply_along_axis(lambda z: pd.Series(z).rank().to_numpy(), 1, X)
    n1, n0 = (Y == 1).sum(), (Y == 0).sum()
    if n1 == 0 or n0 == 0:
        return np.full(len(X), np.nan)
    return (r[:, Y == 1].sum(1) - n1 * (n1 + 1) / 2) / (n1 * n0) - 0.5


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
        d = dict(kind="surv", arm=arm[idx].to_numpy(int), t=t[idx].to_numpy(float), e=e[idx].to_numpy(int))
    else:
        M = pd.read_csv(f"{V}/cohorts/{trial}.expr.tsv.gz", sep="\t", index_col=0)
        cl = pd.read_csv(f"{V}/cohorts/{trial}.clin.tsv", sep="\t")
        cl = cl[cl.y.notna() & cl.arm.notna()].copy(); cl["sample"] = cl["sample"].astype(str)
        idx = [s for s in M.columns if s in set(cl["sample"])]
        cl = cl.set_index("sample").loc[idx]
        Z = ((M[idx].T - M[idx].mean(1)) / M[idx].std(1).replace(0, 1)).T
        d = dict(kind="pcr", arm=cl.arm.to_numpy(int), y=cl.y.to_numpy(int))
    SIGS = L.signatures()
    pool = [g for g in Z.index if np.isfinite(Z.loc[g].to_numpy(float)).all()]
    rng = np.random.default_rng(sum(map(ord, trial)))
    scores, sizes, fams = [], [], []
    for nm, gl in SIGS.items():
        g = [x for x in gl if x in pool]
        if len(g) < 3:
            continue
        fam = "Immune" if nm in L.IMM else "Control" if nm.startswith("CTRL_") else "Other"
        if fam == "Other":
            continue
        scores.append(Z.loc[g].to_numpy(float).mean(0)); sizes.append(len(g)); fams.append(fam)
    d["S"] = np.vstack(scores); d["fam"] = np.array(fams)
    d["rung"] = np.array([min(LADDER, key=lambda k: abs(k - s)) for s in sizes])
    Zp = Z.loc[pool].to_numpy(float)
    d["null"] = {k: np.vstack([Zp[rng.choice(len(Zp), min(k, len(Zp)), replace=False)].mean(0)
                               for _ in range(NRAND)]) for k in sorted(set(d["rung"]))}
    return d


def stats_rows(d, X, sel):
    """The statistic the manuscript uses, benefit-positive: a Cox coefficient per standard deviation
    of the score inside the arm (v8_within_arm_null.py standardizes the score within the arm before
    fitting, and so does this), or a directional AUC for pathological complete response."""
    if d["kind"] == "pcr":
        return auc_rows(d["y"][sel], X[:, sel])
    t, e = d["t"][sel], d["e"][sel]
    out = np.empty(len(X))
    for i in range(len(X)):
        v = X[i, sel]
        sd = v.std()
        if not np.isfinite(sd) or sd == 0:
            out[i] = np.nan; continue
        b, _, ok = L.cox_fit(((v - v.mean()) / sd)[:, None], t, e)
        out[i] = -b[0] if ok else np.nan          # benefit is a negative log hazard ratio
    return out


def one(args):
    trial, seed = args
    d = STATE[trial]
    arm = d["arm"] if seed is None else np.random.default_rng(seed).permutation(d["arm"])
    res = {}
    for label, val in (("checkpoint", 1), ("comparator", 0)):
        sel = np.where(arm == val)[0]
        cuts = {k: np.nanpercentile(stats_rows(d, d["null"][k], sel), 95) for k in d["null"]}
        obs = stats_rows(d, d["S"], sel)
        clear = np.array([o > cuts[r] for o, r in zip(obs, d["rung"])])
        for fam in ("Immune", "Control"):
            m = d["fam"] == fam
            res[(label, fam)] = float(np.nanmean(clear[m])) if m.any() else np.nan
    return {fam: res[("checkpoint", fam)] - res[("comparator", fam)] for fam in ("Immune", "Control")} | \
           {f"{lab}_{fam}": res[(lab, fam)] for lab in ("checkpoint", "comparator") for fam in ("Immune", "Control")}


TRIALS = ["IMmotion150", "PRINCE", "ISPY2_pembro", "ISPY2_durva"]


def main():
    for tr in TRIALS:
        STATE[tr] = load(tr)
        d = STATE[tr]
        print(f"{tr}: {len(d['arm'])} patients, {int((d['arm']==1).sum())} on checkpoint, "
              f"{int((d['fam']=='Immune').sum())} immune and {int((d['fam']=='Control').sum())} control signatures, "
              f"null rungs {sorted(d['null'])}")
    # the observed values must be the ones the manuscript already reports, or the test is meaningless
    PUBLISHED = {("IMmotion150", "checkpoint"): 9, ("IMmotion150", "comparator"): 0,
                 ("PRINCE", "checkpoint"): 6, ("PRINCE", "comparator"): 0,
                 ("ISPY2_pembro", "checkpoint"): 14, ("ISPY2_pembro", "comparator"): 15,
                 ("ISPY2_durva", "checkpoint"): 6, ("ISPY2_durva", "comparator"): 4}
    rows = []
    for tr in TRIALS:
        obs = one((tr, None))
        for lab in ("checkpoint", "comparator"):
            got = round(17 * obs[f"{lab}_Immune"])
            want = PUBLISHED[(tr, lab)]
            assert got == want, f"{tr} {lab}: reproduced {got} of 17, the manuscript reports {want}"
        print(f"  {tr}: observed counts reproduce the manuscript", flush=True)
        with Pool(6) as pl:
            null = pl.map(one, [(tr, s) for s in range(NPERM)], chunksize=8)
        for fam in ("Immune", "Control"):
            nd = np.array([x[fam] for x in null])
            p_ = float((np.sum(nd >= obs[fam]) + 1) / (len(nd) + 1))
            rows.append(dict(trial=tr, family=fam, checkpoint=obs[f"checkpoint_{fam}"],
                             comparator=obs[f"comparator_{fam}"], difference=obs[fam],
                             null_mean=float(np.nanmean(nd)), p_one_sided=p_))
            print(f"  {tr:14s} {fam:8s} {obs[f'checkpoint_{fam}']:.3f} against {obs[f'comparator_{fam}']:.3f}, "
                  f"difference {obs[fam]:+.3f} (permutation mean {np.nanmean(nd):+.3f})  one-sided P = {p_:.4f}")
    R = pd.DataFrame(rows)
    R.to_csv(f"{D}/wa_perm.csv", index=False)
    from scipy.stats import chi2
    for fam in ("Immune", "Control"):
        d_ = R[R.family == fam]
        c = -2 * np.sum(np.log(d_.p_one_sided))
        print(f"\n{fam}: Fisher across the four comparisons chi2 = {c:.1f}, P = {chi2.sf(c, 2*len(d_)):.4g}")
    sv = R[(R.family == "Immune") & (R.trial.isin(["IMmotion150", "PRINCE"]))]
    c = -2 * np.sum(np.log(sv.p_one_sided))
    print(f"the two survival trials alone (pre-specified secondary): chi2 = {c:.1f}, P = {chi2.sf(c, 2*len(sv)):.4g}")


if __name__ == "__main__":
    mp.set_start_method("fork", force=True)
    main()
