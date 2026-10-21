"""The treatment interaction pooled across the randomized trials
(validation/PLAN_pooled_interaction.md, SHA-256 75a039d77a299db0…, deposited 2026-09-23 22:36 UTC).

    python pi1_pooled.py   -> icb/pi_pooled.csv

One pre-specified score, the manuscript's 12-gene inflamed score, standardised within each trial, so
there is no multiplicity and no loss to the correlation between signatures. Model 1 pools the two
survival trials in a Cox model stratified by trial; model 2 pools the two I-SPY2 contrasts in a
logistic model with a trial indicator.
"""
import sys, functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
import statsmodels.api as sm
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L
from scipy.stats import norm

P = L.P
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
INFLAMED = ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "TAP1", "PSMB9", "HLA-DRA", "CD3D", "CD3E"]
CONTROLS = {k: v for k, v in L.signatures().items() if k.startswith("CTRL_")}
NPERM = 2000


def score_of(Z, genes):
    """Mean of gene-wise z over the genes the matrix carries, ignoring the few missing values
    (two cells in ISPY2_pembro), then standardised within the trial."""
    g = [x for x in genes if x in Z.index]
    s = np.nanmean(Z.loc[g].to_numpy(float), axis=0)
    if not np.isfinite(s).all():                      # a patient with no measured gene of the set
        s = np.where(np.isfinite(s), s, np.nanmean(s))
    return (s - s.mean()) / (s.std() if s.std() else 1), len(g)


def surv(trial):
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
    return Z, pd.DataFrame(dict(trial=trial, arm=arm[idx].to_numpy(int), t=t[idx].to_numpy(float),
                                e=e[idx].to_numpy(int)))


def pcr(trial):
    M = pd.read_csv(f"{V}/cohorts/{trial}.expr.tsv.gz", sep="\t", index_col=0)
    cl = pd.read_csv(f"{V}/cohorts/{trial}.clin.tsv", sep="\t")
    cl = cl[cl.y.notna() & cl.arm.notna()].copy(); cl["sample"] = cl["sample"].astype(str)
    idx = [s for s in M.columns if s in set(cl["sample"])]
    cl = cl.set_index("sample").loc[idx]
    Z = ((M[idx].T - M[idx].mean(1)) / M[idx].std(1).replace(0, 1)).T
    return Z, pd.DataFrame(dict(trial=trial, arm=cl.arm.to_numpy(int), y=cl.y.to_numpy(int)))


def cox_pool(frames, label):
    d = pd.concat(frames, ignore_index=True)
    d["score_arm"] = d.score * d.arm
    cph = CoxPHFitter().fit(d[["t", "e", "score", "arm", "score_arm", "trial"]], "t", "e",
                            strata=["trial"])
    b = cph.params_["score_arm"]; se = cph.standard_errors_["score_arm"]
    p1 = norm.cdf(b / se)                                     # benefit is a negative interaction
    print(f"{label}: n = {len(d)}, events {int(d.e.sum())} | score x arm {b:+.3f} (se {se:.3f}), "
          f"HR {np.exp(b):.2f} per s.d., one-sided P = {p1:.4g}")
    return d, b, se, p1


rows = []
FR, ZS = [], {}
for tr in ("IMmotion150", "PRINCE"):
    Z, d = surv(tr); s, k = score_of(Z, INFLAMED)
    d["score"] = s; FR.append(d); ZS[tr] = Z
    print(f"  {tr}: {len(d)} patients, {k} of 12 inflamed genes present")
d1, b1, se1, p1 = cox_pool(FR, "\nPRIMARY  two survival trials, Cox stratified by trial")

# permutation of the arm label within trial, as a check on the Wald P
rng = np.random.default_rng(23)
null = []
for _ in range(NPERM):
    dd = d1.copy()
    dd["arm"] = dd.groupby("trial").arm.transform(lambda a: rng.permutation(a.to_numpy()))
    dd["score_arm"] = dd.score * dd.arm
    try:
        c = CoxPHFitter().fit(dd[["t", "e", "score", "arm", "score_arm", "trial"]], "t", "e", strata=["trial"])
        null.append(c.params_["score_arm"])
    except Exception:
        null.append(np.nan)
null = np.array(null, float)
pperm = float((np.nansum(null <= b1) + 1) / (np.isfinite(null).sum() + 1))
print(f"         permutation of the arm label within trial: P = {pperm:.4f} ({np.isfinite(null).sum()} draws)")

FR2 = []
for tr in ("ISPY2_pembro", "ISPY2_durva"):
    Z, d = pcr(tr); s, k = score_of(Z, INFLAMED)
    d["score"] = s; FR2.append(d)
    print(f"  {tr}: {len(d)} patients, {k} of 12 inflamed genes present")
d2 = pd.concat(FR2, ignore_index=True)
X = pd.DataFrame(dict(score=d2.score, arm=d2.arm, score_arm=d2.score * d2.arm,
                      durva=(d2.trial == "ISPY2_durva").astype(int)))
ok2 = np.isfinite(X.to_numpy(float)).all(1) & np.isfinite(d2.y.to_numpy(float))
X, d2 = X[ok2], d2[ok2]
m2 = sm.Logit(d2.y.to_numpy(int), sm.add_constant(X)).fit(disp=0)
b2, se2 = m2.params["score_arm"], m2.bse["score_arm"]
p2 = norm.sf(b2 / se2)
print(f"\nSECONDARY  two I-SPY2 contrasts, logistic with a trial term: score x arm {b2:+.3f} "
      f"(se {se2:.3f}), one-sided P = {p2:.4g}")

z1, z2 = norm.isf(p1), norm.isf(p2)
w1, w2 = np.sqrt(len(d1)), np.sqrt(len(d2))
zc = (w1 * z1 + w2 * z2) / np.sqrt(w1 ** 2 + w2 ** 2)
print(f"           Stouffer over all four comparisons ({len(d1)+len(d2)} patients): z = {zc:.2f}, "
      f"one-sided P = {norm.sf(zc):.4g}")

print("\nCONTROL scores, where no interaction is expected")
for nm, genes in CONTROLS.items():
    fr = []
    for tr in ("IMmotion150", "PRINCE"):
        Z = ZS[tr]; d = FR[0 if tr == "IMmotion150" else 1].copy()
        s, k = score_of(Z, genes); d["score"] = s; fr.append(d)
    _, b, se, p = cox_pool(fr, f"  {nm}")
    rows.append(dict(model="control", score=nm, coef=b, se=se, p_one_sided=p))

rows += [dict(model="survival_pooled", score="inflamed12", coef=b1, se=se1, p_one_sided=p1, p_perm=pperm,
              n=len(d1)),
         dict(model="pcr_pooled", score="inflamed12", coef=b2, se=se2, p_one_sided=p2, n=len(d2)),
         dict(model="stouffer", score="inflamed12", coef=np.nan, se=np.nan, p_one_sided=float(norm.sf(zc)),
              n=len(d1) + len(d2))]
pd.DataFrame(rows).to_csv(f"{D}/pi_pooled.csv", index=False)
print("\nwrote pi_pooled.csv")
