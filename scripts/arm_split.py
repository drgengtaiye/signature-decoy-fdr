"""IMmotion150's checkpoint arm is two regimens: split them.

    python arm_split.py

The trial randomised to atezolizumab alone (86 with progression data), atezolizumab with bevacizumab
(88) or sunitinib (89). The rest of this paper pools the two atezolizumab arms, as the trial's own
biomarker report does; this script asks whether the within-arm null result and the treatment
interaction hold in each of them separately. Added after external review pointed out that they need
not, and they do not: the interaction is carried by the bevacizumab combination."""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L
from lifelines import CoxPHFitter
from scipy.stats import norm

P = L.P
LADDER = [4, 7, 12, 20, 32, 50, 80, 130, 260, 440]
INFLAMED = ["CD8A","GZMB","PRF1","IFNG","CXCL9","CXCL10","STAT1","TAP1","PSMB9","HLA-DRA","CD3D","CD3E"]
C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]
C, M = C.loc[keep], M[keep]
t, e = L._pfs(C)
ici = C["ICI_RX"].fillna("none").astype(str); non = C["NON_ICI_RX"].fillna("none").astype(str)
grp = pd.Series(np.where(ici.eq("Atezolizumab") & non.eq("Bevacizumab"), "atezo_bev",
                np.where(ici.eq("Atezolizumab"), "atezo_mono", "sunitinib")), index=C.index)
ok = t.notna() & (t > 0)
idx = C.index[ok.reindex(C.index).fillna(False)]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T[idx]
grp, t, e = grp[idx], t[idx], e[idx]
print({k: int(v) for k, v in grp.value_counts().items()}, "| events", int(e.sum()))

def coef(v, tt, ee):
    v = np.asarray(v, float); sd = v.std()
    if not np.isfinite(sd) or sd == 0: return np.nan
    b, _, o = L.cox_fit(((v - v.mean())/sd)[:, None], np.asarray(tt, float), np.asarray(ee, int))
    return float(b[0]) if o else np.nan

SIGS = L.signatures()
pool = list(Z.index); rng = np.random.default_rng(8)
rows_out = []
print("\nwithin-arm null, immune signatures clearing the 5th percentile in the benefit direction")
for g in ("atezo_mono", "atezo_bev", "sunitinib"):
    m = (grp == g).values
    tt, ee = t.values[m], e.values[m]
    null = {k: np.array([coef(Z.loc[rng.choice(pool, k, replace=False)].values[:, m].mean(0), tt, ee)
                         for _ in range(200)]) for k in (4, 7, 12, 20)}
    hits = tot = 0; hrs = []
    for nm, gl in SIGS.items():
        if nm not in L.IMM: continue
        gg = [x for x in gl if x in Z.index]
        if len(gg) < 3: continue
        k = min([4,7,12,20], key=lambda x: abs(x - len(gg)))
        c = coef(Z.loc[gg].values[:, m].mean(0), tt, ee)
        v = null[k]; tot += 1; hrs.append(c)
        hits += int(c < np.nanquantile(v, .05))
    print(f"  {g:11s} n={int(m.sum()):3d} events={int(ee.sum()):3d}  {hits:2d} of {tot}  median HR {np.exp(np.nanmedian(hrs)):.2f}")
    rows_out.append(dict(analysis="within_arm_null", arm=g, n=int(m.sum()), events=int(ee.sum()), clears=hits, tested=tot, median_hr=float(np.exp(np.nanmedian(hrs)))))

s = Z.loc[[x for x in INFLAMED if x in Z.index]].values.mean(0)
s = (s - s.mean()) / s.std()
print("\npooled 12-gene inflamed score, score x arm against sunitinib")
for g in ("atezo_mono", "atezo_bev", "atezo_mono+atezo_bev"):
    sel = grp.isin(g.split("+") + ["sunitinib"]).values
    d = pd.DataFrame(dict(t=t.values[sel], e=e.values[sel], score=s[sel],
                          arm=(~grp[sel].eq("sunitinib")).astype(int).values))
    d["score_arm"] = d.score * d.arm
    c = CoxPHFitter().fit(d, "t", "e")
    b, se = c.params_["score_arm"], c.standard_errors_["score_arm"]
    print(f"  {g:22s} n={len(d):3d}  score x arm {b:+.3f} (se {se:.3f})  HR {np.exp(b):.2f}  one-sided P = {norm.cdf(b/se):.4f}")
    rows_out.append(dict(analysis="interaction_vs_sunitinib", arm=g, n=len(d), coef=float(b), se=float(se), hr=float(np.exp(b)), p_one_sided=float(norm.cdf(b/se))))

pd.DataFrame(rows_out).to_csv("/groups/adv2105_gp/yichen/Yi/com/gastric/icb/arm_split.csv", index=False)
print("\nwrote arm_split.csv")
