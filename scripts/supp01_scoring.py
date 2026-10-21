"""Supplementary Figure 1 - the shift does not depend on how signatures were scored.

a  each signature's shift under its own publication's scoring against uniform z-mean scoring
b  median shift by scoring engine, inflamed family against all others
c  share of signatures clearing the matched null in each checkpoint cohort
d  median response gain per signature against the earlier fixed 48-gene null and the matched null
e  the shift against gene-set size
f  the inflamed-minus-confounder difference with each checkpoint cohort or cancer type removed
"""
import warnings; warnings.filterwarnings("ignore")
from statsmodels.stats.proportion import proportion_confint
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.formula.api as smf
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D
X = pd.read_csv(f"{D}/icb/v3_dissociation.csv")
X = X.merge(pd.read_csv(f"{D}/icb/v4_standardised_effects.csv")[["signature", "z_resp", "z_prog", "z_shift"]], on="signature")
Cc = pd.read_csv(f"{D}/icb/dissociation_curated.csv")[["signature", "shift"]]
mg = X[["signature", "shift", "family", "method"]].merge(Cc, on="signature", suffixes=("_native", "_zmean"))
R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med; R["clear"] = R.AUC > R.rand_p95; R["z"] = (R.AUC - R.rand_med) / R.rand_sd
O = pd.read_csv(f"{D}/icb/native_response.csv"); O["delta_old"] = O.AUC - O.rand_med
PP = pd.read_csv(f"{D}/icb/native_prognosis_matched.csv"); PP = PP[PP.endpoint == "OS"].copy()
PP["z"] = (PP.C - PP.rand_med) / PP.rand_sd
famd = X.set_index("signature").family
ORD3 = ["Other", "Confounder", "Inflamed/IFN"]

C = K.Canvas(100)

# ---- a  native against uniform scoring
ax = C.ax(12, 8, 34, 34)
for f in ORD3:
    g = mg[mg.family == f]
    ax.scatter(g.shift_zmean, g.shift_native, s=12, c=fs.FAM[f], lw=.3, ec="white", zorder=3)
lim = [-.065, .17]
ax.plot(lim, lim, ls=(0, (3, 2)), color=fs.INK, lw=.5, zorder=1)
ax.axhline(0, color=fs.LIGHT, lw=.6); ax.axvline(0, color=fs.LIGHT, lw=.6)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xticks([0, .1]); ax.set_yticks([0, .1])
ax.set_xlabel("shift, uniform z-mean"); ax.set_ylabel("shift, each paper's own method")
ra, pa = spearmanr(mg.shift_native, mg.shift_zmean)
K.key(ax, [(fs.BLUE, "inflamed/IFN"), (fs.GREY, "other"), (fs.VERM, "confounder control")], x=.05, y=.95, dy=.1, s=10)
K.frame(ax)
C.letter(1.5, 3, "a")

# ---- b  every scoring engine
ax = C.ax(70, 8, 36, 34)
meths = ["GSVA", "ssGSEA", "Weighted Mean"]
inf = [mg[(mg.method == m) & (mg.family == "Inflamed/IFN")].shift_native.median() for m in meths]
oth = [mg[(mg.method == m) & (mg.family != "Inflamed/IFN")].shift_native.median() for m in meths]
w = .36; x = np.arange(3)
for i_, m_ in enumerate(meths):                      # every signature, with the family medians
    for k_, (sel, col_) in enumerate([(mg.family == "Inflamed/IFN", fs.BLUE), (mg.family != "Inflamed/IFN", fs.GREY)]):
        v_ = mg[(mg.method == m_) & sel].shift_native.to_numpy(); xc = i_ + (k_ - .5) * w
        ax.scatter(xc + K.swarm(v_, .13, 9) * .8, v_, s=5, c=col_, lw=.2, ec="white", alpha=.8, zorder=2)
        ax.plot([xc - .14, xc + .14], [np.median(v_)] * 2, color=fs.INK if k_ == 0 else fs.SLATE, lw=1.2, zorder=4)
ax.axhline(0, color=fs.INK, lw=.5)
ax.set_xticks(x); ax.set_xticklabels(["GSVA", "ssGSEA", "weighted\nmean"], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylabel("shift (native units)"); ax.set_ylim(-.03, .13)
K.key(ax, [(fs.BLUE, "inflamed/IFN"), (fs.GREY, "all others")], x=.55, y=.95, dy=.1, marker="s", s=10)
K.frame(ax)
C.letter(56, 3, "b")

# ---- c  per checkpoint cohort
ax = C.ax(140, 8, 38, 34)
pc = (R.groupby("cohort").apply(lambda g: pd.Series({"n": g.n.iloc[0], "beat": g.clear.mean() * 100}), include_groups=False)
        .reset_index().sort_values("beat"))
y = np.arange(len(pc))
kn = R.groupby("cohort").clear.agg(["sum", "size"]).loc[pc.cohort]
ci_ = np.array([proportion_confint(int(a_), int(b_), method="wilson") for a_, b_ in zip(kn["sum"], kn["size"])]) * 100
ax.hlines(y, ci_[:, 0], ci_[:, 1], color=[fs.TUMCOL[fs.ICB_TUM[c]] for c in pc.cohort], lw=.9)
ax.scatter(pc.beat, y, s=14, c=[fs.TUMCOL[fs.ICB_TUM[c]] for c in pc.cohort], lw=.4, ec="white", zorder=3)
ax.axvline(5, color=fs.INK, ls=(0, (3, 2)), lw=.55)
ax.set_yticks(y); ax.set_yticklabels([f"{r.cohort.replace('ICB_', '').replace('_', ' ')} ({int(r.n)})" for _, r in pc.iterrows()],
                                     fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0); ax.set_xlabel("signatures clearing the null (%)")
K.frame(ax)
C.letter(112, 3, "c")

# ---- d  the shift under every scoring engine, in null s.d.
ax = C.ax(40, 66, 70, 26)
meths = [("GSVA", "GSVA"), ("ssGSEA", "ssGSEA"), ("Weighted Mean", "weighted mean")]
rng = np.random.default_rng(3)
eng = {}
for i, (m_, ml) in enumerate(meths):
    g = X[X.method == m_]
    for f in ORD3:
        gg = g[g.family == f]
        ax.scatter(gg.z_shift, np.full(len(gg), i) + rng.uniform(-.2, .2, len(gg)), s=10, c=fs.FAM[f], lw=.3, ec="white", zorder=3)
    ax.plot([g.z_shift.median()] * 2, [i - .33, i + .33], color=fs.INK, lw=1.2, zorder=4)
    from scipy.stats import wilcoxon
    eng[ml] = (len(g), round(g.z_shift.median(), 2), wilcoxon(g.z_shift).pvalue)
ax.axvline(0, color=fs.LIGHT, lw=.7, zorder=1)
ax.set_yticks(range(3)); ax.set_yticklabels([f"{ml} ({eng[ml][0]})" for _, ml in meths], fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0); ax.set_ylim(2.6, -.6); ax.set_xlim(-1.3, 2.6); ax.set_xticks([-1, 0, 1, 2])
ax.set_xlabel("shift, response − prognosis (null s.d.)"); K.frame(ax, left=False)
C.letter(1.5, 60, "d")

C.save("SF_scoring")
print(f"S1 | a n {len(mg)} rho {ra:+.2f} P {pa:.1e} | b inflamed {[round(v, 3) for v in inf]} others {[round(v, 3) for v in oth]}")
print("  d engines", eng)
print(pc.round(1).to_string(index=False))
