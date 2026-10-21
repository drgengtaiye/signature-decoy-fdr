"""Supplementary figure - the gastric discovery cohorts in detail.

a  pooled hazard ratio per s.d. across the three gastric cohorts, the immune axis and its baselines
b  discrimination added to stage plus age in each cohort
c  AUC for pembrolizumab response in the gastric checkpoint cohort, all markers
d  how much of a combined model's AUC is optimism: apparent against five-fold cross-validated
e  TCGA-STAD: each signature's C-index against its own matched null and against stage plus age
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D; RR = f"{D}/revision/results"
meta = pd.read_csv(f"{RR}/B4_meta.csv"); ci = pd.read_csv(f"{RR}/B4_cindex.csv")
auc = pd.read_csv(f"{RR}/A2_auc_ci.csv"); cv = pd.read_csv(f"{RR}/A3_cv_combo.csv")
SHOW = {"ESTIMATE(imm-str)": "ESTIMATE (immune − stromal)", "NetImmune(ssGSEA)": "NetImmune (ssGSEA)", "Barrier": "BarrierIndex"}
C = K.Canvas(178)

# ---- a  pooled hazard ratios
ax = C.ax(42, 10, 44, 30)
order = ["ESTIMATE(imm-str)", "NetImmune", "NetImmune(ssGSEA)", "Immune", "CD8A", "Barrier"]
mm = meta.set_index("score").loc[order].reset_index(); y = np.arange(len(mm))[::-1]
for i, (_, r) in enumerate(mm.iterrows()):
    c = fs.VERM if r.HR > 1 else (fs.BLUE if r["score"] == "NetImmune" else fs.SLATE)
    ax.plot([r.lo, r.hi], [y[i]] * 2, color=c, lw=1.0); ax.scatter(r.HR, y[i], s=16 if r["score"] == "NetImmune" else 11, c=c, zorder=3, lw=0)
ax.axvline(1, color=fs.INK, ls=(0, (3, 2)), lw=.55)
ax.set_yticks(y); ax.set_yticklabels([SHOW.get(s, s) for s in mm["score"]], fontsize=K.T_SMALL); ax.tick_params(axis="y", length=0)
for t, s in zip(ax.get_yticklabels(), mm["score"]):
    if s == "NetImmune": t.set_fontweight("bold")
ax.set_xlabel("pooled HR per s.d. (95% CI)"); ax.set_xlim(.6, 1.45); ax.set_ylim(-.6, len(mm) - .4); K.frame(ax, left=False)
C.letter(1.5, 4, "a")

# ---- b  increment over stage + age
ax = C.ax(110, 10, 66, 30)
coh = ["TCGA", "GSE84437", "GSE62254"]; scs = ["NetImmune", "ESTIMATE(imm-str)", "CD8A"]; cols = [fs.BLUE, fs.TEAL, fs.GREY]
w = .25
for j, s in enumerate(scs):
    g = ci[ci.score == s].set_index("cohort").loc[coh]; x = np.arange(3) + (j - 1) * w
    ax.vlines(x, g.dC_lo, g.dC_hi, color=cols[j], lw=1.1, zorder=2)
    ax.scatter(x, g.dC, s=22, c=cols[j], marker=["o", "D", "s"][j], lw=.5, ec="white", zorder=3)
ax.axhline(0, color=fs.INK, lw=.5)
ax.set_xticks(range(3)); ax.set_xticklabels(["TCGA-STAD", "GSE84437", "GSE62254"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("ΔC-index over stage + age"); lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi * 1.45)
for j, (c, s) in enumerate(zip(cols, scs)): K.key(ax, [(c, SHOW.get(s, s))], x=.40, y=.95 - .12 * j, marker=["o", "D", "s"][j], s=12)
K.frame(ax)
C.letter(96, 4, "b")

# ---- c  response in the gastric checkpoint cohort
ax = C.ax(30, 58, 48, 48)
a2 = auc.sort_values("AUC"); y = np.arange(len(a2))
for i, (_, r) in enumerate(a2.iterrows()):
    isn = r.Marker == "NetImmune"; c = fs.BLUE if isn else fs.SLATE
    ax.plot([r.CI_low, r.CI_high], [y[i]] * 2, color=c, lw=1.0); ax.scatter(r.AUC, y[i], s=16 if isn else 10, c=c, zorder=3, lw=0)
ax.axvline(.5, color=fs.INK, ls=(0, (3, 2)), lw=.55)
MK = {"TAM_SUP": "TAM module", "T_CELL": "T cell module", "TLS_B": "TLS/B module", "Barrier": "BarrierIndex", "TIDE(neg)": "TIDE (negated)", "GEP18": "GEP-18", "IFN": "IFN module", "APM": "APM module"}
ax.set_yticks(y); ax.set_yticklabels([MK.get(m_, m_) for m_ in a2.Marker], fontsize=K.T_SMALL); ax.tick_params(axis="y", length=0)
for t, s in zip(ax.get_yticklabels(), a2.Marker):
    if s == "NetImmune": t.set_fontweight("bold")
ax.set_xlabel("AUC, pembrolizumab response (95% CI)"); ax.set_xlim(.25, 1.0); ax.set_ylim(-.7, len(a2) - .3); K.frame(ax, left=False)
C.letter(1.5, 52, "c")

# ---- d  optimism
ax = C.ax(110, 58, 66, 48)
x = np.arange(len(cv)); w = .34
for i_ in x:                                           # optimism: apparent to cross-validated, with the CV interval
    ax.annotate("", (i_ + w / 2, cv.CV5_mean.iloc[i_]), (i_ - w / 2, cv.Apparent.iloc[i_]),
                arrowprops=dict(arrowstyle="-|>", color=fs.SLATE, lw=.6, mutation_scale=5, shrinkA=3, shrinkB=3), zorder=2)
ax.scatter(x - w / 2, cv.Apparent, s=24, facecolor="white", edgecolor=fs.SLATE, lw=.9, zorder=3)
ax.vlines(x + w / 2, cv.CV5_lo, cv.CV5_hi, color=fs.BLUE, lw=1.2, zorder=3)
ax.scatter(x + w / 2, cv.CV5_mean, s=24, c=fs.BLUE, lw=.4, ec="white", zorder=4)
ax.set_xticks(x); ax.set_xticklabels(["MSI/EBV", "NetImmune", "combined"]); ax.tick_params(axis="x", length=0)
ax.set_ylim(.5, 1.05); ax.set_yticks([.5, .6, .7, .8, .9, 1.0]); ax.set_ylabel("AUC")
ax.scatter([.05], [.95], s=14, facecolor="white", edgecolor=fs.SLATE, lw=.8, transform=ax.transAxes, clip_on=False)
ax.text(.09, .95, "apparent", transform=ax.transAxes, fontsize=K.T_SMALL, va="center")
ax.scatter([.05], [.86], s=14, c=fs.BLUE, lw=0, transform=ax.transAxes, clip_on=False)
ax.text(.09, .86, "five-fold cross-validated (95% CI)", transform=ax.transAxes, fontsize=K.T_SMALL, va="center"); K.frame(ax)
C.letter(96, 52, "d")

# ---- e  TCGA-STAD per signature
PR = pd.read_csv(f"{D}/icb/native_prognosis_matched.csv"); PR["beats"] = PR.C > PR.rand_p95
V3 = pd.read_csv(f"{D}/icb/v3_dissociation.csv").set_index("signature"); famd = V3.family
R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"]
rmed = (R.AUC - R.rand_med).groupby(R.signature).median()
order = []
for f in ["Inflamed/IFN", "Other", "Confounder"]: order += sorted([s for s in rmed.index if famd[s] == f], key=lambda s: -rmed[s])
st = PR[(PR.endpoint == "OS") & (PR.cancer == "STAD")].set_index("signature").reindex(order)
ax = C.ax(14, 124, 164, 26)
for i, (_, r) in enumerate(st.iterrows()):
    ax.plot([i, i], [r.rand_med, r.rand_p95], color=fs.LIGHT, lw=2.2, solid_capstyle="butt", zorder=2)
ax.scatter(range(len(st)), st.C, s=11, c=[fs.FAM[famd[s]] for s in st.index], lw=.3, ec="white", zorder=4)
clin = pd.read_csv(f"{D}/pancan/ceiling2_results.csv").query("endpoint=='OS' and cancer=='STAD'").C_clinical.iloc[0]
ax.axhline(clin, color=fs.TEAL, lw=.9); ax.text(len(st) - .5, clin + .003, "stage + age", color=fs.TEAL, fontsize=K.T_SMALL, ha="right", va="bottom")
ax.axhline(.5, color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(-.7, len(st) - .3); ax.set_ylim(.48, .67); ax.set_yticks([.5, .55, .6, .65])
ax.set_xticks(range(len(st))); ax.set_xticklabels([K.lab(s) for s in st.index], rotation=90, fontsize=5.0)
for t, s in zip(ax.get_xticklabels(), st.index): t.set_color(fs.FAM[famd[s]] if famd[s] != "Other" else fs.SLATE)
ax.tick_params(axis="x", length=0, pad=1.5); ax.set_ylabel("C-index, TCGA-STAD"); K.frame(ax)
C.letter(1.5, 118, "e")
C.save("SF_gastric")
print(f"S-gastric | STAD n {int(st.n.iloc[0])} events {int(st.events.iloc[0])} beat {int(st.beats.sum())} {list(st.index[st.beats == 1])} clin {clin:.3f}")
print(meta.round(3).to_string(index=False)); print(cv.round(3).to_string(index=False))
