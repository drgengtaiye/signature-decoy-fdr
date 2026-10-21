"""Additional file 1, Fig. S5 (Genome Biology) - the third randomized trial, I-SPY2.

Panels g-i of the v6 Figure 8, moved here when the locked-prediction validation went to a
separate paper; panel b is coloured by signature family instead of by locked triage call.

a  I-SPY2: score x arm interaction on pCR for the 17 immune and 3 control gene sets
b  I-SPY2 pembrolizumab contrast: each signature's AUC in the two arms
c  every randomized comparison: immune gene sets favoring the checkpoint arm, and hits
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"); sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v4kit as K
import figstyle as fs
from ispy2_interaction import logit_fit
from v7_lib import signatures, IMM
K.apply()
D = K.D; V = f"{D}/validation"; I = f"{D}/icb"
FAM = pd.read_csv(f"{I}/v3_dissociation.csv").set_index("signature").family
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
# trial gene-set names exactly as in Fig. 5 (fig4_v6.py)
SHORT = {"GEP18_TIS (Ayers 2017)": "GEP-18/TIS (Ayers)", "IFNG_6gene (Ayers 2017)": "IFNG 6-gene (Ayers)",
         "Chemokine12 (Messina 2012)": "Chemokine-12 (Messina)", "TLS (Cabrita 2020)": "TLS (Cabrita)",
         "TLS_9gene (Prabhakaran)": "TLS 9-gene (Prabhakaran)", "CTRL_Hypoxia (Buffa)": "Hypoxia (Buffa)",
         "CTRL_Proliferation": "Proliferation", "CTRL_Ribosome": "Ribosome", "barriernet_Immune": "NetImmune",
         "CD8_Tcells": "CD8 T cells", "T_cells": "T cells", "NK_cells": "NK cells", "B_cells": "B cells", "Plasma_cells": "Plasma cells"}
short = lambda s: SHORT.get(s, s.replace("_", " "))
C = K.Canvas(92)

# ============================================================ a  I-SPY2 interaction, trial library
def contrast(nm):
    X = pd.read_csv(f"{V}/cohorts/{nm}.expr.tsv.gz", sep="\t", index_col=0); X.columns = [str(c) for c in X.columns]
    Cc = pd.read_csv(f"{V}/cohorts/{nm}.clin.tsv", sep="\t"); Cc["sample"] = Cc["sample"].astype(str)
    return X[Cc["sample"]], Cc
SIGS = signatures(); IS = {}
for nm in ["ISPY2_pembro", "ISPY2_durva"]:
    X, Cc = contrast(nm); Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T; arm, y = Cc.arm.to_numpy(float), Cc.y.to_numpy(float)
    rows = []
    for s, gl in SIGS.items():
        gg = [x for x in gl if x in Z.index]
        fam = "Immune" if s in IMM else "Control" if s.startswith("CTRL_") else "Other"
        if len(gg) < 3 or fam == "Other": continue
        v = Z.loc[gg].mean(0).to_numpy(); v = (v - v.mean()) / v.std()
        b, se, _ = logit_fit(np.column_stack([np.ones_like(v), v, arm, v * arm]), y)
        rows.append(dict(set=s, fam=fam, beta=b[3], se=se[3]))
    IS[nm] = pd.DataFrame(rows).set_index("set")
order_g = IS["ISPY2_pembro"].sort_values(["fam", "beta"], ascending=[False, False]).index.tolist()
ax = C.ax(30, 12, 34, 70)
for k_, (nm, mk, lab_) in enumerate([("ISPY2_pembro", "o", "pembrolizumab"), ("ISPY2_durva", "D", "durvalumab + olaparib")]):
    t = IS[nm].reindex(order_g); y_ = np.arange(len(t)) + (k_ - .5) * .32
    cols = [fs.BLUE if f_ == "Immune" else fs.VERM for f_ in t.fam]
    ax.hlines(y_, t.beta - 1.96 * t.se, t.beta + 1.96 * t.se, color=cols, lw=.7, alpha=.9 if k_ == 0 else .5)
    ax.scatter(t.beta, y_, s=11, marker=mk, c=cols if k_ == 0 else "white", ec=cols, lw=.7, zorder=3)
ax.axvline(0, color=fs.INK, lw=.5)
ax.set_yticks(range(len(order_g))); ax.set_yticklabels([short(s) for s in order_g], fontsize=5.0)
ax.tick_params(axis="y", length=0); ax.set_ylim(len(order_g) - .4, -.8); ax.set_xlim(-2.2, 2.6)
ax.set_xlabel("interaction, log odds of pCR per s.d.\n(positive favors the checkpoint arm)", fontsize=K.T_SMALL); K.frame(ax, left=False)
ax.scatter([.02], [1.06], s=10, c=fs.INK, transform=ax.transAxes, clip_on=False); ax.text(.07, 1.06, "pembrolizumab", transform=ax.transAxes, fontsize=5.0, va="center")
ax.scatter([.02], [1.025], s=10, marker="D", c="white", ec=fs.INK, lw=.7, transform=ax.transAxes, clip_on=False); ax.text(.07, 1.025, "durvalumab + olaparib", transform=ax.transAxes, fontsize=5.0, va="center")
C.title(2, 4, "Interaction with arm", fontweight="bold")
C.letter(1.5, 1, "a")

# ============================================================ b  I-SPY2 pembrolizumab: AUC by arm
N = pd.read_csv(f"{V}/cohorts/ISPY2_pembro.native_scores.csv", index_col=0); N.index = N.index.astype(str)
_, Cc = contrast("ISPY2_pembro"); N = N.loc[Cc["sample"]]
A = []
for s in [x for x in FAM.index if x in N.columns]:
    v = N[s].to_numpy()
    A.append(dict(sig=s, fam=FAM[s], icb=roc_auc_score(Cc.y[Cc.arm == 1], v[Cc.arm == 1]), ctl=roc_auc_score(Cc.y[Cc.arm == 0], v[Cc.arm == 0])))
A = pd.DataFrame(A)
ax = C.ax(84, 16, 44, 44)
for f_ in ["Other", "Confounder", "Inflamed/IFN"]:
    h = A[A.fam == f_]; ax.scatter(h.ctl, h.icb, s=14, c=fs.FAM[f_], lw=.35, ec="white", zorder=3)
ax.plot([.3, .85], [.3, .85], color=fs.INK, lw=.5, ls=":")
mt = A.groupby("fam")[["icb", "ctl"]].median()
ax.text(.03, .97, f"inflamed/IFN signatures, median AUC\n{mt.loc['Inflamed/IFN', 'icb']:.2f} + pembrolizumab\n{mt.loc['Inflamed/IFN', 'ctl']:.2f} paclitaxel",
        transform=ax.transAxes, va="top", fontsize=5.0)
K.key(ax, [(fs.FAM[f_], fs.FAMLAB[f_]) for f_ in FAMORD], .60, .22, dy=.075, size=5.0, s=12, xpad=.035)
ax.set_xlim(.32, .85); ax.set_ylim(.32, .85); ax.set_xticks([.4, .6, .8]); ax.set_yticks([.4, .6, .8])
ax.set_xlabel("AUC for pCR, paclitaxel", fontsize=K.T_SMALL); ax.set_ylabel("AUC, paclitaxel + pembrolizumab", fontsize=K.T_SMALL); K.frame(ax)
C.title(74, 4, "Pembrolizumab contrast, both arms", fontweight="bold")
C.letter(73.5, 1, "b")

# ============================================================ c  all randomized comparisons
T7 = pd.read_csv(f"{D}/pancan/v7_persignature_exact.csv"); T7 = T7[T7.family == "Immune"]
cmp_ = []
for tr, lab_ in [("IMmotion150", "IMmotion150, renal"), ("PRINCE", "PRINCE, pancreatic")]:
    t = T7[T7.trial == tr]; cmp_.append((lab_, int((t.beta < 0).sum()), int(((t.wald_p < .05) & (t.beta < 0)).sum()), len(t)))
for nm, lab_ in [("ISPY2_pembro", "I-SPY2, pembrolizumab"), ("ISPY2_durva", "I-SPY2, durvalumab + olaparib")]:
    t = IS[nm][IS[nm].fam == "Immune"]; p_ = 2 * norm.sf(np.abs(t.beta / t.se))
    cmp_.append((lab_, int((t.beta > 0).sum()), int(((p_ < .05) & (t.beta > 0)).sum()), len(t)))
ax = C.ax(138, 16, 38, 40)
for i, (lab_, fav, hit, n_) in enumerate(cmp_):
    ax.barh([i], [n_], height=.42, color=fs.LIGHT, lw=0, zorder=1)
    ax.barh([i], [fav], height=.42, color=fs.BLUE, alpha=.55, lw=0, zorder=2)
    ax.scatter([hit], [i], marker="D", s=16, c=fs.INK, zorder=3)
    ax.text(n_ + .6, i, f"{fav}/{n_}; {hit}", va="center", fontsize=5.0)
    ax.text(0, i - .3, lab_, va="bottom", ha="left", fontsize=5.0)          # name above its bar, inside the panel
ax.set_yticks([]); ax.set_ylim(len(cmp_) - .5, -.95); ax.set_xlim(-.5, 23); ax.set_xticks([0, 5, 10, 15])
ax.set_xlabel("immune gene sets (of 17)", fontsize=K.T_SMALL); K.frame(ax, left=False)
for yy, mk, c, lab_ in [(-.30, "s", fs.BLUE, "favor the checkpoint arm"), (-.37, "D", fs.INK, "and P < 0.05")]:
    ax.scatter([.03], [yy], s=12, marker=mk, c=c, alpha=.55 if mk == "s" else 1, lw=0, transform=ax.transAxes, clip_on=False)
    ax.text(.08, yy, lab_, transform=ax.transAxes, va="center", fontsize=5.0)
C.title(132, 4, "All four randomized comparisons", fontweight="bold")
C.letter(131.5, 1, "c")

C.save("FigS_ISPY2")
print("pembro interactions >0:", int((IS["ISPY2_pembro"].query("fam == 'Immune'").beta > 0).sum()), "| comparisons", cmp_)
print(mt.round(3)); print("AUC signatures:", len(A))
