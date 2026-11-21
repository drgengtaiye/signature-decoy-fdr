"""Supplementary Figure 7: the breast cancer library in detail (ten decoy sets; panels a and b, set 1).
python run_restyled.py figS7.py

a  share of variance on the first principal component, each signature against the median of its co-expressed and
   random-gene decoys (one point per signature x cohort; applications/v11_breast_pc1_strength.R, post hoc)
b  decoys matched on the mean absolute pairwise correlation (post hoc): their internal |r| against their target's
c  signatures called and the estimated FDR in each analysis: first principal component with decoys matched on the
   signed and on the absolute correlation, and the mean z
d  each signature's median delta C scored by the mean z, against its co-expressed decoys
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v4kit as K
import figstyle as fs
from fig67_lib import bars2, per_signature, C_FOLD, C_DIR, C_SIG, C_DEC
K.apply()
D = K.D + "/icb/"; P = D + "v12_pkg/"; P11 = D + "v11_pkg/"
TT = 6.9
ST = pd.read_csv(P11 + "breast_pc1_strength.csv")
IA = pd.read_csv(P + "breast_abs_decoy_info.csv.gz"); IA = IA[IA.set == 1]
CM = pd.read_csv(P + "breast_mean_calls.csv.gz")
F = {k: pd.read_csv(P + f) for k, f in (("pc1", "breast_pc1_fdr.csv"), ("abs", "breast_abs_pc1_fdr.csv"), ("mean", "breast_mean_fdr.csv"))}
C = K.Canvas(150)

ax = C.ax(15, 15, 45, 40)
ax.scatter(ST.pc1_target, ST.pc1_coex_med, s=3, color=C_SIG, alpha=.55, lw=0, zorder=3, label="co-expressed decoys")
ax.scatter(ST.pc1_target, ST.pc1_rand_med, s=3, color=C_DEC, alpha=.45, lw=0, zorder=2, label="random-gene decoys")
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.5, ls=(0, (3, 2)))
lim = max(ST.pc1_target.max(), ST.pc1_coex_med.max()) * 1.05
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("signature: variance on PC1"); ax.set_ylabel("its decoys (median): variance on PC1")
ax.legend(fontsize=K.T_SMALL, loc="upper left", frameon=False, handletextpad=.2, markerscale=2)
K.frame(ax)
C.title(15, 9, "Matched on mean r: a weaker PC1", size=TT, fontweight="bold"); C.letter(7.5, 8.5, "a")

ax = C.ax(80, 15, 38, 40)
for cls, cc in (("rand", C_DEC), ("coex", C_SIG)):
    g = IA[IA["class"] == cls]
    ax.scatter(g.r_target_med, g.r_achieved_med, s=3, color=cc, alpha=.5, lw=0, zorder=3)
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(0, .6); ax.set_ylim(0, .6)
ax.set_xlabel("signature mean |$r$|"); ax.set_ylabel("decoy mean |$r$|")
K.key(ax, [(C_SIG, "co-expressed"), (C_DEC, "random-gene")], x=.04, y=.97, s=9, size=K.T_SMALL)
K.frame(ax)
C.title(80, 9, "Matched on |r| instead", size=TT, fontweight="bold"); C.letter(72.5, 8.5, "b")

ax = C.ax(135, 15, 35, 40)
def row(df, cls):
    return df[(df.stat == "fold") & (df.decoys == cls) & (df.boot == "clusters")].iloc[0]
combos = [("PC1,\nsigned r", row(F["pc1"], "coex"), C_FOLD), ("PC1,\n|r|", row(F["abs"], "coex"), C_DIR),
          ("mean z", row(F["mean"], "coex"), fs.AMBER)]
bars2(ax, combos, max(c[1].called for c in combos) * 1.45 + 2)
ax.set_ylabel("signatures called (of 45)")
C.title(135, 9, "Calls and the estimate", size=TT, fontweight="bold"); C.letter(127.5, 8.5, "c")

ax = C.ax(40, 72, 50, 74)
per_signature(ax, CM, "fold", fam=None)
ax.set_xlabel("median ΔC over cohorts (folded, mean z)")
C.title(15, 66, "Scored by the mean z, against its decoys", size=TT, fontweight="bold"); C.letter(7.5, 65.5, "d")
C.save("SupplementaryFigure7")
