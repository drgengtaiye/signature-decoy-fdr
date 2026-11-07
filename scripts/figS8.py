"""Additional file 1: Fig. S6 (GB v7) - the decoy pool at library scale (analysis 8 of
validation/PLAN_round7_estimand_and_bias.md), on the 20 TCGA cancer types.

a  genes left for decoys (measured in all 20 types and in no set) as the library grows: the 52 immune signatures, MSigDB
   Hallmark, Hallmark + C2 canonical pathways, Hallmark + C2 CP + C5 GO biological process
b  for each of the 50 Hallmark sets, the internal correlation its decoys reach (median of 4) against its own, from each
   MSigDB pool
c  the Hallmark library on overall survival (mean z, folded C-index, ten decoy sets): calls and the estimated FDR
   (plain, and pi0-corrected beneath), with the global null
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
from fig67_lib import bars2, C_FOLD, C_SIG, C_DEC
K.apply()
P = K.D + "/icb/v12_pkg/"
TT = 6.9
PS = pd.read_csv(P + "hallmark_pool_sizes.csv")
PD = pd.read_csv(P + "hallmark_pool_decoys.csv.gz")
F = pd.read_csv(P + "hallmark_OS_fdr.csv")
GN = pd.read_csv(P + "hallmark_OS_globalnull_summary.csv")
LAB = {"immune": "52 immune\nsignatures", "hallmark": "Hallmark", "hallmark_cp": "+ C2\npathways",
       "hallmark_cp_gobp": "+ GO\nprocesses"}
PCOL = {"hallmark": C_SIG, "hallmark_cp": fs.AMBER, "hallmark_cp_gobp": fs.VERM}
C = K.Canvas(62)

ax = C.ax(15, 13, 42, 38)
for i, r in enumerate(PS.itertuples()):
    ax.bar(i, r.pool / 1000, .65, color=PCOL.get(r.library, C_DEC), lw=0, zorder=3)
    ax.text(i, r.pool / 1000 + .3, f"{r.pool:,}\n({r.sets:,} sets)", ha="center", va="bottom", fontsize=K.T_SMALL - .3)
ax.set_xticks(range(len(PS))); ax.set_xticklabels([LAB[x] for x in PS.library], fontsize=K.T_SMALL - .3)
ax.tick_params(axis="x", length=0); ax.set_ylabel("genes left for decoys (thousands)"); ax.set_ylim(0, 23)
fs.note(ax, .98, .97, f"of {PS.universe.iloc[0]:,} genes measured\nin all 20 types", size=K.T_SMALL, color=fs.SLATE, ha="right")
K.frame(ax)
C.title(15, 7, "The pool shrinks as the library grows", size=TT, fontweight="bold"); C.letter(7.5, 6.5, "a")

ax = C.ax(75, 13, 38, 38)
for pool, g in PD.groupby("pool", sort=False):
    m = g.groupby("target").agg(t=("r_target_med", "first"), a=("r_achieved_med", "median"))
    ax.scatter(m.t, m.a, s=7, color=PCOL[pool], lw=0, alpha=.85, zorder=3,
               label={"hallmark": "Hallmark pool", "hallmark_cp": "+ C2 canonical pathways", "hallmark_cp_gobp": "+ GO biological processes"}[pool])
lim = max(PD.r_target_med.max(), PD.r_achieved_med.max()) * 1.08
ax.plot([0, lim], [0, lim], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("Hallmark set internal $r$"); ax.set_ylabel("its decoys' internal $r$ (median of 4)")
ax.legend(fontsize=K.T_SMALL - .3, loc="upper left", frameon=False, handletextpad=.2, markerscale=1.4)
K.frame(ax)
C.title(75, 7, "Decoys still reach the sets' coherence", size=TT, fontweight="bold"); C.letter(67.5, 6.5, "b")

ax = C.ax(130, 13, 26, 38)
row = lambda d, b="clusters": F[(F["class"] == d) & (F.boot == b)].iloc[0]
combos = [("co-expressed", row("coex"), C_FOLD), ("random-gene", row("rand"), C_DEC)]
bars2(ax, combos, max(c[1].called for c in combos) * 1.45 + 2)
ax.set_ylabel("Hallmark sets called (of 50)")
g = GN[GN.stat == "fold"].iloc[0]
fs.note(ax, .03, .97, f"global null: ratio\n{g.ratio:.2f} ({g.ratio_lo:.2f}–{g.ratio_hi:.2f})", size=K.T_SMALL, color=fs.SLATE)
C.title(130, 7, "Hallmark, overall survival", size=TT, fontweight="bold"); C.letter(122.5, 6.5, "c")
C.save("FigS6_scale_GB7")
