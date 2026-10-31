"""Additional file 1: Fig. S4 (GB v7) - prognosis in 20 TCGA cancer types, in detail (ten decoy sets; panel a, set 1,
which is the decoy set of applications/v11_decoy_r_by_cohort.R).

a  decoys' internal correlation against their target's in each cancer type (median over a target's co-expressed
   decoys; one point per signature x cancer type)
b  under 200 within-type permutations of (time, event): the share of units called, by kind and statistic (nominal 2.5%)
c  each signature's median delta C, directed statistic (overall survival), against its co-expressed decoys
d  signatures called by family, for each endpoint and statistic
"""
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
from fig67_lib import per_signature, C_FOLD, C_DIR, C_SIG, C_DEC
K.apply()
D = K.D + "/icb/"; P = D + "v12_pkg/"; P11 = D + "v11_pkg/"
TT = 6.9
FAM = pd.read_csv(D + "v5_signature_clusters.csv").set_index("signature").family.to_dict()
CL = pd.read_csv(P + "tcga_OS_calls.csv.gz")
GN = pd.read_csv(P + "tcga_OS_globalnull_summary.csv").set_index("stat")
RP = pd.read_csv(P11 + "tcga_decoy_r_by_type.csv")          # written by applications/v11_decoy_r_by_cohort.R
C = K.Canvas(150)

ax = C.ax(15, 15, 45, 40)
ax.scatter(RP.r_target, RP.r_decoy_med, s=3, color=C_SIG, alpha=.5, lw=0, zorder=3)
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(-.05, 1); ax.set_ylim(-.05, 1)
ax.set_xlabel("signature internal $r$"); ax.set_ylabel("its co-expressed decoys' internal $r$\n(median)")
fs.note(ax, .04, .95, "one point per signature\n× cancer type", size=K.T_SMALL, color=fs.SLATE)
K.frame(ax)
C.title(15, 9, "Decoys fall short of coherent signatures", size=TT, fontweight="bold"); C.letter(7.5, 8.5, "a")

ax = C.ax(80, 15, 85, 40)
groups = [("signatures", "rate_targets"), ("co-expressed\ndecoys", "rate_coex"), ("random-gene\ndecoys", "rate_rand")]
for gi, (lab, col) in enumerate(groups):
    for si, (st, cc) in enumerate([("fold", C_FOLD), ("dir", C_DIR)]):
        v = GN.loc[st, col]; x = gi + (si - .5) * .36
        ax.bar(x, v, .34, color=cc, lw=0, zorder=3)
        ax.text(x, v + .08, f"{v:.1f}", ha="center", va="bottom", fontsize=K.T_SMALL, color=cc, fontweight="bold")
ax.axhline(2.5, color=fs.INK, lw=.55, ls=(0, (3, 2))); ax.text(2.6, 2.62, "2.5%", fontsize=K.T_SMALL, va="bottom", ha="right")
ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylim(0, max(GN.rate_targets.max(), GN.rate_coex.max()) * 1.3); ax.set_xlim(-.6, 2.7)
ax.set_ylabel("called per unit (%)")
K.key(ax, [(C_FOLD, "folded"), (C_DIR, "directed")], x=.78, y=.96, s=9, size=K.T_SMALL, marker="s")
K.frame(ax)
C.title(80, 9, "Permuted survival: the folded call runs above nominal for every unit", size=TT, fontweight="bold")
C.letter(72.5, 8.5, "b")

ax = C.ax(40, 72, 50, 74)
per_signature(ax, CL, "dir", fam=FAM)
ax.set_xlabel("median ΔC over 20 cancer types (directed, OS)")
C.title(15, 66, "Each signature against its decoys, directed", size=TT, fontweight="bold"); C.letter(7.5, 65.5, "c")

ax = C.ax(115, 74, 50, 40)
rows = []
for ep in ("OS", "DSS", "PFI"):
    f = P + f"tcga_{ep}_calls.csv.gz"
    if not os.path.exists(f): continue
    c = pd.read_csv(f); c = c[(c.kind == "target") & c.call]
    for st in ("fold", "dir"):
        cc = c[c.stat == st]
        rows.append((f"{ep}\n{'folded' if st == 'fold' else 'directed'}",
                     [int((cc.target.map(FAM) == f).sum()) for f in ("Inflamed/IFN", "Other", "Confounder")]))
for k, (lab, v) in enumerate(rows):
    b = 0
    for f, n in zip(("Inflamed/IFN", "Other", "Confounder"), v):
        ax.bar(k, n, .7, bottom=b, color=fs.FAM[f], lw=0, zorder=3); b += n
ax.set_xticks(range(len(rows))); ax.set_xticklabels([r[0] for r in rows], fontsize=4.8); ax.tick_params(axis="x", length=0)
ax.set_ylabel("signatures called")
K.key(ax, [(fs.FAM[f], fs.FAMLAB[f]) for f in ("Inflamed/IFN", "Other", "Confounder")], x=.02, y=.97, s=9, size=K.T_SMALL, marker="s")
K.frame(ax)
C.title(115, 66, "What is called, by family", size=TT, fontweight="bold"); C.letter(107.5, 65.5, "d")
C.save("FigS4_tcga_GB7")
