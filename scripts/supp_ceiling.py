"""Supplementary figure - the information ceiling.

a  per cancer type, C-index of a cross-validated full-transcriptome Cox model, of the 48-gene
   NetImmune score and of matched random gene sets
b  the medians across 20 cancer types
c  in the gastric checkpoint cohort, AUC against the number of genes used
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
c = pd.read_csv(f"{K.D}/pancan/ceiling_results.csv").dropna(subset=["C_full_transcriptome"]).sort_values("C_full_transcriptome")
C = K.Canvas(62)
ax = C.ax(14, 12, 86, 32); x = np.arange(len(c))
for i in x: ax.plot([i, i], [c.C_NetImmune.iloc[i], c.C_full_transcriptome.iloc[i]], color=fs.LIGHT, lw=.8, zorder=1)
ax.scatter(x, c.C_full_transcriptome, s=14, marker="s", c=fs.SLATE, lw=0, zorder=3)
ax.scatter(x, c.C_NetImmune, s=13, c=fs.BLUE, lw=0, zorder=3)
ax.scatter(x, c.C_rand_dir_med, s=26, marker="_", c=fs.VERM, lw=1.0, zorder=2)
ax.axhline(.5, color=fs.INK, ls=(0, (3, 2)), lw=.5)
ax.set_xticks(x); ax.set_xticklabels(c.cancer, rotation=90, fontsize=5.0); ax.tick_params(axis="x", length=0)
ax.set_xlim(-.7, len(c) - .3); ax.set_ylabel("C-index (OS)"); K.frame(ax)
K.key(ax, [(fs.SLATE, "full transcriptome (cross-validated)")], x=.02, y=1.10, marker="s", s=12)
K.key(ax, [(fs.BLUE, "48-gene NetImmune")], x=.45, y=1.10, s=12)
K.key(ax, [(fs.VERM, "random gene sets")], x=.72, y=1.10, marker="s", s=10)
C.letter(1.5, 3, "a")
ax = C.ax(116, 12, 24, 32)
vals = [c.C_full_transcriptome.median(), c.C_NetImmune.median(), c.C_rand_dir_med.median()]
V3c = c[["C_full_transcriptome", "C_NetImmune", "C_rand_dir_med"]].to_numpy()
for r_ in V3c: ax.plot(range(3), r_, color=fs.LIGHT, lw=.5, zorder=1)          # one line per cancer type
for k_, col_ in enumerate([fs.SLATE, fs.BLUE, fs.VERM]):
    ax.scatter(np.full(len(V3c), k_), V3c[:, k_], s=5, c=col_, lw=0, alpha=.7, zorder=2)
    ax.plot([k_ - .22, k_ + .22], [vals[k_]] * 2, color=col_, lw=1.6, zorder=3)
ax.set_xticks(range(3)); ax.set_xticklabels(["all\ngenes", "48\ngenes", "random"], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_xlim(-.45, 2.45); ax.set_ylim(.4, .88); ax.set_yticks([.4, .5, .6, .7, .8]); ax.axhline(.5, color=fs.INK, ls=(0, (3, 2)), lw=.5)
ax.set_ylabel("C-index (OS)"); K.frame(ax)
C.letter(104, 3, "b")
ax = C.ax(154, 12, 24, 32)
av = [.754, .797, .826, .854]
ax.plot(range(4), av, color=fs.BLUE, lw=.9, zorder=2)
ax.scatter(range(4), av, s=22, c=[fs.SLATE, fs.BLUE, fs.BLUE, fs.BLUE], lw=.5, ec="white", zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels(["all", "48", "24", "5"], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_xlabel("genes used", fontsize=K.T_SMALL); ax.set_ylim(.5, .9); ax.set_yticks([.5, .6, .7, .8, .9]); ax.set_ylabel("AUC, pembrolizumab"); K.frame(ax)
C.letter(142, 3, "c")
C.save("SF_ceiling")
print("S-ceiling | medians", [round(v, 3) for v in vals], "| n", len(c))
