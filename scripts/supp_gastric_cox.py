"""Supplementary figure - the gastric cohorts at closer inspection (redrawn).

a  TCGA-STAD overall survival: hazard ratios for NetImmune, BarrierIndex, age and stage,
   each alone (open) and all four together (filled)
b  the genes of the six modules across all 383 TCGA-STAD tumors, ordered by NetImmune
(a uses the 314 patients with stage and age recorded; b needs neither)
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from lifelines import CoxPHFitter
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D
d = pd.read_csv(f"{D}/revision/results/B2_scores_TCGA.csv")
d = d[d.time.notna() & (d.time > 0) & d.event.notna() & d.NetImmune.notna() & d.age.notna() & d.stage_ord.notna()].copy()
for c in ["NetImmune", "Barrier"]: d[c + "_z"] = (d[c] - d[c].mean()) / d[c].std()
VARS = [("NetImmune_z", "NetImmune (per s.d.)"), ("Barrier_z", "BarrierIndex (per s.d.)"),
        ("age", "age (per year)"), ("stage_ord", "stage (per step)")]
uni = {v: CoxPHFitter().fit(d[["time", "event", v]], "time", "event").summary.loc[v] for v, _ in VARS}
mv = CoxPHFitter().fit(d[["time", "event"] + [v for v, _ in VARS]], "time", "event").summary
C = K.Canvas(162)
ax = C.ax(40, 10, 50, 26)
for i, (v, lab_) in enumerate(VARS):
    for k, (s, mk, fc, off) in enumerate([(uni[v], "o", "white", -.14), (mv.loc[v], "o", None, .14)]):
        hr, lo, hi = np.exp(s["coef"]), np.exp(s["coef lower 95%"]), np.exp(s["coef upper 95%"])
        c = fs.BLUE if v == "NetImmune_z" else (fs.VERM if v == "Barrier_z" else fs.SLATE)
        ax.plot([lo, hi], [i + off] * 2, color=c, lw=.9)
        ax.scatter(hr, i + off, s=16, facecolor=fc or c, edgecolor=c, lw=.8, zorder=4)
        ax.text(1.02, i + off, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", transform=ax.get_yaxis_transform(),
                fontsize=K.T_SMALL, va="center", color=fs.INK if k else fs.SLATE)
ax.axvline(1, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_xscale("log"); ax.minorticks_off(); ax.set_xticks([.6, .8, 1, 1.5, 2]); ax.set_xticklabels(["0.6", "0.8", "1", "1.5", "2"])
ax.set_xlim(.55, 2.3); ax.set_yticks(range(4)); ax.set_yticklabels([l for _, l in VARS], fontsize=K.T_SMALL)
ax.set_ylim(3.5, -.5); ax.tick_params(axis="y", length=0)
ax.set_xlabel("hazard ratio, overall survival (95% CI)")
K.frame(ax, left=False)
K.key(ax, [(fs.SLATE, "alone")], x=1.02, y=1.16, s=12)
ax.scatter([1.02], [1.16], s=12, facecolor="white", edgecolor=fs.SLATE, lw=.8, transform=ax.transAxes, clip_on=False, zorder=11)
K.key(ax, [(fs.SLATE, "all four together")], x=1.32, y=1.16, s=12)
C.letter(1.5, 4, "a")

G = pd.read_csv(f"{D}/icb/v4_tcga_stad_modgenes.csv", index_col=0)
dall = pd.read_csv(f"{D}/revision/results/B2_scores_TCGA.csv").set_index("sample")
G = G.loc[G.index.intersection(dall.index)]
Z = (G - G.mean()) / G.std()
order = dall.loc[Z.index].NetImmune.sort_values().index
MODS = [("APM", ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"]),
        ("IFN", ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"]),
        ("TLS/B", ["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"]),
        ("T cell", ["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"]),
        ("CAF/ECM", ["COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1"]),
        ("TAM", ["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"])]
rows, bounds = [], []
for m, gs_ in MODS:
    gs_ = [g for g in gs_ if g in Z.columns]; bounds.append((m, len(rows), len(gs_))); rows += gs_
M = Z.loc[order, rows].T.values
ax = C.ax(30, 56, 140, 98)
ax.imshow(np.clip(M, -2.5, 2.5), aspect="auto", cmap=fs.DIV, norm=Normalize(-2.5, 2.5), interpolation="nearest")
ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=5.0); ax.set_xticks([])
ax.tick_params(length=0, pad=1)
for m, s0, n in bounds[1:]: ax.axhline(s0 - .5, color="white", lw=1.2)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel(f"TCGA-STAD tumors (n = {len(order)}), ordered by NetImmune from low to high", fontsize=K.T_SMALL)
for m, s0, n in bounds:
    c = fs.BLUE if m in ("APM", "IFN", "TLS/B", "T cell") else fs.VERM
    C.text(3, 56 + 98 * (s0 + n / 2) / len(rows), m, size=K.T_SMALL, va="center", color=c, fontweight="bold")
tr = C.ax(30, 50.5, 140, 4)
nv = dall.loc[order].NetImmune.values
tr.bar(range(len(nv)), nv, width=1, color=np.where(nv > 0, fs.BLUE, fs.GREY), lw=0)
tr.set_xlim(-.5, len(nv) - .5); fs.clean(tr)
C.text(28.5, 52.5, "NetImmune", size=K.T_SMALL, ha="right", va="center")
cax = C.ax(172, 86, 1.3, 36)
cb = C.fig.colorbar(plt.cm.ScalarMappable(Normalize(-2.5, 2.5), fs.DIV), cax=cax, ticks=[-2, 0, 2])
cb.outline.set_linewidth(.3); cax.tick_params(labelsize=5.0, length=1.3, pad=1); cax.set_title("z", fontsize=K.T_SMALL, pad=2)
C.letter(1.5, 46, "b")
C.save("SF_gcox")
print(f"S-gastric-cox | n {len(d)} events {int(d.event.sum())}")
for v, lab_ in VARS:
    u = uni[v]; m_ = mv.loc[v]
    print(f"  {lab_}: uni {np.exp(u['coef']):.2f} ({np.exp(u['coef lower 95%']):.2f}-{np.exp(u['coef upper 95%']):.2f}) P {u['p']:.2g} | "
          f"multi {np.exp(m_['coef']):.2f} ({np.exp(m_['coef lower 95%']):.2f}-{np.exp(m_['coef upper 95%']):.2f}) P {m_['p']:.2g}")
print(f"  heatmap genes {len(rows)} samples {len(order)}")
