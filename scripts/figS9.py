"""Additional file 1 Fig. S9 (GB v15): how reproducible a cancer type's value is
(validation/PLAN_round12_reproducibility_by_type.md with amendments 1-3; icb/v18_split from
applications/v18_split_half.R and decoys/v18_reproducibility.py).

a  agreement between random halves: the Pearson correlation over types between two halves' family values, over 500
   splits, by family, TCGA and the independent cohorts; Spearman-Brown reliability at full size
b  each TCGA type's inflamed value with its sampling noise: +-1.96 SD from the halves (pre-specified) and from a
   bootstrap of patients (post hoc); deaths
c  TCGA against the independent cohort of each type, with the same noise on both axes; the discordance test
d  why the halves can move together: the correlation of two halves' values over splits against the whole-cohort value
e  one signature in one cohort: how often a signature above 95% of its decoys in one half is so in the other
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
V = K.D + "/icb/v18_split/"
FAMD = pd.read_csv(K.D + "/icb/v5_signature_clusters.csv").set_index("signature").family
LAB = {"Inflamed/IFN": "inflamed", "Confounder": "confounder", "Other": "other"}
REL = pd.read_csv(V + "reliability.csv").set_index(["setting", "family"])
TV = pd.read_csv(V + "type_values_boot.csv")
TR = pd.read_csv(V + "transport_test.csv").set_index("family")
TRB = pd.read_csv(V + "transport_test_boot.csv").set_index("family")
SL = pd.read_csv(V + "signature_level.csv").set_index("setting")
SET = {"tcga": ("TCGA", K.SLATE, "o", False), "external": ("independent", K.INK, "D", True)}


def halves(setting):
    r = pd.read_csv(V + f"rho_split_{setting}.csv.gz")
    r["family"] = r.target.map(FAMD).map(LAB)
    x = r[r.half != "full"].groupby(["family", "half", "split", "cohort"]).rho.mean()
    return x


H = {s: halves(s) for s in ("tcga", "external")}


def split_r(setting, fam):
    x = H[setting].loc[fam]
    A = x.loc["A"].unstack("cohort"); B = x.loc["B"].unstack("cohort")[A.columns]
    a = A.values - A.values.mean(1, keepdims=True); b = B.values - B.values.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1)), A, B


C = K.Canvas(133)

# ================================================================== a  agreement between halves
ax = C.panel("a", 22, 12, 46, 44, head="Agreement between random halves")
rows = [("tcga", "inflamed"), ("tcga", "confounder"), ("tcga", "other"),
        ("external", "inflamed"), ("external", "confounder"), ("external", "other")]
FC = {"inflamed": K.BLUE, "confounder": K.VERM, "other": K.GRAY}
for i, (s, f) in enumerate(rows):
    y0 = (len(rows) - 1 - i) + (0.45 if s == "tcga" else 0)
    r, _, _ = split_r(s, f)
    K.ridge(ax, y0, r, FC[f], height=0.85, lo=-1, hi=1, alpha=0.75)
    med = np.median(r)
    ax.plot([med, med], [y0, y0 + 0.55], color=K.INK, lw=0.8, zorder=6)
    q = REL.loc[(s, f)]
    ax.text(1.04, y0 + 0.18, f"{q.R_spearman_brown:.2f}", ha="left", va="center", fontsize=K.T_SMALL, color=K.INK)
    ax.text(-1.04, y0 + 0.18, f, ha="right", va="center", fontsize=K.T_SMALL, color=K.INK2)
ax.text(-1.0, 6.25, "TCGA, 20 types", ha="left", va="bottom", fontsize=K.T_NOTE, fontweight="bold", color=K.INK)
ax.text(-1.0, 2.8, "independent, 11 types", ha="left", va="bottom", fontsize=K.T_NOTE, fontweight="bold",
        color=K.INK)
ax.text(1.04, 6.25, "R", ha="left", va="bottom", fontsize=K.T_SMALL, color=K.INK, style="italic")
K.ref(ax, 0, axis="v")
ax.set_xlim(-1, 1); ax.set_ylim(-0.1, 6.9); ax.set_yticks([]); K.frame(ax, left=False)
ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_xticklabels(["−1", "−0.5", "0", "0.5", "1"])
ax.set_xlabel("r between halves over types (500 splits)")

# ================================================================== b  TCGA types with sampling noise
ax = C.panel("b", 92, 12, 34, 104, head="Each TCGA type and its noise")
t = TV[(TV.setting == "tcga") & (TV.family == "inflamed")].sort_values("value")
for j, q in enumerate(t.itertuples()):
    y = j
    ax.plot([q.value - 1.96 * q.sd_boot, q.value + 1.96 * q.sd_boot], [y, y], color=K.tint(K.BLUE, 0.6), lw=2.6,
            solid_capstyle="butt", zorder=2)
    ax.plot([q.value - 1.96 * q.sd_full, q.value + 1.96 * q.sd_full], [y, y], color=K.BLUE, lw=0.9,
            solid_capstyle="butt", zorder=3)
    K.dots(ax, [q.value], [y], K.INK, s=12, z=5)
    ax.text(1.32, y, f"{int(q.events)}", ha="right", va="center", fontsize=K.T_SMALL, color=K.INK2)
ax.text(1.32, len(t) - 0.2, "deaths", ha="right", va="bottom", fontsize=K.T_SMALL, color=K.INK2)
ax.set_yticks(range(len(t))); ax.set_yticklabels(t.cohort, fontsize=K.T_SMALL)
K.ref(ax, 0.5, axis="v")
ax.set_xlim(-0.2, 1.33); ax.set_ylim(-0.7, len(t) - 0.1)
ax.set_xticks([0, 0.5, 1]); ax.set_xlabel("inflamed family, mean ρ")
K.legend(ax, [(K.BLUE, "±1.96 SD, halves", "line"), (K.tint(K.BLUE, 0.6), "±1.96 SD, bootstrap", "line")],
         x=0.06, y=-0.085, dy=0.022, size=K.T_SMALL)

# ================================================================== c  TCGA against independent
ax = C.panel("c", 141, 12, 37, 50, head="TCGA against independent")
x = TV[(TV.setting == "tcga") & (TV.family == "inflamed")].set_index("type")
y = TV[(TV.setting == "external") & (TV.family == "inflamed")].set_index("type")
ty = y.index
K.identity(ax, -0.1, 1.1)
off = TR.loc["inflamed", "offset"]
ax.plot([-0.1, 1.1], [-0.1 + off, 1.1 + off], color=K.RULE, lw=0.5, ls=(0, (3, 2)), zorder=1)
for tp in ty:
    xv, yv = x.loc[tp, "value"], y.loc[tp, "value"]
    sd_x, sd_y = x.loc[tp, "sd_boot"], y.loc[tp, "sd_boot"]
    ax.plot([xv - sd_x, xv + sd_x], [yv, yv], color=K.tint(K.BLUE, 0.45), lw=0.9, solid_capstyle="butt", zorder=2)
    ax.plot([xv, xv], [yv - sd_y, yv + sd_y], color=K.tint(K.BLUE, 0.45), lw=0.9, solid_capstyle="butt", zorder=2)
    K.dots(ax, [xv], [yv], K.INK, s=11, z=5)
LBL = {"SKCM": (-4, 4), "LUSC": (-4, 3), "OV": (4, 4), "STAD": (-4, -4), "CESC": (4, 4), "SARC": (4, -4),
       "LUAD": (4, -4), "HNSC": (-4, -4), "BRCA": (-4, 4), "BLCA": (4, 4), "COAD": (-4, -4)}
for tp in ty:
    dx, dy = LBL.get(tp, (4, 3))
    ax.annotate(tp, (x.loc[tp, "value"], y.loc[tp, "value"]), xytext=(dx, dy), textcoords="offset points",
                fontsize=K.T_SMALL - 0.3, color=K.INK2, ha="left" if dx > 0 else "right", va="center", zorder=6)
ax.set_xlim(-0.05, 1.12); ax.set_ylim(-0.05, 1.12); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
ax.set_xlabel("TCGA, inflamed mean ρ"); ax.set_ylabel("independent cohort")
pfmt = lambda p: "< 0.0001" if p < 1e-4 else f"= {p:.2f}"
ax.text(0.98, 0.5, f"Q = {TR.loc['inflamed', 'Q']:.0f} ({TRB.loc['inflamed', 'Q']:.0f})\n10 df\n"
        f"P {pfmt(TR.loc['inflamed', 'P'])}", transform=ax.transAxes, ha="right", va="center",
        fontsize=K.T_SMALL - 0.2, color=K.INK, linespacing=1.05)

# ================================================================== d  halves moving together
ax = C.panel("d", 22, 76, 46, 40, head="Where halves move together")
for s in ("tcga", "external"):
    lab, col, mk, fill = SET[s]
    x_ = H[s].loc["inflamed"]
    A = x_.loc["A"].unstack("cohort"); B = x_.loc["B"].unstack("cohort")[A.columns]
    cr = pd.Series({c: np.corrcoef(A[c], B[c])[0, 1] for c in A.columns})
    v = TV[(TV.setting == s) & (TV.family == "inflamed")].set_index("cohort").value
    K.dots(ax, v.loc[cr.index], cr, col, s=13, marker=mk, filled=fill, z=5)
K.ref(ax, 0); K.ref(ax, 0.5, axis="v")
ax.set_xlim(0, 1); ax.set_ylim(-1, 1); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticklabels(["−1", "−0.5", "0", "0.5", "1"])
ax.set_xlabel("inflamed family, mean ρ (whole cohort)"); ax.set_ylabel("r of the two halves over splits")
K.legend(ax, [(K.SLATE, "TCGA type", "o", False), (K.INK, "independent cohort", "D", True)], x=0.55, y=0.93,
         dy=0.085, s=8, size=K.T_SMALL)

# ================================================================== e  single-signature claims
ax = C.panel("e", 141, 82, 37, 34, head="One signature, one half")
for i, s in enumerate(("tcga", "external")):
    q = SL.loc[s]
    ax.bar(i, 100 * q.repeat_rate, 0.55, color=K.tint(K.INK, 0.35), edgecolor="white", lw=0.6, zorder=3)
    K.vbar_interval(ax, i, 100 * q.repeat_rate_lo, 100 * q.repeat_rate_hi, K.INK, lw=0.8, z=4)
    ax.plot([i - 0.36, i + 0.36], [100 * q.share_cells_ge95_half] * 2, color=K.VERM, lw=1.1, zorder=5)
    ax.text(i, 100 * q.repeat_rate_hi + 2, f"{100 * q.repeat_rate:.0f}%", ha="center", va="bottom",
            fontsize=K.T_SMALL, color=K.INK)
ax.set_xticks([0, 1]); ax.set_xticklabels(["TCGA", "independent"])
ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 70); ax.set_yticks([0, 20, 40, 60])
ax.set_ylabel("above 95% of decoys in the\nother half too (%)")
K.legend(ax, [(K.VERM, "share of all cells", "line")], x=0.08, y=0.95, size=K.T_SMALL)

C.save("FigS9_reproducibility_GB15")
print("FigS9_reproducibility_GB15")
