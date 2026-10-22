"""Figure 7 (GB v10) - a library outside immunology: the breast-cancer prognostic signatures of Venet et al.
(folded C-index, first principal component unless marked, ten decoy sets).

a  the matched null is itself prognostic: per cohort, the median folded C of each signature's null bank
b  the 13 cohorts: patients and events
c  signatures called and the estimate (pi0-corrected open) with intervals over clusters, across analyses
d  leaving each signature's development cohorts out: median delta before and after, for the 11 signatures affected
e  what a first-principal-component score uses: variance on PC1, signatures against their decoys
f  under 200 within-cohort permutations of the outcome: realized false calls and the estimate
g  one decoy set against ten
h  each signature's median delta against its 320 co-expressed decoys
i  provenance of each signature in each cohort, compiled from the original publications before it was applied
j  each signature's gain over its matched null in each cohort (same columns as h and i)
k  Kaplan-Meier curves in METABRIC split at the median of SOTIRIOU-GGI and of its median co-expressed decoy (an
   illustration chosen by effect size; validation/PLAN_round8_figure_displays.md)
l  the same signature's gain in every cohort against its co-expressed decoys' (first decoy set)
Data: icb/v12_pkg/breast_* (v12_rerun.R, v12_breast_origin.R), v11_pkg/breast_pc1_strength.csv (set 1),
breast/cohorts/manifest.csv, validation/venet_provenance*.csv.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from matplotlib.patches import Rectangle
import cns10 as K
K.apply()
D = K.D + "/icb/"; P = D + "v12_pkg/"; VAL = K.D + "/signature-decoy-fdr/validation/"
F = {k: pd.read_csv(P + f"breast_{k}_fdr.csv") for k in ("pc1", "mean", "abs_pc1")}
OR = pd.read_csv(P + "breast_origin_fdr.csv")
B = {k: pd.read_csv(P + f"breast_{k}_fdr_by_set.csv") for k in ("pc1", "mean", "abs_pc1")}
CL = pd.read_csv(P + "breast_pc1_calls.csv.gz"); CO = pd.read_csv(P + "breast_pc1_origin_calls.csv.gz")
U = pd.read_csv(P + "breast_pc1_units.csv.gz"); U = U[U.kind == "target"]
GN = pd.read_csv(P + "breast_pc1_globalnull_summary.csv"); GP = pd.read_csv(P + "breast_pc1_globalnull_per_perm.csv")
ST = pd.read_csv(D + "v11_pkg/breast_pc1_strength.csv")
MAN = pd.read_csv(K.D + "/breast/cohorts/manifest.csv"); MAN = MAN[MAN.included.astype(str) == "True"].set_index("cohort")
PV = pd.read_csv(VAL + "venet_provenance.csv").set_index("signature"); BL = pd.read_csv(VAL + "venet_provenance_borderline.csv")
EXC = pd.read_csv(P + "breast_origin_excluded_tests.csv")
COH = ["NKI", "Loi_GSE6532", "TCGA_BRCA", "GSE7390", "GSE11121", "GSE1456", "GSE4922", "GSE20685", "GSE21653", "GSE42568",
       "GSE17705", "GSE96058", "METABRIC"]
CLAB = {"Loi_GSE6532": "Loi", "TCGA_BRCA": "TCGA", "GSE96058": "SCAN-B"}
cl_ = lambda c: CLAB.get(c, c)


def row(df, decoys="coex", boot="clusters"):
    return df[(df.stat == "fold") & (df.decoys == decoys) & (df.boot == boot)].iloc[0]


def orow(score, rule, boot="clusters"):
    return OR[(OR.score == score) & (OR["rule"] == rule) & (OR.boot == boot)].iloc[0]


C = K.Canvas(231)

# ================================================================== a  the null is prognostic
ax = C.panel("a", 14, 12, 44, 28, head="Random gene sets are prognostic")
order = U.groupby("cohort").fold_med.median().sort_values().index
for i, co in enumerate(order):
    v = U[U.cohort == co].fold_med
    ax.scatter(i - .12 + K.swarm(v.to_numpy(), .22, 14), v, s=2.0, color=K.tint(K.SLATE, .35), lw=0, zorder=2)
    K.band(ax, i + .25, v=v, color=K.INK, width=.55, dot_s=6, z=3)
K.ref(ax, 0.5)
ax.set_xticks(range(len(order))); ax.set_xticklabels([cl_(c) for c in order], rotation=60, ha="right", fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0, pad=0.5); ax.set_ylabel("null median folded C")

# ================================================================== b  cohorts
ax = C.panel("b", 78, 12, 22, 28, head="13 cohorts")
mo = MAN.loc[COH].sort_values("n_samples")
y = np.arange(len(mo))
ax.barh(y, mo.n_samples, 0.62, color=K.LGRAY, edgecolor="white", linewidth=0.4, label="patients")
ax.barh(y, mo.n_events, 0.62, color=K.SLATE, edgecolor="white", linewidth=0.4, label="events")
ax.set_xscale("log"); ax.set_xlim(20, 5000); ax.set_xticks([30, 300, 3000]); ax.set_xticklabels(["30", "300", "3,000"]); ax.minorticks_off()
ax.set_yticks(y); ax.set_yticklabels([f"{cl_(c)} ({mo.loc[c, 'endpoint']})" for c in mo.index], fontsize=4.6)
K.frame(ax, left=False); ax.set_ylim(-0.6, len(mo) - 0.4)
K.legend(ax, [(K.LGRAY, "patients", "box"), (K.SLATE, "events", "box")], x=0.55, y=0.18, dy=0.1, s=7, size=K.T_SMALL)

# ================================================================== c  analyses
ax = C.panel("c", 136, 12, 42, 28, head="Calls and the estimate")
rows = [("PC1, every cohort", row(F["pc1"]), K.SLATE), ("PC1, development out", orow("pc1", "D and C"), K.SLATE),
        ("PC1, + borderline out", orow("pc1", "D, C and every borderline"), K.SLATE),
        ("PC1, closest decoys", row(F["pc1"], "coex closest"), K.SLATE), ("PC1, abs. r decoys (post hoc)", row(F["abs_pc1"]), K.AQUA),
        ("mean z", row(F["mean"]), K.ORANGE)]
for k, (lab, r, cc) in enumerate(rows):
    K.hbar_interval(ax, k, r.fdr_lo, min(r.fdr_hi, 70), K.tint(cc, .45), lw=1.4)
    K.dots(ax, [r.fdr], [k], cc, s=11, z=5); K.dots(ax, [r.fdr_pi0], [k], cc, s=11, filled=False, z=5)
    ax.text(71, k, f"{int(r.called)}", va="center", fontsize=K.T_SMALL, color=K.INK2)
ax.text(71, -0.9, "called", fontsize=K.T_SMALL, color=K.INK2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=K.T_SMALL)
ax.set_ylim(len(rows) - 0.4, -0.6); K.frame(ax, left=False); ax.set_xlim(0, 70); ax.set_xlabel("estimated FDR (%)")

# ================================================================== d  development cohorts out
ax = C.panel("d", 24, 56, 28, 28, head="Development cohorts out")
aff = sorted(EXC.target.unique())
a = CL[(CL.stat == "fold") & (CL.kind == "target")].set_index("target")
b = CO[(CO.stat == "fold") & (CO.kind == "target")].set_index("target")
aff = sorted(aff, key=lambda s: a.loc[s, "med_delta"])
for i, s in enumerate(aff):
    ax.plot([a.loc[s, "med_delta"], b.loc[s, "med_delta"]], [i, i], color=K.RULE, lw=0.8, zorder=1)
    K.dots(ax, [a.loc[s, "med_delta"]], [i], K.GRAY, s=9, filled=bool(a.loc[s, "call"]))
    K.dots(ax, [b.loc[s, "med_delta"]], [i], K.ORANGE, s=9, filled=bool(b.loc[s, "call"]))
K.ref(ax, 0, axis="v")
ax.set_yticks(range(len(aff))); ax.set_yticklabels([f"{s} ({(EXC.target == s).sum()})" for s in aff], fontsize=K.T_SMALL)
K.frame(ax, left=False); ax.set_ylim(-0.6, len(aff) - 0.4); ax.set_xlabel("median ΔC")
K.legend(ax, [(K.GRAY, "all cohorts"), (K.ORANGE, "development out")], x=0.5, y=0.16, dy=0.09, s=7, size=K.T_SMALL)

# ================================================================== e  PC1 variance
ax = C.panel("e", 68, 56, 28, 28, head="Variance on PC1")
lim = max(ST.pc1_target.max(), ST.pc1_coex_med.max()) * 1.05
for col_, cc in (("pc1_rand_med", K.SLATE), ("pc1_coex_med", K.AQUA)):
    K.contour2d(ax, ST.pc1_target, ST[col_], cc, (0, lim), (0, lim), levels=(0.5, 0.75, 0.9, 0.97), bw=0.18)
    ax.scatter(ST.pc1_target, ST[col_], s=0.7, color=cc, alpha=0.45, lw=0, zorder=3)
K.identity(ax, 0, lim); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("signature"); ax.set_ylabel("decoys (median)")
K.legend(ax, [(K.AQUA, "co-expressed"), (K.SLATE, "random-gene")], x=0.05, y=0.96, dy=0.1, s=7, size=K.T_SMALL)

# ================================================================== f  global null
ax = C.panel("f", 112, 56, 26, 28, head="Permuted outcome")
K.globalnull_hist(ax, GP, "fold", GN); ax.set_ylabel("permutations")
for p_ in ax.patches:
    p_.set_facecolor(K.tint(K.SLATE, .55)); p_.set_edgecolor("white")

# ================================================================== g  ten sets
ax = C.panel("g", 152, 56, 26, 28, head="Ten decoy sets")
for i, (k, cc, lab) in enumerate((("pc1", K.SLATE, "PC1"), ("abs_pc1", K.AQUA, "PC1 |r|"), ("mean", K.ORANGE, "mean z"))):
    v = B[k][B[k].stat == "fold"].fdr.to_numpy()
    K.dots(ax, i + np.linspace(-.22, .22, len(v)), v, cc, s=6)
    ax.plot([i - .3, i + .3], [row(F[k]).fdr] * 2, color=K.INK, lw=1.0, zorder=5)
ax.set_xticks(range(3)); ax.set_xticklabels(["PC1", "abs. r", "mean z"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("estimated FDR (%)"); ax.set_ylim(0, 55); ax.set_xlim(-0.6, 2.6)

# ================================================================== h + i  forest and provenance, aligned
ax = C.panel("h", 32, 97, 140, 22, head="Each signature against its 320 co-expressed decoys (PC1)")
dev = set(EXC.target)
c = CL[CL.stat == "fold"]; t = c[c.kind == "target"].set_index("target")
d = c[(c.kind == "decoy") & (c["class"] == "coex")]
qd = d.groupby("target").med_delta.quantile([.05, .25, .5, .75, .95]).unstack()
order = t.med_delta.sort_values(ascending=False).index
for i, s in enumerate(order):
    if i % 2 == 0:
        ax.axvspan(i - .5, i + .5, color=K.GRID, lw=0, zorder=0)
    K.band(ax, i, q=qd.loc[s].to_numpy(), color=K.AQUA, width=0.9, dot=False)
    ax.plot([i - .28, i + .28], [qd.loc[s, .5]] * 2, color=K.AQUA, lw=0.8, zorder=4)
    K.dots(ax, [i], [t.loc[s, "med_delta"]], K.ORANGE if s in dev else K.INK, s=11 if t.loc[s, "call"] else 8,
           filled=bool(t.loc[s, "call"]), z=5)
K.ref(ax, 0); ax.set_xlim(-0.6, len(order) - 0.4); ax.set_xticks([]); ax.set_ylabel("median ΔC")
K.hlegend(C, 104, 94.0, [(K.INK, "called", "o", True), (K.INK, "not called", "o", False),
                          (K.ORANGE, "developed in a cohort here", "o", True), (K.tint(K.AQUA, .45), "decoys 5–95%", "s"), (K.tint(K.AQUA, .18), "25–75%", "s")],
          size=K.T_SMALL, s=7, gap=1.6)
ax2 = C.ax(32, 121, 140, 19)
C.letter(6.5, 122.5, "i")
code_col = {"D": K.ORANGE, "C": K.tint(K.ORANGE, 0.45), "V": K.LGRAY}
bl = set(zip(BL.target, BL.cohort))
for i, s in enumerate(order):
    for j, co in enumerate(COH):
        code = str(PV.loc[s, co]) if s in PV.index and pd.notna(PV.loc[s, co]) else ""
        if code in code_col:
            ax2.add_patch(Rectangle((i - .45, j - .45), .9, .9, fc=code_col[code], ec="none"))
        if (s, co) in bl:
            ax2.scatter([i], [j], s=1.6, color=K.INK, lw=0, zorder=4)
ax2.set_xlim(-0.6, len(order) - 0.4); ax2.set_ylim(len(COH) - 0.5, -0.5)
ax2.set_yticks(range(len(COH))); ax2.set_yticklabels([cl_(c) for c in COH], fontsize=4.6)
ax2.set_xticks([]); ax2.tick_params(length=0, pad=1)
[sp.set_visible(False) for sp in ax2.spines.values()]
for j in range(len(COH) + 1):
    ax2.axhline(j - 0.5, color=K.GRID, lw=0.4, zorder=0)
K.hlegend(C, 32, 181.0, [(K.ORANGE, "developed (D)", "s"), (K.tint(K.ORANGE, .45), "shaped (C)", "s"), (K.LGRAY, "tested only (V)", "s"),
                         (K.INK, "borderline", "o")], size=K.T_SMALL, s=7, gap=1.8)

# ================================================================== j  gain per cohort, same columns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ax3 = C.ax(32, 142, 140, 19)
C.letter(6.5, 143.5, "j")
UU = U.drop_duplicates(["cohort", "target"]).assign(dC=lambda d_: d_.fold_obs - d_.fold_med)
HM = UU.pivot_table(index="cohort", columns="target", values="dC").reindex(index=COH, columns=order)
HL = 0.1
ax3.imshow(HM.to_numpy(), cmap=K.DIV, norm=Normalize(-HL, HL), aspect="auto", interpolation="nearest")
ax3.set_xlim(-0.6, len(order) - 0.4)
ax3.set_yticks(range(len(COH))); ax3.set_yticklabels([cl_(c) for c in COH], fontsize=4.6)
ax3.set_xticks(range(len(order))); ax3.set_xticklabels(order, rotation=90, fontsize=4.4); ax3.tick_params(length=0, pad=1)
[sp.set_visible(False) for sp in ax3.spines.values()]
K.cbar(C, 174.5, 143, 1.5, 17, K.DIV, -HL, HL, [-.1, 0, .1], label="ΔC", horizontal=False)

# ================================================================== k  Kaplan-Meier example
KM_ = pd.read_csv(D + "v13_fig/km_breast.csv")
C.letter(6.5, 189.0, "k"); C.head(11.5, 189.2, "SOTIRIOU-GGI in METABRIC (RFS), split at the median")
kmax = [K.km_risk(C, 16, 195, 32, 19, KM_, "target", "signature", K.INK, 300, [0, 100, 200, 300], "months"),
        K.km_risk(C, 56, 195, 32, 19, KM_, "decoy", "median decoy", K.AQUA, 300, [0, 100, 200, 300], "months", show_y=False)]

# ================================================================== l  the example across cohorts
ax = C.panel("l", 106, 195, 72, 19, head="SOTIRIOU-GGI and its 32 co-expressed decoys, every cohort")
EXS = "SOTIRIOU-GGI"
Ud = pd.read_csv(P + "breast_pc1_units.csv.gz")
Ud = Ud[(Ud.target == EXS) & (((Ud.kind == "decoy") & (Ud["class"] == "coex") & (Ud.set == 1)) | (Ud.kind == "target"))].copy()
Ud["dC"] = Ud.fold_obs - Ud.fold_med
oc = Ud[Ud.kind == "target"].set_index("cohort").dC.reindex(COH).sort_values(ascending=False).index
for i, co in enumerate(oc):
    dv = Ud[(Ud.cohort == co) & (Ud.kind == "decoy")].dC.to_numpy()
    if i % 2 == 0:
        ax.axvspan(i - .5, i + .5, color=K.GRID, lw=0, zorder=0)
    ax.scatter(i + K.swarm(dv, 0.36, 12), dv, s=1.8, color=K.AQUA, alpha=0.85, lw=0, zorder=2)
    K.dots(ax, [i], [Ud[(Ud.cohort == co) & (Ud.kind == "target")].dC.iloc[0]], K.INK, s=12, z=4)
K.ref(ax, 0)
ax.set_xticks(range(len(oc))); ax.set_xticklabels([cl_(c) for c in oc], rotation=40, ha="right", fontsize=4.6)
ax.tick_params(axis="x", length=0); ax.set_xlim(-0.5, len(oc) - 0.5); ax.set_ylabel("ΔC")
C.save("Figure7_GB10")
print("Figure7_GB10 |", [(r[0], int(r[1].called), round(r[1].fdr, 1), round(r[1].fdr_pi0, 1)) for r in rows], "| affected", len(aff),
      "| null", GN[["stat", "ratio"]].round(2).values.tolist(), "| km P", [round(p_, 4) for _, p_ in kmax],
      "| l target dC", Ud[Ud.kind == "target"].set_index("cohort").dC.round(3).to_dict())
