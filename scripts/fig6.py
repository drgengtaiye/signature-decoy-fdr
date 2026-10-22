"""Figure 6 (GB v10) - the same procedure for prognosis in 20 TCGA cancer types (folded C-index, ten decoy sets).

a  signatures called per endpoint, split into expected true and false discoveries; estimate and pi0-corrected estimate
b  the estimate with intervals over correlation clusters, by endpoint and decoy class
c  decoys' internal correlation against their target's (set 1)
d  the shortfall of decoy co-expression in each cancer type
e  under 200 permutations of (time, event) within type: realized false calls and the estimate
f  one decoy set against ten, per endpoint
g  response and prognosis are different claims: each signature's median delta for both
h  what is called, by family, per endpoint
i  each signature against its 320 co-expressed decoys (overall survival)
j  each signature's gain over its matched null in each cancer type (same columns as i)
k  Kaplan-Meier curves in TCGA KIRP split at the median of CIN70 (Carter) and of its median co-expressed decoy (an
   illustration chosen by effect size; validation/PLAN_round8_figure_displays.md)
l  the same signature's gain in every cancer type against its co-expressed decoys' (first decoy set)
Data: icb/v12_pkg/tcga_* (v12_rerun.R), v11_pkg/tcga_decoy_r_by_type.csv (set 1 decoys), calib_v3_units.csv (response).
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
D = K.D + "/icb/"; P = D + "v12_pkg/"
FAM = pd.read_csv(D + "v5_signature_clusters.csv").set_index("signature").family.to_dict()
EP = ["OS", "DSS", "PFI"]; EPL = {"OS": "overall survival", "DSS": "disease-specific", "PFI": "progression-free"}
F = {e: pd.read_csv(P + f"tcga_{e}_fdr.csv") for e in EP}
B = {e: pd.read_csv(P + f"tcga_{e}_fdr_by_set.csv") for e in EP}
CL = {e: pd.read_csv(P + f"tcga_{e}_calls.csv.gz") for e in EP}
INFO = pd.read_csv(P + "tcga_decoy_info.csv.gz")
RP = pd.read_csv(D + "v11_pkg/tcga_decoy_r_by_type.csv")
GN = pd.read_csv(P + "tcga_OS_globalnull_summary.csv"); GP = pd.read_csv(P + "tcga_OS_globalnull_per_perm.csv")
U0 = pd.read_csv(D + "calib_v3_units.csv"); U0 = U0[(U0["null"] == "N0") & (U0.kind == "target") & (U0.stat == "fold")]


def row(df, decoys="coex", boot="clusters", st="fold"):
    return df[(df.stat == st) & (df.decoys == decoys) & (df.boot == boot)].iloc[0]


C = K.Canvas(224)

# ================================================================== a  calls per endpoint
ax = C.panel("a", 14, 12, 35, 30, head="Calls per endpoint")
for i, e in enumerate(EP):
    r = row(F[e])
    K.stacked_calls(ax, i, r.called, r.exp_false, K.SLATE, width=0.6, text=f"{r.fdr:.0f}%\n({r.fdr_pi0:.0f}%)")
ax.set_xticks(range(3)); ax.set_xticklabels(EP); ax.tick_params(axis="x", length=0)
ax.set_ylabel("signatures called (of 52)"); ax.set_ylim(0, 32); ax.set_xlim(-0.6, 2.6)
K.legend(ax, [(K.SLATE, "expected true", "box"), (K.tint(K.SLATE, .62), "expected false", "box")], x=0.03, y=0.98, dy=0.09, s=8,
         size=K.T_SMALL)
ax.text(0.03, 0.78, "% estimate (π0-corrected)", transform=ax.transAxes, fontsize=K.T_SMALL, color=K.INK2, va="top")

# ================================================================== b  estimate with intervals
ax = C.panel("b", 62, 12, 34, 35, head="Estimate and interval")
rows = [(e, dc) for e in EP for dc in ("coex", "coex closest", "rand")]
DL = {"coex": "co-expr.", "coex closest": "closest", "rand": "random"}
for k, (e, dc) in enumerate(rows):
    y = k + 0.6 * EP.index(e)
    r = row(F[e], dc)
    K.hbar_interval(ax, y, r.fdr_lo, min(r.fdr_hi, 150), K.tint(K.SLATE, .45), lw=1.4)
    K.dots(ax, [r.fdr], [y], K.SLATE, s=10, z=5); K.dots(ax, [r.fdr_pi0], [y], K.SLATE, s=10, filled=False, z=5)
    if r.fdr_hi > 150:
        ax.text(151, y, "›", va="center", fontsize=K.T_NOTE, color=K.INK2)
ys = [k + 0.6 * EP.index(e) for k, (e, dc) in enumerate(rows)]
ax.set_yticks(ys); ax.set_yticklabels([f"{e} {DL[dc]}" for e, dc in rows], fontsize=K.T_SMALL)
ax.set_ylim(max(ys) + 0.6, -0.6); K.frame(ax, left=False); ax.set_xlim(0, 150); ax.set_xticks([0, 50, 100, 150])
ax.set_xlabel("estimated FDR (%)")

# ================================================================== c  decoy co-expression
ax = C.panel("c", 112, 12, 26, 35, head="Decoy co-expression")
for cls, cc in (("rand", K.SLATE), ("coex", K.AQUA)):
    g = INFO[(INFO["class"] == cls) & (INFO.set == 1)]
    K.contour2d(ax, g.r_target_med, g.r_achieved_med, cc, (0, 1), (-0.05, 1), levels=(0.5, 0.75, 0.9, 0.97), bw=0.15)
    ax.scatter(g.r_target_med, g.r_achieved_med, s=0.7, color=cc, alpha=0.45, lw=0, zorder=3)
K.identity(ax, 0, 1)
ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1); ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
ax.set_xlabel("signature r"); ax.set_ylabel("decoy r")
K.legend(ax, [(K.AQUA, "co-expressed"), (K.SLATE, "random-gene")], x=0.05, y=0.96, dy=0.1, s=7, size=K.T_SMALL)

# ================================================================== d  shortfall by type
ax = C.panel("d", 150, 12, 28, 35, head="Shortfall by type")
RP["gap"] = RP.r_target - RP.r_decoy_med
order = RP.groupby("cohort").gap.median().sort_values().index
for i, co in enumerate(order):
    v = RP[RP.cohort == co].gap
    K.band(ax, i, v=v, color=K.SLATE, orient="h", width=0.55, dot_s=6)
K.ref(ax, 0, axis="v")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=4.8); K.frame(ax, left=False)
ax.set_ylim(-0.6, len(order) - 0.4); ax.set_xlabel("signature − decoy r")

# ================================================================== e  global null
ax = C.panel("e", 14, 62, 30, 28, head="Permuted survival")
K.globalnull_hist(ax, GP, "fold", GN)
for p_ in ax.patches:
    p_.set_facecolor(K.tint(K.SLATE, .55)); p_.set_edgecolor("white")
ax.set_ylabel("permutations")
K.legend(ax, [(K.INK, "realized", "line"), (K.AQUA, "estimate", "line")], x=0.55, y=0.62, dy=0.1, size=K.T_SMALL)

# ================================================================== f  ten sets
ax = C.panel("f", 62, 62, 30, 28, head="Ten decoy sets")
for i, e in enumerate(EP):
    b = B[e][B[e].stat == "fold"].fdr.to_numpy()
    K.dots(ax, i + np.linspace(-.22, .22, len(b)), b, K.SLATE, s=6)
    ax.plot([i - .3, i + .3], [row(F[e]).fdr] * 2, color=K.INK, lw=1.0, zorder=5)
ax.set_xticks(range(3)); ax.set_xticklabels(EP); ax.tick_params(axis="x", length=0)
ax.set_ylabel("estimated FDR (%)"); ax.set_ylim(0, 75); ax.set_xlim(-0.6, 2.6)

# ================================================================== g  response vs prognosis
ax = C.panel("g", 108, 62, 30, 28, head="Response vs survival")
ost = CL["OS"][(CL["OS"].stat == "fold") & (CL["OS"].kind == "target")].set_index("target")
rsp = U0.set_index("unit")
j = pd.DataFrame({"resp": rsp.med_delta, "os": ost.med_delta, "rc": rsp.call.astype(str) == "True", "oc": ost.call}).dropna()
for s, r in j.iterrows():
    filled = bool(r.rc) and bool(r.oc)
    K.dots(ax, [r.resp], [r.os], K.FAM[FAM[s]], s=10 if (r.rc or r.oc) else 6, filled=filled or not (r.rc or r.oc),
           alpha=1 if (r.rc or r.oc) else 0.45)
K.ref(ax, 0); K.ref(ax, 0, axis="v")
ax.set_xlabel("median ΔAUC, response"); ax.set_ylabel("median ΔC, survival")
both = int((j.rc & j.oc).sum())
ax.text(0.97, 0.97, f"filled: called for both ({both})\nopen: for one; faint: neither", transform=ax.transAxes, ha="right",
        va="top", fontsize=K.T_SMALL - 0.3, color=K.INK2, linespacing=1.05)

# ================================================================== h  by family
ax = C.panel("h", 152, 62, 26, 28, head="Called, by family")
for i, e in enumerate(EP):
    t = CL[e][(CL[e].stat == "fold") & (CL[e].kind == "target") & CL[e].call]
    bottom = 0
    for f in ("Inflamed/IFN", "Other", "Confounder"):
        n = int((t.target.map(FAM) == f).sum())
        ax.bar(i, n, 0.6, bottom=bottom, color=K.FAM[f], edgecolor="white", linewidth=0.5, zorder=3); bottom += n
ax.set_xticks(range(3)); ax.set_xticklabels(EP); ax.tick_params(axis="x", length=0)
ax.set_ylabel("signatures called"); ax.set_ylim(0, 27); ax.set_xlim(-0.6, 2.6)
K.legend(ax, [(K.FAM[f], K.FAMLAB[f], "box") for f in ("Inflamed/IFN", "Other", "Confounder")], x=0.03, y=0.98, dy=0.09, s=7,
         size=K.T_SMALL)

# ================================================================== i  forest
ax = C.panel("i", 22, 106, 150, 24, head="Each signature against its 320 co-expressed decoys (folded, OS)")
order = K.forest_bands(ax, CL["OS"], "fold", fam=FAM)
ax.set_xticklabels([]); ax.set_ylabel("median ΔC")
K.hlegend(C, 96, 103.0, [(K.BLUE, "inflamed/IFN"), (K.GRAY, "other"), (K.ORANGE, "confounder"), (K.INK, "called", "o", True),
                          (K.INK, "not called", "o", False), (K.tint(K.AQUA, .45), "decoys 5–95%", "s"), (K.tint(K.AQUA, .18), "25–75%", "s")],
          size=K.T_SMALL, s=7, gap=1.5)
# ================================================================== j  gain per cancer type, same columns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
UO = pd.read_csv(P + "tcga_OS_units.csv.gz"); UO = UO[UO.kind == "target"].drop_duplicates(["cohort", "target"])
UO["dC"] = UO.fold_obs - UO.fold_med
HM = UO.pivot_table(index="cohort", columns="target", values="dC")
tord = HM.median(axis=1).sort_values(ascending=False).index
HM = HM.reindex(index=tord, columns=order)
ax3 = C.ax(22, 132, 150, 26); C.letter(6.5, 133.5, "j")
HL = 0.1
ax3.imshow(HM.to_numpy(), cmap=K.DIV, norm=Normalize(-HL, HL), aspect="auto", interpolation="nearest")
ax3.set_xlim(-0.6, len(order) - 0.4)
ax3.set_yticks(range(len(tord))); ax3.set_yticklabels(tord, fontsize=4.2)
ax3.set_xticks(range(len(order))); ax3.set_xticklabels([K.lab(x_) for x_ in order], rotation=90, fontsize=4.2)
ax3.tick_params(length=0, pad=1); [sp.set_visible(False) for sp in ax3.spines.values()]
K.cbar(C, 174.5, 134, 1.5, 20, K.DIV, -HL, HL, [-.1, 0, .1], label="ΔC", horizontal=False)

# ================================================================== k  Kaplan-Meier example
KM_ = pd.read_csv(D + "v13_fig/km_tcga.csv"); KM_["time"] = KM_.time / 365.25
C.letter(6.5, 184.0, "k"); C.head(11.5, 184.2, "CIN70 (Carter) in kidney papillary (KIRP), overall survival")
kmax = [K.km_risk(C, 16, 190, 32, 19, KM_, "target", "signature", K.ORANGE, 10, [0, 2.5, 5, 7.5, 10], "years"),
        K.km_risk(C, 56, 190, 32, 19, KM_, "decoy", "median decoy", K.AQUA, 10, [0, 2.5, 5, 7.5, 10], "years", show_y=False)]

# ================================================================== l  the example across cancer types
ax = C.panel("l", 106, 190, 72, 19, head="CIN70 (Carter) and its 32 co-expressed decoys, every cancer type")
EXS = "CIN70_Carter"
Ud = pd.read_csv(P + "tcga_OS_units.csv.gz")
Ud = Ud[(Ud.target == EXS) & (((Ud.kind == "decoy") & (Ud["class"] == "coex") & (Ud.set == 1)) | (Ud.kind == "target"))].copy()
Ud["dC"] = Ud.fold_obs - Ud.fold_med
oc = Ud[Ud.kind == "target"].drop_duplicates("cohort").set_index("cohort").dC.sort_values(ascending=False).index
for i, co in enumerate(oc):
    dv = Ud[(Ud.cohort == co) & (Ud.kind == "decoy")].dC.to_numpy()
    if i % 2 == 0:
        ax.axvspan(i - .5, i + .5, color=K.GRID, lw=0, zorder=0)
    ax.scatter(i + K.swarm(dv, 0.36, 12), dv, s=1.6, color=K.AQUA, alpha=0.85, lw=0, zorder=2)
    K.dots(ax, [i], [Ud[(Ud.cohort == co) & (Ud.kind == "target")].dC.iloc[0]], K.ORANGE, s=12, z=4)
K.ref(ax, 0)
ax.set_xticks(range(len(oc))); ax.set_xticklabels(oc, rotation=90, fontsize=4.2)
ax.tick_params(axis="x", length=0, pad=1); ax.set_xlim(-0.5, len(oc) - 0.5); ax.set_ylabel("ΔC")
C.save("Figure6_GB10")
print("Figure6_GB10 |", {e: (int(row(F[e]).called), round(row(F[e]).fdr, 1), round(row(F[e]).fdr_pi0, 1)) for e in EP},
      "| both", both, "| null", GN[["stat", "ratio"]].round(2).values.tolist(), "| km P", [round(p_, 5) for _, p_ in kmax],
      "| l n types", len(oc))
