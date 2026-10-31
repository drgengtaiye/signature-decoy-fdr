"""Figure 5 (GB v11) - what the estimate estimates, how far it sits from it, and the standard estimates beside it.

a  where every signature is null (1,000 within-cohort label permutations): realized false calls per permutation, the
   plain and pi0-corrected decoy estimates and the nominal expectation m x 0.025, folded and directed
b  per signature under the same permutations: its own call rate against its decoys'
c  where every call is false, in every analysis: estimated over realized false calls, nominal and decoy estimates
d  planted signal in the response cohorts, truth by gene content (>= 10% of a signature's genes causal)
e  the same replicates against the formal truth (called more often than its decoys over 200 further label draws)
f  planted survival signal in 20 TCGA cancer types with new causal gene sets, formal truth
g  where the plain estimate's excess comes from
h  replicate by replicate: realized against estimated FDR (formal truth; folded, beta = 1, 12 cohorts)
i  the estimate minus the realized FDR in every replicate, by scenario and statistic
j  the nominal, permutation and decoy estimates over the realized FDR in each of the 22 planted settings
Ratios of totals over replicates (200 per response setting, 100 per TCGA setting). validation/PLAN_round7_estimand_and_bias.md;
c and j: validation/PLAN_round9_standard_fdr_comparison.md. Panels of the v10 figure moved to Fig. S2 (figS2_v11.py).
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import glob, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
D = K.D + "/icb/"; T = D + "v12_truth/"; S14 = D + "v14_standard/"
P0 = pd.read_csv(T + "null_pi0.csv"); SN = pd.read_csv(T + "summary_null.csv").set_index("stat")
NS = pd.read_csv(D + "v10_calib_null_summary.csv").set_index(["stat", "estimator"])
PT = pd.read_csv(D + "v10_calib_null_per_target.csv", index_col=0)
SG = pd.read_csv(D + "v10_calib_signal.csv"); SG = SG[SG.truth == "10%"]
SR = pd.read_csv(T + "summary_response.csv"); ST = pd.read_csv(T + "summary_tcga.csv")
SIG = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(T + "signal_*.csv"))])
GN = pd.read_csv(S14 + "standard_globalnull.csv"); SP = pd.read_csv(S14 + "standard_planted.csv")
SCOL = {"P1 inflamed": K.BLUE, "P2 non-inflamed": K.ORANGE, "P3 outside signatures": K.AQUA}
SLAB = {"P1 inflamed": "3 inflamed signatures", "P2 non-inflamed": "3 other signatures", "P3 outside signatures": "outside every signature"}
TCOL = {"T1 library": K.BLUE, "T2 outside": K.AQUA, "T3 mixed": K.ORANGE}
TLAB = {"T1 library": "2 signatures", "T2 outside": "1 Hallmark set (out of library)", "T3 mixed": "both"}
NOM = K.GRAY                                                       # the nominal (Benjamini-Hochberg) expectation

C = K.Canvas(210)
OV = C.fig.add_axes([0, 0, 1, 1], facecolor="none"); OV.set_xlim(0, C.w); OV.set_ylim(C.h, 0); OV.axis("off")
OV.set_zorder(20)


def lkey(x, y, color, ls, label, lw=0.9):
    """One legend entry drawn as a short line in figure millimetres, in the line style it stands for."""
    OV.plot([x, x + 3.2], [y, y], color=color, lw=lw, ls=ls)
    OV.text(x + 4.0, y, label, fontsize=K.T_SMALL, va="center", ha="left", color=K.INK)


# ================================================================== a  global null
for gi, st in enumerate(("fold", "dir")):
    x = P0[P0.stat == st]
    ax = C.panel("a", 14, 12, 28, 27, head="All signatures null") if gi == 0 else C.ax(48, 12, 28, 27)
    bins = np.arange(-0.5, max(x.R.max(), 8) + 1.5, 1)
    K.hist_steps(ax, x.R, bins, K.SLATE, alpha=0.35)
    ax.axvline(x.E.mean(), color=K.AQUA, lw=1.0, zorder=4)
    ax.axvline(x.E0.mean(), color=K.AQUA, lw=0.8, ls=(0, (1.2, 1.2)), zorder=4)
    ax.axvline(52 * 0.025, color=NOM, lw=0.9, ls=(0, (3, 1.5)), zorder=4)
    ax.axvline(x.R.mean(), color=K.INK, lw=1.0, zorder=3)
    ax.set_xlim(-0.6, 8.6); ax.set_xticks(range(0, 9, 2))
    ax.set_xlabel("false calls per library")
    if gi == 0:
        ax.set_ylabel("permutations (share)")
    ax.text(0.97, 0.97, K.STATLAB[st], transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL)
lkey(14, 49.0, K.INK, "-", "realized (mean)")
lkey(38, 49.0, K.AQUA, "-", "decoy estimate")
lkey(14, 52.2, K.AQUA, (0, (1.2, 1.2)), "π0-corrected")
lkey(38, 52.2, NOM, (0, (3, 1.5)), "nominal, 52 × 2.5%")

# ================================================================== b  per signature under permutation
ax = C.panel("b", 92, 12, 26, 27, head="Per signature")
for st, mk in (("fold", "o"), ("dir", "^")):
    K.dots(ax, 100 * PT[f"{st}_decoy_callrate"], 100 * PT[f"{st}_target_callrate"], K.STAT[st], s=7, marker=mk)
K.identity(ax, 0, 9)
ax.set_xlim(0, 8); ax.set_ylim(0, 8); ax.set_xticks([0, 4, 8]); ax.set_yticks([0, 4, 8])
ax.set_xlabel("decoys called (%)"); ax.set_ylabel("signature called (%)")
rf, rd = NS.loc[("fold", "coex A"), "per_target_r"], NS.loc[("dir", "coex A"), "per_target_r"]
K.legend(ax, [(K.SLATE, f"folded, r = {rf:.2f}"), (K.BLUE, f"directed, r = {rd:.2f}", "^")], x=0.06, y=0.95, dy=0.1, s=7,
         size=K.T_SMALL)

# ================================================================== c  global null, every analysis: nominal and decoy
ax = C.ax(146, 12, 32, 27); C.letter(124.5, 6.8, "c"); C.head(129.5, 7.0, "Where every call is false")
ALAB = {("response", "fold"): "response, folded", ("response", "dir"): "response, directed",
        ("tcga_OS", "fold"): "TCGA survival", ("breast_pc1", "fold"): "breast, PC1", ("hallmark_OS", "fold"): "Hallmark"}
lg = np.log10
ax.axvspan(lg(0.8), lg(1.25), color=K.GRID, lw=0, zorder=0)
ax.axvline(0, color=K.RULE, lw=0.6, zorder=1)
for i, r in enumerate(GN.itertuples()):
    for off, v, lo, hi, col, filled, mk in ((-0.17, r.decoy_ratio, r.decoy_lo, r.decoy_hi, K.AQUA, True, "o"),
                                            (0.17, r.nominal_ratio, r.nominal_lo, r.nominal_hi, NOM, False, "D")):
        K.hbar_interval(ax, i + off, lg(lo), lg(hi), col, lw=0.9)
        K.dots(ax, [lg(v)], [i + off], col, s=10 if mk == "o" else 7, filled=filled, marker=mk, z=5)
ax.set_yticks(range(len(GN))); ax.set_yticklabels([ALAB[(a, s)] for a, s in zip(GN.analysis, GN.stat)])
ax.set_ylim(len(GN) - 0.45, -0.55); K.frame(ax, left=False); ax.tick_params(axis="y", length=0)
tk = [0.5, 1, 2]
ax.set_xlim(lg(0.38), lg(2.1)); ax.set_xticks(lg(tk)); ax.set_xticklabels(["0.5", "1", "2"]); ax.minorticks_off()
ax.set_xlabel("estimated / realized false calls")
for xx, (col, lab_, mk, filled) in ((136.0, (K.AQUA, "decoy estimate", "o", True)), (157.0, (NOM, "nominal", "D", False))):
    if filled:
        OV.scatter([xx], [49.0], s=9, marker=mk, c=[col], edgecolors="white", linewidths=0.3)
    else:
        OV.scatter([xx], [49.0], s=7, marker=mk, facecolors="white", edgecolors=col, linewidths=0.7)
    OV.plot([xx - 1.6, xx + 1.6], [49.0, 49.0], color=col, lw=0.9, zorder=-1)
    OV.text(xx + 2.4, 49.0, lab_, fontsize=K.T_SMALL, va="center", ha="left", color=K.INK)
OV.text(136.0 - 1.6, 52.2, "line, 95% interval; band, 0.8–1.25", fontsize=K.T_SMALL, va="center", ha="left", color=K.INK2)
OV.text(146 + 32 * (lg(1.25) - lg(0.38)) / (lg(2.1) - lg(0.38)) + 0.6, 13.2, "set in\nadvance", fontsize=K.T_SMALL - 0.4,
        color=K.INK2, va="top", ha="left", linespacing=1.0)


# ================================================================== d, e, f  planted signal, log axes
def logax(ax, lo=0.05, hi=300, ylab=True):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo * 20, hi)
    ax.set_xticks([0.1, 1, 10, 100]); ax.set_xticklabels(["0.1", "1", "10", "100"])
    ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"])
    ax.minorticks_off()
    K.identity(ax, lo, hi)
    ax.set_xlabel("realized FDR (%)")
    if ylab:
        ax.set_ylabel("estimate (%)")


ax = C.panel("d", 14, 64, 34, 31, head="Truth by gene content")
for (sc, b, c, st), g in SG.groupby(["scenario", "beta", "cohorts", "stat"]):
    xr, ye = 100 * g.V.sum() / g.R.sum(), 100 * g["E coex A"].sum() / g["R coex A"].sum()
    K.dots(ax, [xr], [ye], SCOL[sc], s=13 if c == 12 else 7, marker="o" if st == "fold" else "^", alpha=1 if b == 1 else 0.55)
logax(ax)
ax.text(0.04, 0.96, "above the line:\nestimate high", transform=ax.transAxes, fontsize=K.T_SMALL - 0.3, color=K.INK2, va="top")

ax = C.panel("e", 62, 64, 34, 31, head="Formal truth")
for r in SR.itertuples():
    for v, filled in ((r.estimate, True), (r.estimate_pi0, False)):
        K.dots(ax, [max(r.realized, 0.06)], [v], SCOL[r.scenario], s=13 if r.cohorts == 12 else 7, filled=filled,
               marker="o" if r.stat == "fold" else "^", alpha=1 if r.beta == 1 else 0.6)
logax(ax, ylab=False)

ax = C.panel("f", 110, 64, 34, 31, head="TCGA survival, new causal sets")
for r in ST.itertuples():
    for v, filled in ((r.estimate, True), (r.estimate_pi0, False)):
        K.dots(ax, [r.realized], [v], TCOL[r.kind], s=14 if r.beta == 0.6 else 8, filled=filled, marker="s")
logax(ax, ylab=False)
ax.text(0.04, 0.96, "TCGA decoys, 0.34 vs 0.47\nco-expression", transform=ax.transAxes, fontsize=K.T_SMALL - 0.3,
        color=K.INK2, va="top", linespacing=1.05)

# one vertical key for the scenario colours and markers of d-j
KEY = ([(SCOL[k], f"response: {SLAB[k]}", "o", True) for k in SCOL]
       + [(TCOL[k], f"TCGA: {TLAB[k]}", "s", True) for k in TCOL]
       + [(K.INK, "plain estimate", "o", True), (K.INK, "π0-corrected", "o", False), (K.INK, "directed", "^", True)])
for i, (c, lab_, mk, filled) in enumerate(KEY):
    yy = 63.5 + i * 3.45 + (1.2 if i >= 3 else 0) + (1.2 if i >= 6 else 0)
    if filled:
        OV.scatter([151.5], [yy], s=8, marker=mk, c=[c], edgecolors="white", linewidths=0.3)
    else:
        OV.scatter([151.5], [yy], s=8, marker=mk, facecolors="white", edgecolors=c, linewidths=0.7)
    OV.text(153.3, yy, lab_, fontsize=K.T_SMALL - 0.2, va="center", ha="left", color=K.INK)

# ================================================================== g  decomposition
ax = C.panel("g", 27, 118, 82, 34, head="Where the excess comes from")
rows = []
for r in SR[(SR.beta == 1.0) & (SR.cohorts == 12)].itertuples():
    rows.append((f"{r.scenario.split()[0]} {K.STATLAB[r.stat]}", r, SCOL[r.scenario]))
for r in ST[ST.beta == .6].itertuples():
    rows.append((f"TCGA {r.kind.split()[0]}", r, TCOL[r.kind]))
for i, (lab_, r, cc) in enumerate(rows):
    left = 0
    for v, col in ((r.exp_V, K.INK), (r.excess_null, cc), (r.part_true, K.tint(cc, 0.6))):
        ax.barh(i, v, 0.62, left=left, color=col, edgecolor="white", linewidth=0.5, zorder=3)
        left += v
    K.dots(ax, [r.estimate_pi0], [i], K.INK, s=11, filled=False, marker="D", z=6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([x[0] for x in rows]); ax.set_ylim(len(rows) - 0.4, -0.6)
K.frame(ax, left=False); ax.set_xlabel("% of calls (expected)"); ax.set_xlim(0, 240)
# key: the three parts in one example hue (P3), as each bar carries its own scenario's hue
ex = K.AQUA
K.legend(ax, [(K.INK, "false calls of null signatures", "box"),
              (ex, "decoys' excess over null signatures (scenario hue)", "box"),
              (K.tint(ex, 0.6), "decoy rates of true signatures (light)", "box"),
              (K.INK, "π0-corrected estimate", "D", False)], x=0.60, y=0.31, dy=0.085, s=8, size=K.T_SMALL)

# ================================================================== h  replicate level
ax = C.ax(133, 121, 36, 31); C.letter(125.5, 112.8, "h"); C.head(130.5, 113.0, "Replicate by replicate")
sub = SIG[(SIG.stat == "fold") & (SIG.beta == 1.0) & (SIG.cohorts == 12) & (SIG.R > 0)].copy()
sub["est"] = 100 * sub.E / sub.R; sub["real"] = 100 * sub.V_cp / sub.R
FL, TOP = 0.5, 500.0                                               # log axes; a realized FDR of 0 is drawn at 0.5%
lgc = lambda v: np.log10(np.clip(v, FL, TOP))
for sc in SCOL:
    g = sub[sub.scenario == sc]
    ax.scatter(lgc(g.est), lgc(g.real) + np.where(g.real <= 0, np.random.default_rng(2).uniform(-0.06, 0.06, len(g)), 0),
               s=2.2, color=SCOL[sc], alpha=0.45, lw=0, zorder=3)
K.identity(ax, lgc(FL), lgc(TOP))
tk = np.array([1, 10, 100])
ax.set_xlim(lgc(FL) - 0.1, lgc(TOP)); ax.set_ylim(lgc(FL) - 0.12, lgc(TOP))
ax.set_xticks(np.log10(tk)); ax.set_xticklabels(["1", "10", "100"])
ax.set_yticks(np.r_[lgc(FL), np.log10(tk)]); ax.set_yticklabels(["0", "1", "10", "100"])
ax.axhline(lgc(FL) + 0.15, color=K.GRID, lw=0.4, zorder=1)
ax.set_xlabel("estimated FDR (%)"); ax.set_ylabel("realized FDR (%)")
K.marginal(C, 133, 121, 36, 31, [lgc(sub[sub.scenario == sc].est) for sc in SCOL],
           [lgc(sub[sub.scenario == sc].real) for sc in SCOL], [SCOL[sc] for sc in SCOL], size=3.5,
           xlim=(lgc(FL) - 0.1, lgc(TOP)), ylim=(lgc(FL) - 0.12, lgc(TOP)))
above = 100 * (sub.est >= sub.real).mean()
ax.text(0.04, 0.97, f"estimate ≥ realized in\n{above:.0f}% of {len(sub)} replicates", transform=ax.transAxes, ha="left",
        va="top", fontsize=K.T_SMALL, color=K.INK2, linespacing=1.1)

# ================================================================== i  estimate minus realized
ax = C.panel("i", 14, 172, 68, 27, head="Estimate minus realized, every replicate")
allr = SIG[SIG.R > 0].copy(); allr["gap"] = 100 * (allr.E - allr.V_cp) / allr.R
CAP = 250
pos_ = 0; ticks, tlabs = [], []
for sc in SCOL:
    for st in ("fold", "dir"):
        v = allr[(allr.scenario == sc) & (allr.stat == st)].gap
        K.violin(ax, np.clip(v, None, CAP), pos_, SCOL[sc], width=0.4, alpha=0.85 if st == "fold" else 0.4)
        ax.text(pos_, CAP + 8, f"{100 * (v >= 0).mean():.0f}%", ha="center", va="bottom", fontsize=K.T_SMALL, color=K.INK2)
        ticks.append(pos_); tlabs.append(f"{sc.split()[0]} {'F' if st == 'fold' else 'D'}"); pos_ += 1
    pos_ += 0.4
K.ref(ax, 0)
ax.set_xticks(ticks); ax.set_xticklabels(tlabs, fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylabel("percentage points"); ax.set_ylim(-60, CAP + 32); ax.set_yticks([0, 100, 200])
ax.set_xlim(-0.6, pos_ - 0.8)

# ================================================================== j  standard estimates against the realized FDR
ax = C.panel("j", 102, 172, 76, 27, head="Standard estimates beside the decoys, 22 planted settings")
EST = [("nominal", "nominal\n(BH)"), ("permutation", "permutation\n(SAM)"), ("decoy", "decoys"), ("decoy_pi0", "decoys,\nπ0-corrected")]
ax.axhspan(lg(0.25), 0, color=K.tint(K.ORANGE, 0.9), lw=0, zorder=0)
ax.axhline(0, color=K.RULE, lw=0.6, zorder=1)
for xi, (col, lab_) in enumerate(EST):
    s = SP[SP[col].notna()].copy()
    s["ratio"] = s[col] / s.realized
    rng = np.random.default_rng(31 + xi)
    for r in s.itertuples():
        hue = SCOL.get(r.scenario) or TCOL[r.scenario.replace("TCGA ", "")]
        mk = "s" if r.data == "tcga" else ("o" if r.stat == "fold" else "^")
        K.dots(ax, [xi + rng.uniform(-0.22, 0.22)], [lg(r.ratio)], hue, s=9 if r.cohorts != 6 else 5, marker=mk,
               alpha=1 if r.beta in (1.0, 0.6) else 0.6, z=4)
    nb = int((s.ratio < 1).sum())
    ax.text(xi, lg(900), f"{nb} of {len(s)}\nbelow", ha="center", va="top", fontsize=K.T_SMALL, color=K.INK2 if nb == 0 else K.ORANGE,
            linespacing=1.0)
ax.set_xticks(range(len(EST))); ax.set_xticklabels([e[1] for e in EST], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.6, len(EST) - 0.4)
ax.set_ylim(lg(0.25), lg(900)); yt = [0.3, 1, 10, 100]
ax.set_yticks(lg(yt)); ax.set_yticklabels(["0.3", "1", "10", "100"]); ax.minorticks_off()
ax.set_ylabel("estimate / realized FDR")
ax.text(len(EST) - 0.45, lg(0.29), "estimate too low", ha="right", va="bottom", fontsize=K.T_SMALL - 0.3, color=K.ORANGE)

C.save("Figure5_GB11")
print("Figure5_GB11 | null", SN[["ratio", "ratio_pi0"]].round(3).to_dict(), "| h above", round(above, 1), len(sub),
      "| i", allr.groupby(["scenario", "stat"]).gap.apply(lambda v: round(100 * (v >= 0).mean(), 1)).to_dict(),
      "| c", GN[["analysis", "stat", "nominal_ratio", "decoy_ratio"]].round(2).values.tolist(),
      "| j below", {c: int((SP[c] < SP.realized).sum()) for c, _ in EST})
