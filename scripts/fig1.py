"""Figure 1 (GB v10) - a matched null for 52 signatures, and the target-decoy procedure built on it.

a  the procedure (schematic; the matrix, densities and points in it are illustrative, not data)
b  the 52 signatures by size and scoring engine, with the sizes of the benchmark's null ladder
c  the null's 95th percentile moves with set size and engine
d  matched against a fixed 48-gene null, signature by signature
e  internal co-expression per cohort: signatures, their co-expressed and random-gene decoys (first package set)
f  one signature in every cohort: the density of its 1,000-draw null bank, its 32 co-expressed decoys and its own
   folded AUC; Gide marked with its delta
g  the same signature's delta in every cohort with its decoys', and a confounder signature's
h  the same signature's ROC curve in Gide against its decoys' and its null bank's (5th-95th percentile)
i  the width of the matched null against cohort size, with the Mann-Whitney expectation
j  per cohort, the share of tests clearing a fixed 48-gene null against the matched null
Data as in fig1_cns.py.
"""
import sys, json
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/decoys")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy.stats import spearmanr, gaussian_kde
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Polygon, Rectangle
import cns10 as K
from v10_calib_lib import signed_rank_p
K.apply()
D = K.D

R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med
O = pd.read_csv(f"{D}/icb/native_response.csv"); O = O[O.cohort != "ICB_VanDenEnde"]
FIX = O.groupby("cohort")[["rand_med", "rand_p95"]].first()
X = pd.read_csv(f"{D}/icb/v3_dissociation.csv"); famd = X.set_index("signature").family
RB = pd.read_csv(f"{D}/icb/v13_fig/r_by_cohort_response.csv")
EX = pd.read_csv(f"{D}/icb/v13_fig/example_units.csv")

C = K.Canvas(198)

# ================================================================== a  schematic
ax = C.ax(4, 3, 174, 46); ax.set_xlim(0, 174); ax.set_ylim(46, 0); ax.axis("off")
C.letter(1.5, 2.5, "a")
STG = [(0, 45), (49, 40), (93, 38), (135, 39)]
HDR = ["Targets and matched decoys", "Read against a matched null", "Called over cohorts", "FDR of the claim"]
for k, ((x, w), h) in enumerate(zip(STG, HDR)):
    ax.add_patch(FancyBboxPatch((x, 1), w, 44, boxstyle="round,pad=0,rounding_size=1.6", fc="#f7f8f9", ec=K.RULE, lw=0.5, zorder=0))
    ax.add_patch(Circle((x + 3.6, 5.0), 2.1, fc=K.INK, ec="none", zorder=3))
    ax.text(x + 3.6, 5.05, str(k + 1), color="white", fontsize=6.0, fontweight="bold", ha="center", va="center", zorder=4)
    ax.text(x + 7.0, 5.0, h, fontsize=K.T_HEAD, fontweight="bold", va="center")
for (x0, w0), (x1, _) in zip(STG[:-1], STG[1:]):
    ax.add_patch(FancyArrowPatch((x0 + w0 + 0.3, 25), (x1 - 0.3, 25), arrowstyle="-|>", mutation_scale=7, lw=0.8,
                                 color=K.SLATE, shrinkA=0, shrinkB=0, zorder=2))

# stage 1: an illustrative expression matrix - rows are genes, columns patients
rng = np.random.default_rng(11)
npat, nr = 20, 5
f1, f2 = rng.normal(0, 1, npat), rng.normal(0, 1, npat)
blocks = [(K.BLUE, "signature", f1, 0.85), (K.AQUA, "co-expressed decoy", f2, 0.80), (K.GRAY, "random-gene decoy", None, 0)]
x0m, y0m, cw, rh = 5.0, 10.5, 1.05, 1.15
for b, (c, name, f, load) in enumerate(blocks):
    yb = y0m + b * (nr * rh + 1.6)
    M = (load * f[None, :] + np.sqrt(1 - load ** 2) * rng.normal(0, 1, (nr, npat))) if f is not None else rng.normal(0, 1, (nr, npat))
    for i in range(nr):
        for j in range(npat):
            ax.add_patch(Rectangle((x0m + j * cw, yb + i * rh), cw, rh, fc=K.DIV((np.clip(M[i, j], -2, 2) + 2) / 4), ec="white",
                                   lw=0.15, zorder=2))
    ax.add_patch(Rectangle((x0m - 1.6, yb), 1.0, nr * rh, fc=c, ec="none", zorder=2))
    ax.text(x0m + npat * cw + 1.2, yb + nr * rh / 2 - 0.9, name, fontsize=K.T_SMALL, va="center", color=K.INK)
    ax.text(x0m + npat * cw + 1.2, yb + nr * rh / 2 + 1.3,
            ["published genes", "× 32, co-expressed", "× 8, random genes"][b], fontsize=K.T_SMALL - 0.6, va="center", color=K.INK2)
ax.text(x0m, 10.0, "patients", fontsize=K.T_SMALL - 0.7, color=K.INK2, va="bottom")
ax.add_patch(FancyArrowPatch((x0m + 7.2, 9.2), (x0m + 12.5, 9.2), arrowstyle="-|>", mutation_scale=3.5, lw=0.4, color=K.INK2))
ax.text(x0m - 2.0, 43.5, "decoy genes come from no signature; size, engine, weights,\nexpression and co-expression matched; 10 sets",
        fontsize=K.T_SMALL - 0.6, color=K.INK2, va="bottom", linespacing=1.0)

# stage 2: matched null with the target and its decoys
xs = np.linspace(0, 1, 200); dens = np.exp(-((xs - 0.40) / 0.15) ** 2) + 0.25 * np.exp(-((xs - 0.62) / 0.12) ** 2)
dens /= dens.max()
bx0, bw_, by0, bh = 52.5, 33, 33, 17
px, py = bx0 + xs * bw_, by0 - dens * bh
for k in range(len(xs) - 1):                                        # density filled as a gradient towards its tail
    ax.add_patch(Polygon([(px[k], by0), (px[k], py[k]), (px[k + 1], py[k + 1]), (px[k + 1], by0)],
                         fc=K.tint(K.SLATE, 0.55 + 0.35 * abs(xs[k] - 0.40) / 0.6), ec="none", zorder=2))
ax.plot(px, py, color=K.SLATE, lw=0.7, zorder=3)
ax.plot([bx0, bx0 + bw_], [by0, by0], color=K.INK, lw=0.5, zorder=3)
xm_ = bx0 + 0.40 * bw_; xt_ = bx0 + 0.86 * bw_
ax.plot([xm_, xm_], [by0 - bh - 1, by0], color=K.INK2, lw=0.5, ls=(0, (2, 1.2)), zorder=3)
ax.plot([xt_, xt_], [by0 - bh * 0.95, by0], color=K.BLUE, lw=1.4, zorder=4)
ax.text(xt_ + 0.8, by0 - bh * 0.9, "signature", fontsize=K.T_SMALL - 0.4, color=K.BLUE, va="center")
for xd in (0.33, 0.38, 0.44, 0.47, 0.51, 0.55, 0.60, 0.69, 0.78):
    ax.plot([bx0 + xd * bw_] * 2, [by0 + 0.6, by0 + 2.4], color=K.AQUA, lw=0.8, zorder=4)
ax.text(bx0 + bw_ + 0.3, by0 + 1.5, "decoys", fontsize=K.T_SMALL - 0.6, color=K.AQUA, va="center")
ax.add_patch(FancyArrowPatch((xm_, by0 - bh - 0.6), (xt_, by0 - bh - 0.6), arrowstyle="<|-|>", mutation_scale=4.5, lw=0.6,
                             color=K.INK, zorder=5))
ax.text((xm_ + xt_) / 2, by0 - bh - 1.2, "Δ", fontsize=K.T_LAB, ha="center", va="bottom")
ax.text(bx0 - 0.5, 37.0, "1,000 random sets of the signature’s size\nand engine; Δ = statistic − null median\n(AUC or "
        "C-index; folded or directed)", fontsize=K.T_SMALL - 0.6, color=K.INK2, va="top", linespacing=1.05)

# stage 3: deltas over 12 cohorts, the signature and its decoys
dv = np.array([0.09, 0.02, 0.12, 0.06, -0.03, 0.08, 0.10, 0.04, 0.07, -0.01, 0.11, 0.05])
cx0, cdx, cy0, sc = 97, 2.75, 27, 75
ax.add_patch(Rectangle((cx0 - 1.5, cy0 - 0.13 * sc), 12 * cdx + 1.5, 0.13 * sc, fc=K.tint(K.BLUE, 0.92), ec="none", zorder=1))
ax.plot([cx0 - 1.5, cx0 + 12 * cdx], [cy0, cy0], color=K.INK, lw=0.5, zorder=2)
for i, v in enumerate(dv):
    for d_ in rng.normal(0, 0.035, 6):
        ax.add_patch(Circle((cx0 + i * cdx + rng.uniform(-.5, .5), cy0 - d_ * sc), 0.32, fc=K.AQUA, ec="none", alpha=0.8, zorder=3))
    ax.add_patch(Circle((cx0 + i * cdx, cy0 - v * sc), 0.9, fc=K.BLUE, ec="white", lw=0.35, zorder=4))
ax.text(cx0 + 12 * cdx + 0.5, cy0, "0", fontsize=K.T_SMALL - 0.4, va="center", color=K.INK2)
ax.text(cx0 - 1.5, 12.3, "Δ in each of 12 cohorts", fontsize=K.T_SMALL - 0.4, color=K.INK2, va="top")
ax.text(cx0 - 1.5, 36.5, "Wilcoxon signed-rank over cohorts:\ncalled at P < 0.05 with a positive median;\nthe same test "
        "for every decoy", fontsize=K.T_SMALL - 0.6, color=K.INK2, va="top", linespacing=1.05)

# stage 4: the estimate
sx0, sy0 = 139, 15.5
ax.text(sx0 + 0.5, 10.5, "R targets called", fontsize=K.T_SMALL - 0.4, va="center", color=K.INK)
for i in range(10):
    c = K.BLUE if i < 8 else K.tint(K.BLUE, 0.55)
    ax.add_patch(Rectangle((sx0 + 0.5 + i * 3.3, sy0 - 1.6), 2.9, 3.2, fc=c, ec="none", zorder=2))
ax.add_patch(FancyArrowPatch((sx0 + 0.5 + 8 * 3.3, sy0 + 2.6), (sx0 + 0.5 + 10 * 3.3 - 0.4, sy0 + 2.6), arrowstyle="-",
                             lw=0.6, color=K.INK, zorder=3))
ax.text(sx0 + 0.5 + 9 * 3.3 - 0.2, sy0 + 3.2, "expected false\n= Σ$_s$ $p_s$", fontsize=K.T_SMALL - 0.7, ha="center", va="top",
        linespacing=1.0)
ax.text(sx0 + 17.5, 26.0, "FDR = Σ$_s$ $p_s$ / R", fontsize=K.T_LAB + 1.2, ha="center", va="center", fontweight="bold")
ax.text(sx0 + 0.5, 31.5, "$p_s$: share of target s’s decoys called\n× $π_0$: share of null targets (corrected)\n"
        "check: outcomes permuted, every call\nfalse; estimated / realized ≈ 1", fontsize=K.T_SMALL - 0.6, color=K.INK2, va="top",
        linespacing=1.1)

# ================================================================== b  the library
SIZ = R.drop_duplicates("signature")[["signature", "n_genes", "engine"]].assign(family=lambda d: d.signature.map(famd))
LADDER = sorted(R.null_size.unique()); MEDSZ = float(SIZ.n_genes.median())
ax = C.panel("b", 31, 60, 50, 27, head="52 signatures, three engines")
for i, (e, name) in enumerate([("gsva", "GSVA"), ("wmean", "weighted mean"), ("ssgsea", "ssGSEA")]):
    g = SIZ[SIZ.engine == e].sort_values("n_genes")
    lv = np.log10(g.n_genes.to_numpy(float))
    off = K.swarm(lv, 0.34, 22)
    K.dots(ax, g.n_genes, i + off, [K.FAM[f] for f in g.family], s=11)
    ax.text(2.4, i, f"{name} ({len(g)})", ha="right", va="center", fontsize=K.T_NOTE)
K.ref(ax, MEDSZ, axis="v")
ax.text(MEDSZ * 1.08, 2.66, f"median {MEDSZ:g}", fontsize=K.T_SMALL, color=K.INK2, va="center")
for r_ in LADDER:
    ax.plot([r_, r_], [-0.95, -0.62], color=K.SLATE, lw=0.8, clip_on=False)
ax.text(2.4, -0.80, "null sizes", ha="right", va="center", fontsize=K.T_SMALL, color=K.SLATE)
ax.set_xscale("log"); ax.set_xlim(2.7, 560); ax.set_ylim(2.9, -1.05); ax.minorticks_off()
ax.set_xticks([3, 10, 30, 100, 300]); ax.set_xticklabels(["3", "10", "30", "100", "300"])
ax.set_yticks([]); K.frame(ax, left=False); ax.set_xlabel("genes in signature")
K.hlegend(C, 33, 95.5, [(K.FAM[f], K.FAMLAB[f]) for f in ("Inflamed/IFN", "Other", "Confounder")], s=8, size=K.T_SMALL)

# ================================================================== c  the null moves
ax = C.panel("c", 98, 60, 32, 27, head="The null moves with size and engine")
ECOL = {"gsva": K.INK, "ssgsea": K.AQUA, "wmean": K.ORANGE}
Q = (R.assign(p95l=R["ladder_p95"] if "ladder_p95" in R else R.rand_p95)
     .groupby(["engine", "null_size"]).agg(p95=("p95l", "median"), lo=("p95l", lambda v: np.percentile(v, 25)),
                                             hi=("p95l", lambda v: np.percentile(v, 75)), n=("p95l", "size")).reset_index())
for e in ("gsva", "ssgsea", "wmean"):
    g = Q[(Q.engine == e) & (Q.n >= 10)].sort_values("null_size")
    ax.plot(g.null_size, g.p95, color=ECOL[e], lw=1.0, zorder=3)
    K.dots(ax, g.null_size, g.p95, ECOL[e], s=8)
old95 = R.cohort.map(FIX.rand_p95).median()
ax.axhline(old95, color=K.INK2, lw=0.6, ls=(0, (2, 1.5)), zorder=2)
ax.text(520, old95 + 0.003, "fixed 48-gene null", fontsize=K.T_SMALL, color=K.INK2, va="bottom", ha="right")
ax.set_xscale("log"); ax.set_xlim(3.2, 560); ax.minorticks_off()
ax.set_xticks([4, 10, 40, 100, 400]); ax.set_xticklabels(["4", "10", "40", "100", "400"])
ax.set_ylim(.58, .77); ax.set_yticks([.6, .65, .7, .75])
ax.set_xlabel("random gene-set size"); ax.set_ylabel("null 95th percentile (AUC)")
for (e, name), yy in zip((("gsva", "GSVA"), ("ssgsea", "ssGSEA"), ("wmean", "weighted mean")), (0.98, 0.89, 0.80)):
    ax.text(0.98, yy, name, transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL, color=ECOL[e], fontweight="bold")

# ================================================================== d  fixed against matched
ax = C.panel("d", 148, 60, 30, 27, head="Matched vs fixed null")
m = R.assign(delta_old=R.AUC - R.cohort.map(FIX.rand_med))
s = m.groupby("signature")[["delta", "delta_old"]].median(); s["family"] = s.index.map(famd)
lo_, hi_ = -.025, .125
ax.fill_between([lo_, hi_], [lo_, hi_], [hi_, hi_], color=K.GRID, lw=0, zorder=0)
K.identity(ax, lo_, hi_)
for f in ("Other", "Confounder", "Inflamed/IFN"):
    g = s[s.family == f]; K.dots(ax, g.delta_old, g.delta, K.FAM[f], s=11)
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_xticks([0, .05, .1]); ax.set_yticks([0, .05, .1])
ax.set_xlabel("median ΔAUC, fixed null"); ax.set_ylabel("median ΔAUC, matched")
rr = spearmanr(s.delta_old, s.delta)[0]
ax.text(0.04, 0.96, f"ρ = {rr:.2f}", transform=ax.transAxes, fontsize=K.T_SMALL, va="top")

# ================================================================== e  co-expression per cohort
ax = C.panel("e", 22, 108, 30, 44, head="Signatures are co-expressed")
tt = RB.drop_duplicates(["target", "cohort"])
order = tt.groupby("cohort").r_target.median().sort_values().index
for i, co in enumerate(order):
    ax.axhspan(i - 0.5, i + 0.5, color=K.GRID if i % 2 == 0 else "white", lw=0, zorder=0)
    for dy, v, c in ((0.25, RB[(RB.cohort == co) & (RB["class"] == "rand")].r_decoy, K.GRAY),
                     (0.0, RB[(RB.cohort == co) & (RB["class"] == "coex")].r_decoy, K.AQUA),
                     (-0.25, tt[tt.cohort == co].r_target, K.INK)):
        K.band(ax, i + dy, v=v, color=c, orient="h", width=0.55, dot_s=5)
ax.set_yticks(range(len(order))); ax.set_yticklabels([K.COHORT[c] for c in order]); K.frame(ax, left=False)
ax.tick_params(axis="y", length=0)
ax.set_xlim(-0.12, 0.9); ax.set_xticks([0, .25, .5, .75]); ax.set_ylim(-0.6, len(order) - 0.4)
ax.set_xlabel("mean pairwise r")
K.hlegend(C, 14, 159.5, [(K.INK, "signatures"), (K.AQUA, "co-expr. decoys"), (K.GRAY, "random decoys")], size=K.T_SMALL, s=7, gap=1.4)

# ================================================================== f  one signature in every cohort
T1, CO = "Chemokine_Messina", "ICB_Gide"
ax = C.panel("f", 70, 108, 38, 44, head=f"{K.lab(T1)}, every cohort")
e1 = EX[EX.target == T1]
cords = e1[e1.kind == "target"].set_index("cohort").fold.sort_values().index
for i, co in enumerate(cords):
    g = e1[e1.cohort == co]
    bank = g[g.kind == "bank"].fold.to_numpy(); dec = g[g.kind == "decoy"].fold.to_numpy(); tv = g[g.kind == "target"].fold.iloc[0]
    z = 20 - i
    hl = co == CO
    K.density(ax, bank, K.SLATE if hl else K.GRAY, lo=0.5, hi=1.0, scale=1.6, base=i, alpha=0.55 if hl else 0.35, lw=0.5, z=z)
    ax.plot([0.5, 1.0], [i, i], color=K.RULE, lw=0.3, zorder=z - 0.5)
    ax.scatter(dec, np.full(len(dec), i + 0.12), marker="|", s=6, color=K.AQUA, lw=0.5, zorder=z + 0.3)
    ax.scatter([tv], [i + 0.12], s=14, color=K.BLUE, edgecolors="white", linewidths=0.4, zorder=z + 0.5)
    if hl:
        med = np.median(bank); tvg, medg = tv, med
        ax.annotate("", xy=(tv, i + 1.55), xytext=(med, i + 1.55), arrowprops=dict(arrowstyle="<|-|>", lw=0.5, color=K.INK,
                    mutation_scale=4), zorder=30)
        ax.plot([med, med], [i, i + 1.55], color=K.INK, lw=0.4, ls=(0, (1.5, 1)), zorder=29)
        ax.text((tv + med) / 2, i + 1.75, f"Δ = {tv - med:+.2f}".replace("-", "−"), ha="center", va="bottom", fontsize=K.T_SMALL,
                zorder=30)
ax.set_yticks(range(len(cords))); ax.set_yticklabels([K.COHORT[c] for c in cords])
for tl, co in zip(ax.get_yticklabels(), cords):
    if co == CO:
        tl.set_fontweight("bold")
ax.tick_params(axis="y", length=0); K.frame(ax, left=False)
ax.set_xlim(0.5, 1.0); ax.set_ylim(-0.4, len(cords) + 1.3); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("folded AUC for response")
K.hlegend(C, 66, 159.5, [(K.GRAY, "null bank", "s"), (K.AQUA, "decoys", "|"), (K.BLUE, "signature")], size=K.T_SMALL, s=7, gap=1.4)


# ================================================================== g  delta over cohorts, two signatures
def per_cohort(ax, T, color, ylabels):
    e = EX[EX.target == T]
    rows = []
    for co, g in e.groupby("cohort"):
        m = np.median(g[g.kind == "bank"].fold)
        rows.append((co, g[g.kind == "target"].fold.iloc[0] - m, (g[g.kind == "decoy"].fold - m).to_numpy(),
                     g[g.kind == "decoy"].unit.to_numpy()))
    rows.sort(key=lambda r: r[1])
    dmat = {}
    for i, (co, dt, dd, du) in enumerate(rows):
        ax.axhspan(i - 0.5, i + 0.5, color=K.GRID if i % 2 == 0 else "white", lw=0, zorder=0)
        ax.scatter(dd, i + K.swarm(dd, 0.36, 12), s=2.0, color=K.AQUA, alpha=0.8, lw=0, zorder=2)
        K.dots(ax, [dt], [i], color, s=15, z=5)
        for u, v in zip(du, dd):
            dmat.setdefault(u, []).append(v)
    ax.axvline(0, color=K.INK2, lw=0.5, zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([K.COHORT[r[0]] for r in rows] if ylabels else [])
    K.frame(ax, left=False); ax.set_ylim(-0.6, len(rows) - 0.4); ax.tick_params(axis="y", length=0)
    tvec = np.array([r[1] for r in rows])
    p = signed_rank_p(tvec[:, None])[0]
    dm = np.array([v for v in dmat.values() if len(v) == len(rows)]).T
    dp = signed_rank_p(dm); dcall = int(((dp < .05) & (np.median(dm, 0) > 0)).sum())
    ax.text(0.0, 1.012, f"P = {p:.3f}\ndecoys called {dcall}/{dm.shape[1]}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=K.T_SMALL - 0.3, color=K.INK, linespacing=1.0)
    ax.set_xlabel("Δ folded AUC")
    return p, dcall


ax = C.panel("g", 128, 108, 23, 44, head=None)
p1 = per_cohort(ax, "Chemokine_Messina", K.BLUE, True)
ax.set_xlim(-0.13, 0.31); ax.set_xticks([0, 0.2])
ax.set_title(K.lab("Chemokine_Messina"), fontsize=K.T_SMALL, color=K.BLUE, pad=13.0, loc="left", fontweight="bold")
ax2 = C.ax(155, 108, 23, 44)
p2 = per_cohort(ax2, "CIN70_Carter", K.ORANGE, False)
ax2.set_xlim(-0.13, 0.31); ax2.set_xticks([0, 0.2])
ax2.set_title(K.lab("CIN70_Carter"), fontsize=K.T_SMALL, color=K.ORANGE, pad=13.0, loc="left", fontweight="bold")

# ================================================================== h  ROC
ROC = pd.read_csv(f"{D}/icb/v13_fig/roc_example.csv")
ax = C.panel("h", 16, 171, 26, 22, head="ROC in Gide")
ax.fill_between(ROC.fpr, ROC.bank_p05, ROC.bank_p95, color=K.tint(K.SLATE, 0.78), lw=0, zorder=1, step="post")
dcols = [c for c in ROC.columns if c.startswith("decoy")]
dq = np.percentile(ROC[dcols].to_numpy(float), [5, 50, 95], axis=1)
ax.fill_between(ROC.fpr, dq[0], dq[2], color=K.tint(K.AQUA, 0.6), lw=0, zorder=2, step="post", alpha=0.8)
ax.step(ROC.fpr, dq[1], where="post", color=K.AQUA, lw=0.8, zorder=3)
ax.step(ROC.fpr, ROC.target, where="post", color=K.BLUE, lw=1.4, zorder=4)
K.identity(ax, 0, 1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
ax.set_xlabel("false-positive rate"); ax.set_ylabel("true-positive rate")
auc_t = np.trapezoid(ROC.target, ROC.fpr)
ax.text(0.97, 0.05, f"signature\nAUC {auc_t:.2f}", transform=ax.transAxes, ha="right", fontsize=K.T_SMALL, color=K.BLUE,
        fontweight="bold", linespacing=1.0)
ax.text(0.97, 0.30, "decoys 5–95%", transform=ax.transAxes, ha="right", fontsize=K.T_SMALL, color=K.AQUA)
ax.text(0.97, 0.40, "null 5–95%", transform=ax.transAxes, ha="right", fontsize=K.T_SMALL, color=K.SLATE)

# ================================================================== i  null width vs cohort size
ax = C.panel("i", 58, 171, 46, 22, head="The null narrows with cohort size")
ST_ = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v10_store"
n1 = {co: int(sum(json.load(open(f"{ST_}/{co}_meta.json"))["y"])) for co in R.cohort.unique()}
nn = R.groupby("cohort").n.first()
ENGC = {"gsva": K.INK, "ssgsea": K.AQUA, "wmean": K.ORANGE}
for co in nn.index:
    g = R[R.cohort == co]
    K.band(ax, nn[co], v=g.rand_sd, color=K.SLATE, orient="v", width=0.6, dot_s=6)
fr = np.median([n1[c] / nn[c] for c in nn.index])
xx = np.linspace(24, 250, 100); hn = np.sqrt(1 - 2 / np.pi)
mw = hn * np.sqrt((xx + 1) / (12 * xx * fr * xx * (1 - fr)))
ax.plot(xx, mw, color=K.ORANGE, lw=0.9, zorder=4)
ax.text(140, hn * np.sqrt(141 / (12 * 140 * fr * 140 * (1 - fr))) + 0.008, "Mann–Whitney,\nfolded", fontsize=K.T_SMALL,
        color=K.ORANGE, linespacing=1.0)
ax.set_xscale("log"); ax.set_xlim(22, 270); ax.set_xticks([30, 60, 120, 240]); ax.set_xticklabels(["30", "60", "120", "240"])
ax.minorticks_off()
ax.set_xlabel("patients in cohort"); ax.set_ylabel("null SD (AUC)")

# ================================================================== j  fixed vs matched, per cohort
ax = C.panel("j", 134, 171, 44, 22, head="A fixed null over-calls")
mm = R.assign(fix95=R.cohort.map(FIX.rand_p95))
pc = mm.groupby("cohort").apply(lambda g: pd.Series(dict(fixed=100 * (g.AUC > g.fix95).mean(), matched=100 * (g.AUC > g.rand_p95).mean())))
pc = pc.sort_values("matched")
for i, co in enumerate(pc.index):
    a_, b_ = pc.loc[co, "matched"], pc.loc[co, "fixed"]
    ax.plot([a_, b_], [i, i], color=K.tint(K.ORANGE, 0.5) if b_ > a_ else K.tint(K.SLATE, 0.6), lw=1.2, zorder=1,
            solid_capstyle="butt")
    K.dots(ax, [b_], [i], K.ORANGE, s=9, z=3); K.dots(ax, [a_], [i], K.INK, s=9, z=4)
ax.set_yticks(range(len(pc))); ax.set_yticklabels([K.COHORT[c] for c in pc.index], fontsize=K.T_SMALL - 0.4)
K.frame(ax, left=False); ax.tick_params(axis="y", length=0)
ax.set_ylim(-0.6, len(pc) - 0.4); ax.set_xlabel("tests clearing the null’s 95th pct (%)"); ax.set_xlim(0, 100)
ax.text(0.98, 0.30, "fixed 48-gene", transform=ax.transAxes, ha="right", fontsize=K.T_SMALL, color=K.ORANGE)
ax.text(0.98, 0.17, "matched", transform=ax.transAxes, ha="right", fontsize=K.T_SMALL, color=K.INK)

C.save("Figure1_GB10")
print("Figure1 j | fixed", round(100 * (mm.AUC > mm.fix95).mean(), 1), "matched", round(100 * (mm.AUC > mm.rand_p95).mean(), 1),
      "fixed higher in", int((pc.fixed > pc.matched).sum()), "| h AUC", round(auc_t, 3), "| e medians",
      tt.groupby("cohort").r_target.median().median().round(2),
      RB[RB["class"] == "coex"].groupby("cohort").r_decoy.median().median().round(2),
      RB[RB["class"] == "rand"].groupby("cohort").r_decoy.median().median().round(2),
      "| i sd range", R.rand_sd.min().round(3), R.rand_sd.max().round(3))
print("Figure1_GB10 |", "median size", MEDSZ, "| rho", round(rr, 2), "| Chemokine", p1, "| CIN70", p2,
      "| f: Gide target", round(tvg, 3), "median", round(medg, 3))
