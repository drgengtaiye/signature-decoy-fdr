"""Figure 5 redrawn at the density and in the look of the AML multi-omic figures (pixl_nc.py).
What the estimate estimates, and how it compares with realized and standard error rates.
Tables: icb/v12_truth (validation/PLAN_round7_estimand_and_bias.md), icb/v14_standard (validation/PLAN_round9_standard_fdr_comparison.md),
icb/v10_calib_* (validation/PLAN_calibration_round5.md). New panels are descriptive reads of those tables
(validation/NOTE_figures_v17_descriptive_only.md)."""
import glob
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt
from matplotlib.patches import Rectangle

I = D + "/icb/"
TR = I + "v12_truth/"
S14 = I + "v14_standard/"
P0 = pd.read_csv(TR + "null_pi0.csv")
SN = pd.read_csv(TR + "summary_null.csv").set_index("stat")
NS = pd.read_csv(I + "v10_calib_null_summary.csv").set_index(["stat", "estimator"])
PT = pd.read_csv(I + "v10_calib_null_per_target.csv", index_col=0)
SR = pd.read_csv(TR + "summary_response.csv")
ST = pd.read_csv(TR + "summary_tcga.csv")
SIGNAL = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(TR + "signal_*.csv"))])
TCGAP = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(TR + "tcga_planted_*.csv"))])
GN = pd.read_csv(S14 + "standard_globalnull.csv")
SP = pd.read_csv(S14 + "standard_planted.csv")
Q = pd.read_csv(S14 + "standard_planted_replicates.csv")
FAMT = pd.read_csv(I + "v5_signature_clusters.csv").set_index("signature")
TINFO = pd.read_csv(I + "v12_pkg/tcga_decoy_info.csv.gz")

check_palette("scenario", SCEN)

# TCGA settings keep the hue the old sheet gave them: library as P1, outside as P3, both as P2.
GROUPS = ["P1", "P2", "P3", "T1", "T2", "T3"]
HUE = {"P1": SCEN["P1"], "P2": SCEN["P2"], "P3": SCEN["P3"], "T1": SCEN["P1"], "T2": SCEN["P3"], "T3": SCEN["P2"]}
GLAB = {"P1": "3 inflamed signatures", "P2": "3 other signatures", "P3": "outside every signature",
        "T1": "2 signatures", "T2": "1 Hallmark set", "T3": "both"}
GROW = {"P1": "inflamed", "P2": "other", "P3": "outside", "T1": "library", "T2": "outside", "T3": "both"}
STAT = {"fold": "F", "dir": "D"}
DOT = (0, (1.0, 1.5))


def grp(s):
    s = s.replace("TCGA ", "")
    return s.split()[0]


def mark(ax, x, y, hue, stat="fold", strong=True, small=False, opened=False, zorder=5):
    """A planted setting: circle folded, triangle directed; small 6 cohorts; light the lower beta; open pi0-corrected."""
    m = "o" if stat == "fold" else "^"
    s = (2.3 if small else 3.5) * (1.1 if m == "^" else 1)
    if opened:
        ax.plot(x, y, marker=m, ms=s, mfc="white", mec=hue if strong else tint(hue, 0.55), mew=0.6, ls="none", zorder=zorder)
    else:
        ax.plot(x, y, marker=m, ms=s, mfc=hue if strong else tint(hue, 0.45), mec="white" if strong else hue, mew=0.3,
                ls="none", zorder=zorder)


def logmini(pg, x, y, w, h, lim, ticks, ylab=False):
    ax = pg.ax(x, y, w, h)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.plot(lim, lim, color=INK, lw=0.5, ls=DOT, zorder=1)
    for f in (ax.set_xticks, ax.set_yticks):
        f(ticks, [f"{t:g}" for t in ticks])
    ax.tick_params(which="minor", length=0)
    ax.xaxis.set_minor_formatter(plt.NullFormatter()); ax.yaxis.set_minor_formatter(plt.NullFormatter())
    if ylab:
        ax.set_ylabel("estimate (%)")
    return ax


def setting_title(pg, x, y, g):
    pg.c.plot(x + 0.7, y, "o", ms=3.0, mfc=HUE[g], mec="white", mew=0.2, zorder=4)
    pg.txt(x + 1.9, y, f"{g}, {GLAB[g]}")


def row_key(pg, x, y, items):
    """value key in millimetres: (colour, line style, label, value), one per line; the value right-aligned at x + 17.2"""
    for i, (c, ls, lab, v) in enumerate(items):
        yy = y + 2.25 * i
        pg.c.plot([x, x + 2.2], [yy, yy], color=c, lw=0.8, ls=ls, solid_capstyle="butt", zorder=6)
        pg.txt(x + 2.8, yy, lab, zorder=6)
        pg.txt(x + 17.2, yy, v, ha="right", zorder=6)


pg = Page(240)

# ---------------------------------------------------------------- a, b  false calls under permutation
for lt, st, xa, lab in (("a", "fold", 8, "folded"), ("b", "dir", 38, "directed")):
    pg.panel(xa - 7 if lt == "a" else xa - 4.5, 3, lt, f"False calls, {lab}")
    ax = pg.ax(xa, 10, 24, 22)
    x = P0[P0.stat == st]
    hist(ax, x.R[x.R <= 8], np.arange(-0.5, 9.5), weights=100 * np.ones((x.R <= 8).sum()) / len(x))
    r, e, e0, nom = x.R.mean(), x.E.mean(), x.E0.mean(), 52 * 0.025
    ax.axvline(nom, color=NOM, lw=0.7, ls=(0, (3, 1.5)), zorder=4)
    ax.axvline(e0, color=PI0, lw=0.7, ls=DOT, zorder=4)
    ax.axvline(e, color=EST, lw=0.7, zorder=4)
    ax.axvline(r, color=INK, lw=0.7, zorder=4)
    ax.set_xlim(-0.6, 8.6); ax.set_xticks([0, 2, 4, 6, 8]); ax.set_ylim(0, 62); ax.set_yticks([0, 20, 40, 60])
    ax.set_xlabel("false calls per library")
    if lt == "a":
        ax.set_ylabel("permutations (%)")
    row_key(pg, xa + 7.5, 11.2, [(INK, "-", "realized", f"{r:.2f}"), (EST, "-", "decoys", f"{e:.2f}"),
                                 (PI0, DOT, "π0-corrected", f"{e0:.2f}"), (NOM, (0, (3, 1.5)), "52 × 2.5%", f"{nom:.2f}")])
    pg.txt(xa + 24, 21.6, f"> 8: {(x.R > 8).sum()} of {len(x):,}", ha="right")
    pg.src(lt, x[["perm", "R", "E", "pi0", "E0"]].assign(nominal=nom))

# ---------------------------------------------------------------- c  estimate against false calls, per permutation
pg.panel(66, 3, "c", "Per permutation")
rng = np.random.default_rng(5)      # display jitter of the whole-number false calls only
for j, (st, m) in enumerate((("fold", "o"), ("dir", "^"))):
    ax = pg.ax(72 + 13.2 * j, 10, 12, 22)
    x = P0[P0.stat == st]
    ax.plot(x.R + rng.uniform(-0.25, 0.25, len(x)), x.E, m, ms=1.5, mfc=EST, mec="none", alpha=0.6, ls="none", rasterized=True, zorder=3)
    ax.set_xscale("symlog", linthresh=2, linscale=0.8)
    ax.plot([0, 8], [0, 8], color=INK, lw=0.5, ls=DOT, zorder=2)
    ax.set_xlim(-0.5, 32); ax.set_xticks([0, 2, 10], ["0", "2", "10"]); ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_ylim(0, 8); ax.set_yticks([0, 2, 4, 6, 8])
    if j == 0:
        ax.set_ylabel("decoy estimate")
    else:
        ax.set_yticklabels([])
    r_ = np.corrcoef(x.R, x.E)[0, 1]
    pg.txt(72 + 13.2 * j + 0.2, 8.4, f"r = {r_:.2f}")
    pg.txt(72 + 13.2 * j + 12, 11.2, STAT[st], ha="right", weight="bold")
pg.txt(72 + 12.6, 37.0, "false calls", size=6, ha="center")
pg.src("c", P0[["perm", "stat", "R", "E"]])

# ---------------------------------------------------------------- d  per signature, permuted labels
pg.panel(98.5, 3, "d", "Call rate per signature")
ax = pg.ax(105, 10, 19.5, 22)
ax.plot([0, 6], [0, 6], color=INK, lw=0.5, ls=DOT, zorder=1)
point(ax, 100 * PT.dir_decoy_callrate, 100 * PT.dir_target_callrate, SIG, "^", ms=2.3)
point(ax, 100 * PT.fold_decoy_callrate, 100 * PT.fold_target_callrate, SIG, "o", ms=2.3)
ax.set_xlim(0, 6); ax.set_ylim(0, 6); ax.set_xticks([0, 2, 4, 6]); ax.set_yticks([0, 2, 4, 6])
ax.set_xlabel("decoys called (%)"); ax.set_ylabel("signature called (%)")
rf, rd = NS.loc[("fold", "coex A"), "per_target_r"], NS.loc[("dir", "coex A"), "per_target_r"]
key(pg, 110.2, 28.2, [(SIG, f"folded, r = {rf:.2f}", "o")], sw=1.2)
key(pg, 110.2, 30.6, [(SIG, f"directed, r = {rd:.2f}", "^")], sw=1.2)
pg.src("d", PT.reset_index(names="signature").assign(r_fold=rf, r_dir=rd))

# ---------------------------------------------------------------- e  estimated over realized, every analysis
pg.panel(128, 3, "e", "Estimated over realized false calls")
ax = pg.ax(148.5, 11.5, 15.5, 20.5)
ALAB = {("response", "fold"): "response, folded", ("response", "dir"): "response, directed", ("tcga_OS", "fold"): "TCGA survival",
        ("breast_pc1", "fold"): "breast, PC1", ("hallmark_OS", "fold"): "Hallmark"}
ax.axvspan(0.8, 1.25, color=EMPTY, lw=0, zorder=0)
zero(ax, 1.0)
for i, r in enumerate(GN.itertuples()):
    yy = len(GN) - 1 - i
    ci(ax, r.decoy_lo, r.decoy_hi, yy + 0.17, EST, lw=0.7)
    ci(ax, r.nominal_lo, r.nominal_hi, yy - 0.17, NOM, lw=0.7)
    point(ax, r.decoy_ratio, yy + 0.17, EST, "o", ms=3.0, zorder=5)
    point(ax, r.nominal_ratio, yy - 0.17, NOM, "D", ms=2.3, filled=False, zorder=5)
    ax.text(1.03, yy, f"{r.decoy_ratio:.2f}", transform=ax.get_yaxis_transform(), va="center", fontsize=5)
    ax.text(1.42, yy, f"{r.nominal_ratio:.2f}", transform=ax.get_yaxis_transform(), va="center", fontsize=5)
ax.text(1.03, 4.72, "dec.", transform=ax.get_yaxis_transform(), va="bottom", fontsize=5, weight="bold")
ax.text(1.42, 4.72, "nom.", transform=ax.get_yaxis_transform(), va="bottom", fontsize=5, weight="bold")
ax.set_xscale("log"); ax.set_xlim(0.38, 2.1); ax.set_xticks([0.5, 1, 2], ["0.5", "1", "2"]); ax.xaxis.set_minor_locator(plt.NullLocator())
ax.set_yticks(range(len(GN))[::-1], [ALAB[(a, s)] for a, s in zip(GN.analysis, GN.stat)], fontsize=5)
ax.set_ylim(-0.6, 4.6); ax.set_xlabel("estimated / realized")
kx = key(pg, 128.5, 8.2, [(EST, "decoys", "o"), (NOM, "nominal", "open-D"), (INK, "95% interval", "line")], gap=1.6, sw=1.6)
key(pg, kx, 8.2, [("#E2E2E2", "0.8–1.25")], sw=1.6)
pg.src("e", GN)

# ---------------------------------------------------------------- f  planted, truth by gene content
LIM_G = {"P1": ((10, 50), [10, 20, 50]), "P2": ((30, 100), [30, 50, 100]), "P3": ((50, 300), [50, 100, 200])}
LIM_F = {"P1": ((0.1, 30), [0.1, 1, 10]), "P2": ((5, 100), [10, 30, 100]), "P3": ((20, 300), [30, 100, 300])}
LIM_T = {"T1": ((1, 50), [1, 3, 10, 30]), "T2": ((5, 200), [10, 30, 100]), "T3": ((5, 100), [10, 30, 100])}
RESP = SP[SP.data == "response"].assign(g=lambda d: d.scenario.map(grp))
RY, RH, RW, RS = 47.5, 21.0, 21.5, 27.5

pg.panel(1, 41, "f", "Planted signal, truth by gene content")
for j, g in enumerate(("P1", "P2", "P3")):
    xa = 10 + RS * j
    setting_title(pg, xa - 1.5, 45.2, g)
    ax = logmini(pg, xa, RY, RW, RH, *LIM_G[g], ylab=j == 0)
    for r in RESP[RESP.g == g].itertuples():
        mark(ax, [r.realized_gc], [r.decoy], HUE[g], r.stat, r.beta == 1, r.cohorts == 6)
    if j == 1:
        ax.set_xlabel("realized FDR (%)")
pg.src("f", RESP[["setting", "stat", "cohorts", "beta", "realized_gc", "decoy"]].rename(columns={"decoy": "estimate"}))

# ---------------------------------------------------------------- g  planted, formal truth
pg.panel(92, 41, "g", "Planted signal, formal truth")
for j, g in enumerate(("P1", "P2", "P3")):
    xa = 101 + RS * j
    setting_title(pg, xa - 1.5, 45.2, g)
    ax = logmini(pg, xa, RY, RW, RH, *LIM_F[g], ylab=j == 0)
    for r in RESP[RESP.g == g].itertuples():
        ax.plot([r.realized, r.realized], [r.decoy_pi0, r.decoy], color=OTHER, lw=0.4, zorder=2)
        mark(ax, [r.realized], [r.decoy_pi0], HUE[g], r.stat, r.beta == 1, r.cohorts == 6, opened=True)
        mark(ax, [r.realized], [r.decoy], HUE[g], r.stat, r.beta == 1, r.cohorts == 6)
    if j == 1:
        ax.set_xlabel("realized FDR (%)")
pg.src("g", RESP[["setting", "stat", "cohorts", "beta", "realized", "decoy", "decoy_pi0"]].rename(columns={"decoy": "estimate", "decoy_pi0": "estimate_pi0"}))

# the marks of f, g, h, in one key
yk = 77.0
kx = key(pg, 10, yk, [(INK, "folded", "o"), (INK, "directed", "^"), (INK, "π0-corrected", "open-o")], gap=2.2)
pg.c.plot(kx + 0.9, yk, "o", ms=2.3, mfc=INK, mec="white", mew=0.2, zorder=3)
pg.txt(kx + 2.4, yk, "6 cohorts")
kx += 11.6
pg.c.plot(kx + 0.9, yk, "o", ms=3.5, mfc=tint(INK, 0.45), mec=INK, mew=0.3, zorder=3)
pg.txt(kx + 2.4, yk, "β = 0.5 (TCGA 0.3)")
kx += 21.5
pg.c.plot(kx + 0.9, yk, "o", ms=3.5, mfc=INK, mec="white", mew=0.3, zorder=3)
pg.txt(kx + 2.4, yk, "β = 1 (TCGA 0.6)")
kx += 19.5
pg.line(kx, yk, kx + 2.0, yk, color=INK, lw=0.5, ls=DOT)
pg.txt(kx + 2.6, yk, "identity; axes differ between panels")

# ---------------------------------------------------------------- h  planted, TCGA
TC = SP[SP.data == "tcga"].assign(g=lambda d: d.scenario.map(grp))
pg.panel(1, 81, "h", "TCGA survival, new causal sets")
c_ = TINFO[(TINFO["class"] == "coex") & (TINFO.set == 1)]
pg.txt(86.8, 85.2 - 4.0, f"co-expression, decoys {c_.r_achieved_med.median():.2f}, signatures {c_.r_target_med.median():.2f}", ha="right")
for j, g in enumerate(("T1", "T2", "T3")):
    xa = 10 + RS * j
    setting_title(pg, xa - 1.5, 85.2, g)
    ax = logmini(pg, xa, 87.2, RW, RH, *LIM_T[g], ylab=j == 0)
    for r in TC[TC.g == g].itertuples():
        ax.plot([r.realized, r.realized], [r.decoy_pi0, r.decoy], color=OTHER, lw=0.4, zorder=2)
        mark(ax, [r.realized], [r.decoy_pi0], HUE[g], "fold", r.beta == 0.6, False, opened=True)
        mark(ax, [r.realized], [r.decoy], HUE[g], "fold", r.beta == 0.6, False)
    if j == 1:
        ax.set_xlabel("realized FDR (%)")
pg.src("h", TC[["setting", "beta", "realized", "decoy", "decoy_pi0"]].rename(columns={"decoy": "estimate", "decoy_pi0": "estimate_pi0"}).assign(
    r_decoys=c_.r_achieved_med.median(), r_signatures=c_.r_target_med.median()))

# ---------------------------------------------------------------- i  signature minus its decoys, by family
pg.panel(92, 81, "i", "Signature minus decoys")
ax = pg.ax(107, 85.4, 20, 22.8)
fam = FAMT.family.reindex(PT.index)
rows_ = [("Inflamed/IFN", "fold"), ("Inflamed/IFN", "dir"), ("Confounder", "fold"), ("Confounder", "dir"), ("Other", "fold"), ("Other", "dir")]
ox = []
for i, (f, st) in enumerate(rows_):
    yy = len(rows_) - 1 - i
    v = 100 * (PT[f"{st}_target_callrate"] - PT[f"{st}_decoy_callrate"])[fam == f]
    strip(ax, yy, v, FAM[f], m="o" if st == "fold" else "^", ms=1.9, spread=0.27, horizontal=True)
    ax.text(1.02, yy, str(len(v)), transform=ax.get_yaxis_transform(), va="center", fontsize=5)
    ox += [dict(signature=s, family=f, stat=st, minus=val) for s, val in v.items()]
ax.text(1.02, 5.62, "n", transform=ax.get_yaxis_transform(), va="bottom", fontsize=5, weight="bold")
zero(ax, 0)
FL = {"Inflamed/IFN": "inflamed", "Confounder": "confounder", "Other": "other"}
ax.set_yticks(range(6)[::-1], [f"{FL[f]}, {STAT[s]}" for f, s in rows_], fontsize=5)
ax.set_ylim(-0.6, 5.6); ax.set_xlim(-2.5, 1.8); ax.set_xticks([-2, 0, 1], ["−2", "0", "1"]); ax.set_xlabel("difference (points)")
pg.src("i", pd.DataFrame(ox))

# ---------------------------------------------------------------- j  estimated over realized under permutation, by decoy class
pg.panel(134, 81, "j", "Null ratio by decoy class")
ax = pg.ax(158, 85.4, 14.5, 22.8)
CL = [("coex A", "co-expressed", COEX, "s", True), ("coex A, closest", "closest-matched", COEX, "s", True), ("coex B", "no TME genes", COEX, "s", True),
      ("random A", "random genes", RAND, "D", True), ("pi0", "π0-corrected", COEX, "s", False)]
ax.axvspan(0.8, 1.25, color=EMPTY, lw=0, zorder=0)
zero(ax, 1.0)
lab, oj = [], []
for si, st in enumerate(("fold", "dir")):
    for ci_, (k, name, c, m, fl) in enumerate(CL):
        yy = 11 - ci_ - (si * 6.0) - (0 if si == 0 else 0.4)
        if k == "pi0":
            v, lo, hi = SN.loc[st, "ratio_pi0"], SN.loc[st, "ratio_pi0_lo"], SN.loc[st, "ratio_pi0_hi"]
        else:
            r = NS.loc[(st, k)]
            v, lo, hi = r.rho, r.rho_lo, r.rho_hi
        ci(ax, lo, hi, yy, c, lw=0.6)
        point(ax, v, yy, c, m, ms=2.2, filled=fl, zorder=5)
        ax.text(1.04, yy, f"{v:.2f}", transform=ax.get_yaxis_transform(), va="center", fontsize=5)
        lab.append((yy, name))
        oj.append(dict(stat=st, decoys=name, ratio=v, lo=lo, hi=hi))
    pg.txt(140.4, 85.4 + 22.8 * (11.55 - (11 - 2 - si * 6.4)) / 12.3, STAT[st], weight="bold")
ax.set_yticks([l[0] for l in lab], [l[1] for l in lab], fontsize=5)
ax.tick_params(axis="y", length=0, pad=1)
ax.set_ylim(-0.75, 11.55); ax.set_xscale("log"); ax.set_xlim(0.66, 1.5)
ax.set_xticks([0.8, 1, 1.25], ["0.8", "1", "1.25"]); ax.xaxis.set_minor_locator(plt.NullLocator())
ax.set_xlabel("estimated / realized")
pg.src("j", pd.DataFrame(oj))

# ---------------------------------------------------------------- k, l  the 22 planted settings, on shared rows
SP["g"] = SP.scenario.map(grp)
SP["o"] = SP.g.map({g: i for i, g in enumerate(GROUPS)})
SP["co"] = (SP.cohorts == 6).astype(int)
SP["so"] = SP.stat.map({"fold": 0, "dir": 1})
S22 = SP.sort_values(["o", "co", "beta", "so"]).reset_index(drop=True)
assert len(S22) == 22
Y0, RH_, GAP = 128.6, 2.38, 0.9
ys, y_ = [], Y0
for i, r in S22.iterrows():
    if i and S22.g[i] != S22.g[i - 1]:
        y_ += GAP
    ys.append(y_); y_ += RH_
Y_END = y_
S22["y"] = ys

pg.panel(1, 119, "k", "Planted settings, FDR (%)")
COLS = [("R", "mean", "calls", None), ("realized_gc", "", "gene", "fdr"), ("realized", "", "formal", "fdr"), ("decoy", "", "plain", "fdr"),
        ("decoy_pi0", "", "π0", "fdr"), ("nominal", "", "nominal", "fdr"), ("permutation", "", "perm.", "fdr")]
x0, cw = 23.4, 6.75


def fdr_c(v):
    return PURPLE(np.clip((np.log10(v) + 1) / (np.log10(300) + 1), 0, 1))


def num(v):
    return f"{v:.1f}" if v < 9.95 else f"{v:.0f}"


for j, (c, h1, h2, kind) in enumerate(COLS):
    if h1:
        pg.txt(x0 + cw * j + cw / 2 - 0.1, 123.8, h1, ha="center", weight="bold")
    pg.txt(x0 + cw * j + cw / 2 - 0.1, 126.6, h2, ha="center", weight="bold")
pg.txt(x0 + cw * 2 - 0.1, 123.8, "realized", ha="center", weight="bold")
pg.txt(x0 + cw * 4 - 0.1, 123.8, "decoys", ha="center", weight="bold")
pg.txt(x0 + cw * 6 - 0.1, 123.8, "standard", ha="center", weight="bold")
for a_, b_ in ((1, 3), (3, 5), (5, 7)):
    pg.line(x0 + cw * a_ + 0.3, 125.2, x0 + cw * b_ - 0.5, 125.2, lw=0.3)
pg.txt(x0 - 0.8, 126.6, "β, cohorts", ha="right", weight="bold")
for i, r in S22.iterrows():
    pg.txt(x0 - 0.8, r.y + RH_ / 2 - 0.05, f"{r.beta:g}, {r.cohorts}, {STAT[r.stat]}", ha="right")
    for j, (c, _, _, kind) in enumerate(COLS):
        v = r[c]
        if not np.isfinite(v):
            tile(pg, x0 + cw * j, r.y, cw, RH_, EMPTY)
        elif kind == "fdr":
            tile(pg, x0 + cw * j, r.y, cw, RH_, fdr_c(v), num(v))
        else:
            tile(pg, x0 + cw * j, r.y, cw, RH_, BLUE(v / 45), f"{v:.1f}")
for g in GROUPS:
    s = S22[S22.g == g]
    pg.box(11.2, s.y.min(), 0.8, s.y.max() + RH_ - s.y.min() - 0.15, HUE[g], z=2)
    ym = (s.y.min() + s.y.max() + RH_) / 2
    pg.txt(1, ym - 1.0, g, weight="bold")
    pg.txt(1, ym + 1.0, GROW[g])
pg.cbar(x0 + 8, Y_END + 1.6, 22, 1.4, PURPLE, -1, np.log10(300), [-1, 0, 1, 2], "FDR (%)", labels=["0.1", "1", "10", "100"])
pg.src("k", S22[["setting", "stat", "reps", "R", "realized_gc", "realized", "decoy", "decoy_pi0", "nominal", "permutation"]])

pg.panel(73, 119, "l", "Parts of the plain estimate")
ax = pg.ax(74.5, Y0, 41, Y_END - Y0)
ax.set_ylim(Y_END, Y0); ax.set_xlim(0, 245)
SR_ = SR.set_index(["scenario", "beta", "cohorts", "stat"])
ST_ = ST.set_index(["kind", "beta"])
ol = []
for r in S22.itertuples():
    if r.data == "response":
        p = SR_.loc[(r.scenario, r.beta, r.cohorts, r.stat)]
    else:
        p = ST_.loc[(r.scenario.replace("TCGA ", ""), r.beta)]
    yy = r.y + RH_ / 2
    left = 0
    for v, col in ((p.exp_V, INK), (p.excess_null, HUE[r.g]), (p.part_true, tint(HUE[r.g], 0.4))):
        ax.barh(yy, v, 1.5, left=left, color=col, edgecolor="white", linewidth=0.25, zorder=3)
        left += v
    point(ax, p.estimate_pi0, yy, INK, "D", ms=2.2, filled=False, mew=0.5, zorder=6)
    ol.append(dict(setting=r.setting, stat=r.stat, false_calls_null=p.exp_V, decoy_excess=p.excess_null, decoy_rates_true=p.part_true, estimate_pi0=p.estimate_pi0))
ax.set_yticks([]); ax.set_xticks([0, 50, 100, 150, 200]); ax.set_xlabel("% of calls (expected)")
ex = SCEN["P3"]
for i, (c, lab_, k) in enumerate(((INK, "false calls of null signatures", "b"), (ex, "decoys' excess over null signatures", "b"),
                                  (tint(ex, 0.4), "decoy rates of true signatures", "b"), (INK, "π0-corrected estimate", "d"))):
    yy = Y0 + 1.4 + 2.5 * i
    if k == "b":
        pg.box(83.0, yy - 0.75, 1.8, 1.5, c, z=2)
    else:
        pg.c.plot(83.9, yy, "D", ms=2.2, mfc="white", mec=c, mew=0.5, zorder=3)
    pg.txt(85.6, yy, lab_)
pg.src("l", pd.DataFrame(ol))

# ---------------------------------------------------------------- m  estimate over realized, four estimates in the 22 settings
pg.panel(119, 119, "m", "Standard estimates beside decoys, 22 settings")
ax = pg.ax(128, 124.5, 50, 26.5)
ESTS = [("nominal", "nominal (BH)"), ("permutation", "permutation (SAM)"), ("decoy", "decoys"), ("decoy_pi0", "decoys, π0-corrected")]
ax.axhspan(0.25, 1, color=EMPTY, lw=0, zorder=0)
zero(ax, 1, axis="y")
LANE = {"P1": -0.34, "P2": -0.22, "P3": -0.10, "T1": 0.12, "T3": 0.24, "T2": 0.36}
om = []
for xi, (col, _) in enumerate(ESTS):
    s = S22[S22[col].notna()]
    ratio = s[col] / s.realized
    for r, q in zip(s.itertuples(), ratio):
        off = 0 if r.data == "tcga" else (-0.025 if r.stat == "fold" else 0.025)
        strong = r.beta in (1.0, 0.6)
        mark(ax, [xi + LANE[r.g] + off], [q], HUE[r.g], r.stat, strong, r.cohorts == 6, zorder=4)
        om.append(dict(setting=r.setting, stat=r.stat, estimate=col, ratio=q))
    nb = int((ratio < 1).sum())
    ax.text(xi, 880, f"{nb} of {len(s)} below", ha="center", va="top", fontsize=5)
ax.set_yscale("log"); ax.set_ylim(0.25, 900); ax.set_yticks([0.3, 1, 10, 100], ["0.3", "1", "10", "100"]); ax.yaxis.set_minor_locator(plt.NullLocator())
ax.set_xlim(-0.5, 3.5); ax.set_ylabel("estimate / realized FDR")
ticks, tl = [], []
for xi, (_, name) in enumerate(ESTS):
    ticks += [xi - 0.22, xi + 0.24]; tl += ["resp.", "TCGA"]
ax.set_xticks(ticks, tl, fontsize=5)
ax.tick_params(axis="x", length=0, pad=1)
for xi, (_, name) in enumerate(ESTS):
    a, b = (name.split(" (") if " (" in name else name.split(", ") if ", " in name else (name, ""))
    b = f"({b}" if " (" in name else b
    ax.text(xi, -0.135, a, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5)
    if b:
        ax.text(xi, -0.215, b, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5)
ax.text(3.47, 0.28, "estimate too low", ha="right", va="bottom", fontsize=5)
pg.src("m", pd.DataFrame(om))

# ---------------------------------------------------------------- n  replicates at or above the realized proportion
pg.panel(119, 162.2, "n", "Replicates at or above realized (%)")
Q["g"] = Q.scenario.map(grp)
Q["lab"] = np.where(Q.data == "tcga", Q.g, Q.g + " " + Q.stat.map(STAT))
LABS = ["P1 F", "P1 D", "P2 F", "P2 D", "P3 F", "P3 D", "T1", "T2", "T3"]
x0n, cwn, y0n, chn = 129.5, 10.6, 169.6, 2.38
HEAD = [("nominal", "nominal"), ("permutation", "perm."), ("decoy", "decoys"), ("decoy_pi0", "π0-corr.")]
for j, (_, h) in enumerate(HEAD):
    pg.txt(x0n + cwn * j + cwn / 2, y0n - 1.4, h, ha="center", weight="bold")
pg.txt(x0n + cwn * 4 + 5.9, y0n - 1.4, "n", ha="right", weight="bold")
on = []
for i, l in enumerate(LABS):
    s = Q[Q.lab == l]
    yy = y0n + chn * i
    pg.txt(x0n - 0.8, yy + chn / 2 - 0.05, l, ha="right")
    for j, (c, _) in enumerate(HEAD):
        y = s[s[c].notna()]
        v = 100 * np.mean(y[c] >= y.realized - 1e-9)
        tile(pg, x0n + cwn * j, yy, cwn, chn, BLUE(np.clip((v - 50) / 50, 0, 1)), f"{v:.0f}")
        on.append(dict(group=l, estimate=c, replicates=len(y), at_or_above_pct=v))
    pg.txt(x0n + cwn * 4 + 5.9, yy + chn / 2 - 0.05, f"{len(s)}", ha="right")
pg.src("n", pd.DataFrame(on))

# ---------------------------------------------------------------- o, p, q  replicate by replicate, and its two margins
pg.panel(1, 196.5, "o", "Replicate by replicate")
ax = pg.ax(9, 202, 28, 30)
sub = SIGNAL[(SIGNAL.stat == "fold") & (SIGNAL.beta == 1.0) & (SIGNAL.cohorts == 12) & (SIGNAL.R > 0)].copy()
sub["est"] = 100 * sub.E / sub.R
sub["real"] = 100 * sub.V_cp / sub.R
sub["g"] = sub.scenario.map(grp)
FL_, TOP = 0.5, 2000.0     # log axes; a realized FDR of 0 is drawn at 0.5%
lgc = lambda v: np.log10(np.clip(v, FL_, TOP))  # noqa: E731
rng = np.random.default_rng(2)
for g in ("P1", "P2", "P3"):
    s = sub[sub.g == g]
    jit = np.where(s.real <= 0, rng.uniform(-0.06, 0.06, len(s)), 0)
    ax.plot(lgc(s.est), lgc(s.real) + jit, "o", ms=1.7, mfc=HUE[g], mec="none", alpha=0.75, ls="none", zorder=3)
ax.plot([lgc(FL_), lgc(TOP)], [lgc(FL_), lgc(TOP)], color=INK, lw=0.5, ls=DOT, zorder=1)
ax.axhline(lgc(FL_) + 0.15, color=OTHER, lw=0.3, zorder=1)
tk = np.array([1, 10, 100, 1000])
ax.set_xlim(lgc(FL_) - 0.1, lgc(TOP)); ax.set_ylim(lgc(FL_) - 0.15, 2.2)
ax.set_xticks(np.log10(tk), ["1", "10", "100", "1,000"])
ax.set_yticks(np.r_[lgc(FL_), np.log10(tk[:3])], ["0", "1", "10", "100"])
ax.set_xlabel("estimated FDR (%)"); ax.set_ylabel("realized FDR (%)")
above = 100 * (sub.est >= sub.real).mean()
ax.text(lgc(FL_) - 0.03, 2.2, f"estimate ≥ realized,\n{above:.0f}% of {len(sub)}", ha="left", va="top", fontsize=5, linespacing=1.15)
key(pg, 11.2, 215, [(HUE["P1"], "P1", "o")], sw=1.2); key(pg, 11.2, 217.4, [(HUE["P2"], "P2", "o")], sw=1.2); key(pg, 11.2, 219.8, [(HUE["P3"], "P3", "o")], sw=1.2)
pg.src("o", sub[["scenario", "rep", "R", "est", "real"]].rename(columns={"est": "estimated_fdr", "real": "realized_fdr"}))


def margin(ax, v, g, bins, top):
    all_ = np.asarray(v)
    hist(ax, all_, bins)
    for k in ("P1", "P2", "P3"):
        h, _ = np.histogram(all_[np.asarray(g) == k], bins)
        ax.step(bins, np.r_[h, h[-1]], where="post", color=HUE[k], lw=0.8, zorder=4)
    ax.set_ylim(0, top)


pg.panel(41, 196.5, "p", "Estimated FDR, replicates")
ax = pg.ax(48, 202, 26, 30)
b = np.arange(1.0, 3.3001, 0.125)
margin(ax, lgc(sub.est), sub.g, b, 104)
ax.set_xlim(1.0, 3.3); ax.set_xticks([1, 2, 3], ["10", "100", "1,000"])
ax.set_xlabel("estimated FDR (%)"); ax.set_ylabel("replicates")
mE = sub.groupby("g").est.median()
ax.text(3.3, 104, "median", ha="right", va="top", fontsize=5, weight="bold")
for i, g in enumerate(("P1", "P2", "P3")):
    ax.text(3.3, 104 - 8.6 * (i + 1), f"{g}, {mE[g]:.0f}", ha="right", va="top", fontsize=5)
pg.src("p", sub[["scenario", "rep", "est"]].rename(columns={"est": "estimated_fdr"}))

pg.panel(78, 196.5, "q", "Realized FDR, replicates")
ax = pg.ax(85, 202, 26, 30)
b = np.r_[lgc(FL_) - 0.09, lgc(FL_) + 0.09, np.arange(0.375, 2.0001 + 0.125, 0.125)]
vq = lgc(sub.real)
margin(ax, vq, sub.g, b, 260)
ax.set_xlim(lgc(FL_) - 0.15, 2.2); ax.set_xticks(np.r_[lgc(FL_), 1, 2], ["0", "10", "100"])
ax.set_xlabel("realized FDR (%)"); ax.set_ylabel("replicates")
mR = sub.groupby("g").real.median()
ax.text(2.2, 260, "median", ha="right", va="top", fontsize=5, weight="bold")
for i, g in enumerate(("P1", "P2", "P3")):
    ax.text(2.2, 260 - 21.5 * (i + 1), f"{g}, {mR[g]:.0f}", ha="right", va="top", fontsize=5)
pg.src("q", sub[["scenario", "rep", "real"]].rename(columns={"real": "realized_fdr"}))

# ---------------------------------------------------------------- r  estimate minus realized, every replicate
pg.panel(115, 196.5, "r", "Estimate minus realized, every replicate")
ax = pg.ax(123, 202, 55, 30)
ALLR = SIGNAL[SIGNAL.R > 0].copy()
ALLR["gap"] = 100 * (ALLR.E - ALLR.V_cp) / ALLR.R
ALLR["g"] = ALLR.scenario.map(grp)
ALLT = TCGAP[TCGAP.R > 0].copy()
ALLT["gap"] = 100 * (ALLT.E - ALLT.V_cp) / ALLT.R
ALLT["g"] = ALLT.kind.map(grp)
GR = [(f"{g} {STAT[st]}", ALLR[(ALLR.g == g) & (ALLR.stat == st)].gap, HUE[g], "o" if st == "fold" else "^") for g in ("P1", "P2", "P3") for st in ("fold", "dir")]
GR += [(g, ALLT[ALLT.g == g].gap, HUE[g], "o") for g in ("T1", "T2", "T3")]
XP = [0, 1, 2.3, 3.3, 4.6, 5.6, 6.9, 7.9, 8.9]
rng = np.random.default_rng(7)      # display jitter only
ax.set_yscale("symlog", linthresh=100, linscale=1.15)
TOPR = 1750
or_ = {}
for xp, (l, v, c, m) in zip(XP, GR):
    v = np.asarray(v, float)
    ax.plot(xp + rng.uniform(-0.3, 0.3, len(v)), v, m, ms=1.25, mfc=c, mec="none", alpha=0.55, ls="none", rasterized=True, zorder=3)
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    ax.add_patch(Rectangle((xp - 0.33, q1), 0.66, q3 - q1, fc="none", ec=INK, lw=0.5, zorder=5))
    ax.plot([xp - 0.33, xp + 0.33], [q2, q2], color=INK, lw=1.0, solid_capstyle="butt", zorder=6)
    ax.text(xp, TOPR * 1.25, f"{100 * (v >= 0).mean():.0f}%", ha="center", va="bottom", fontsize=5)
    or_[l] = dict(group_replicates=len(v), group_share_at_or_above_pct=100 * (v >= 0).mean(), group_q25=q1, group_median=q2, group_q75=q3)
zero(ax, 0, axis="y")
ax.set_xlim(-0.6, 9.5); ax.set_xticks(XP, [g[0] for g in GR], fontsize=5.5); ax.tick_params(axis="x", length=0)
ax.set_ylim(-100, TOPR); ax.set_yticks([-50, 0, 50, 100, 500, 1000], ["−50", "0", "50", "100", "500", "1,000"]); ax.yaxis.set_minor_locator(plt.NullLocator())
ax.set_ylabel("percentage points")
rr = pd.concat([ALLR.assign(group=ALLR.g + " " + ALLR.stat.map(STAT), data="response")[["data", "group", "scenario", "beta", "cohorts", "rep", "stat", "gap"]],
                ALLT.assign(group=ALLT.g, data="tcga", stat="fold", cohorts=20).rename(columns={"kind": "scenario"})[["data", "group", "scenario", "beta", "cohorts", "rep", "stat", "gap"]]])
pg.src("r", rr.join(pd.DataFrame(or_).T, on="group"))

pg.save("Fig5")
