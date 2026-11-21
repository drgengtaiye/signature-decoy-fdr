"""Figure 5: the same procedure for prognosis in 20 TCGA cancer types: folded C-index, ten decoy sets (style: pixl_nc.py).
Parallel to Figure 6 (one library in TCGA, one in breast cohorts).

Data: icb/v12_pkg/tcga_* (applications/v12_rerun.R), v11_pkg/tcga_decoy_r_by_type.csv
(set 1 decoys), calib_v3_units.csv (response), v13_fig/km_tcga.csv (scripts/figprep_km.R), v5_signature_clusters.csv.
Panels outside the figure plan are descriptive reads of those tables (validation/NOTE_figures_v17_descriptive_only.md)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt, fm, TextPath, D

I = D + "/icb/"
P = I + "v12_pkg/"
FAMILY = pd.read_csv(I + "v5_signature_clusters.csv").set_index("signature").family.to_dict()
FLAB = {"Inflamed/IFN": "inflamed/IFN", "Confounder": "confounder", "Other": "other"}
EP = ["OS", "DSS", "PFI"]
F = {e: pd.read_csv(P + f"tcga_{e}_fdr.csv") for e in EP}
B = {e: pd.read_csv(P + f"tcga_{e}_fdr_by_set.csv") for e in EP}
CL = {e: pd.read_csv(P + f"tcga_{e}_calls.csv.gz") for e in EP}
CL = {e: c[c.stat == "fold"] for e, c in CL.items()}
INFO = pd.read_csv(P + "tcga_decoy_info.csv.gz").drop(columns="genes")
RP = pd.read_csv(I + "v11_pkg/tcga_decoy_r_by_type.csv")
GN = pd.read_csv(P + "tcga_OS_globalnull_summary.csv")
GP = pd.read_csv(P + "tcga_OS_globalnull_per_perm.csv")
GT = pd.read_csv(P + "tcga_OS_globalnull_per_target.csv")
U0 = pd.read_csv(I + "calib_v3_units.csv")
U0 = U0[(U0["null"] == "N0") & (U0.kind == "target") & (U0.stat == "fold")]
UO = pd.read_csv(P + "tcga_OS_units.csv.gz", usecols=["cohort", "n", "unit", "kind", "target", "class", "set", "fold_obs", "fold_med", "fold_p95"])
UO["dC"] = UO.fold_obs - UO.fold_med
UT = UO[UO.kind == "target"].drop_duplicates(["cohort", "target"])
KMD = pd.read_csv(I + "v13_fig/km_tcga.csv")
KMD["years"] = KMD.time / 365.25
EXS = "CIN70_Carter"

SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)",
         "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)", "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}


def lab(s):
    if s in SHORT:
        return SHORT[s]
    if "_" in s:
        a, b = s.rsplit("_", 1)
        return f"{a.replace('_', ' ')} ({b})"
    return s


def row(df, decoys="coex", boot="clusters", st="fold"):
    return df[(df.stat == st) & (df.decoys == decoys) & (df.boot == boot)].iloc[0]


def swarm(v, bin_w, step, limit):
    """offsets that lay points of similar value side by side: bins of bin_w, neighbours step apart, at most limit from the centre"""
    v = np.asarray(v, float)
    off = np.zeros(len(v))
    idx = np.floor(v / bin_w).astype(int)
    for b in np.unique(idx):
        k = np.where(idx == b)[0]
        k = k[np.argsort(v[k])]
        off[k] = np.array([0] + [s for i in range(1, len(k)) for s in (i, -i)])[:len(k)] * step
    return np.clip(off, -limit, limit)


def strip2(ax, k, vals, c, m="o", ms=1.8, bin_w=0.01, step=0.04, limit=0.3, horizontal=False, mew=0.12, med_half=0.3):
    """points of one group, side by side, with a black median"""
    vals = np.asarray(vals, float)
    o = k + swarm(vals, bin_w, step, limit)
    a, b = (vals, o) if horizontal else (o, vals)
    point(ax, a, b, c, m=m, ms=ms, mew=mew)
    q = float(np.median(vals))
    (xs, ys) = (([q, q], [k - med_half, k + med_half]) if horizontal else ([k - med_half, k + med_half], [q, q]))
    ax.plot(xs, ys, color=INK, lw=0.7, zorder=5, solid_capstyle="butt")


def half_up(v, nd=2):
    return f"{np.floor(v * 10 ** nd + 0.5) / 10 ** nd:.{nd}f}"


def regions(ax, x, y, color, xlim, ylim, levels=(0.5, 0.75, 0.9, 0.97), n=120, bw=0.15, z=2):
    """regions holding the stated shares of the density: pale fills of one hue and a hairline at each level"""
    from scipy.stats import gaussian_kde
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = gaussian_kde(np.vstack([x, y]), bw_method=bw)
    gx, gy = np.meshgrid(np.linspace(*xlim, n), np.linspace(*ylim, n))
    dz = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    s = np.sort(dz.ravel())[::-1]
    cs = np.cumsum(s) / s.sum()
    th = sorted(s[np.searchsorted(cs, lv)] for lv in levels)
    ax.contourf(gx, gy, dz, levels=list(th) + [dz.max() * 1.01], colors=[tint(color, f) for f in (0.10, 0.17, 0.26, 0.38)], zorder=z)
    ax.contour(gx, gy, dz, levels=th, colors=[color], linewidths=0.3, zorder=z + 0.1)


def km_curve(time, event):
    t, e = np.asarray(time, float), np.asarray(event, int)
    s, xs, ys = 1.0, [0.0], [1.0]
    for u in np.unique(t[e == 1]):
        s *= 1 - ((t == u) & (e == 1)).sum() / (t >= u).sum()
        xs.append(u)
        ys.append(s)
    return np.array(xs), np.array(ys)


def logrank_p(time, event, group):
    """two-group log-rank P, as cns.logrank_p"""
    from scipy.stats import chi2
    t, e, g = np.asarray(time, float), np.asarray(event, int), np.asarray(group, bool)
    oe, v = 0.0, 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u
        nn, n1 = at.sum(), (at & g).sum()
        d, d1 = ((t == u) & (e == 1)).sum(), ((t == u) & (e == 1) & g).sum()
        if nn > 1:
            oe += d1 - d * n1 / nn
            v += d * (n1 / nn) * (1 - n1 / nn) * (nn - d) / (nn - 1)
    return float(chi2.sf(oe ** 2 / v, 1))


def textw(s, size=5.0, weight="normal"):
    """width of a string in mm"""
    return TextPath((0, 0), s, prop=fm.FontProperties(family="Helvetica", size=size, weight=weight)).get_extents().width * 25.4 / 72


def fam_mark(pg, x, y, f, w=1.0, h=1.2):
    pg.box(x, y - h / 2, w, h, FAM[f], z=3)


pg = Page(240)

# ---------------------------------------------------------------- a  calls per endpoint
pg.panel(1, 3, "a", "Calls per endpoint")
ya, ha = 12.0, 19.0
ax = pg.ax(9.5, ya, 15.5, ha)
out = []
for i, e in enumerate(EP):
    r = row(F[e])
    ci(ax, 0, r.exp_true, i, SIG, lw=3.2)
    ci(ax, r.exp_true, r.called, i, COEX, lw=3.2)
    yy = ya + ha * (i + 0.5) / 3
    for xx, s in zip((29.2, 33.8, 38.6, 43.4), (f"{r.called:.0f}", f"{r.exp_false:.1f}", f"{r.fdr:.0f}%", f"{r.fdr_pi0:.0f}%")):
        pg.txt(xx, yy, s, ha="right")
    out.append((e, r.called, r.exp_true, r.exp_false, r.fdr, r.fdr_pi0))
for xx, s in zip((29.2, 33.8, 38.6, 43.4), ("called", "exp.", "FDR", "corr.")):
    pg.txt(xx, ya - 1.4, s, ha="right", weight="bold")
ax.set_yticks(range(3), EP, fontsize=5.5)
ax.set_ylim(2.5, -0.5)
ax.set_xlim(0, 26)
ax.set_xticks([0, 10, 20])
ax.set_xlabel("signatures called (of 52)")
key(pg, 4.6, 7.3, [(SIG, "expected true"), (COEX, "expected false")], gap=1.6)
pg.src("a", pd.DataFrame(out, columns=["endpoint", "called", "expected_true", "expected_false", "fdr_pct", "fdr_pi0_pct"]))

# ---------------------------------------------------------------- b  the estimate with intervals over correlation clusters
pg.panel(47.5, 3, "b", "Estimate and interval")
yb, hb = 9.6, 25.4
ax = pg.ax(66, yb, 24.5, hb)
DL = {"coex": "co-expr.", "coex closest": "closest", "rand": "random"}
rows_ = [(e, dc) for e in EP for dc in DL]
XC = 150
out, ys = [], []
for k, (e, dc) in enumerate(rows_):
    y = k + 0.55 * EP.index(e)
    ys.append(y)
    r = row(F[e], dc)
    ci(ax, r.fdr_lo, min(r.fdr_hi, XC), y, EST, lw=0.7)
    if r.fdr_hi > XC:
        ax.plot(XC, y, marker=">", ms=2.2, mfc=EST, mec="none", clip_on=False, zorder=4)
    point(ax, r.fdr_pi0, y, PI0, ms=2.9, filled=False, mew=0.5)
    point(ax, r.fdr, y, EST, ms=3.3)
    out.append((e, dc, r.fdr, r.fdr_lo, r.fdr_hi, r.fdr_pi0, r.fdr_pi0_lo, r.fdr_pi0_hi))
ymax = max(ys)
for y, o in zip(ys, out):
    ym = yb + hb * (y + 0.6) / (ymax + 1.2)
    pg.txt(102.4, ym, f"{o[2]:.0f} ({o[3]:.0f}–{o[4]:.0f})", ha="right")
    pg.txt(114.6, ym, f"{o[5]:.0f} ({o[6]:.0f}–{o[7]:.0f})", ha="right")
pg.txt(102.4, yb - 1.4, "estimate", ha="right", weight="bold")
pg.txt(114.6, yb - 1.4, "π0-corrected", ha="right", weight="bold")
key(pg, 66.2, yb - 1.4, [(EST, "estimate", "o"), (PI0, "π0-corrected", "open-o")], sw=1.4, gap=1.4)
ax.set_yticks(ys, [f"{e} {DL[dc]}" for e, dc in rows_], fontsize=5)
ax.set_ylim(ymax + 0.6, -0.6)
ax.set_xlim(0, XC)
ax.set_xticks([0, 50, 100, 150])
ax.set_xlabel("estimated FDR (%)")
pg.src("b", pd.DataFrame(out, columns=["endpoint", "decoys", "fdr", "fdr_lo", "fdr_hi", "fdr_pi0", "fdr_pi0_lo", "fdr_pi0_hi"]))

# ---------------------------------------------------------------- c  decoy co-expression
pg.panel(116.5, 3, "c", "Decoy co-expression")
ax = pg.ax(124.5, 9, 22, 26)
SET1 = INFO[INFO.set == 1]
for cls in ("rand", "coex"):
    g = SET1[SET1["class"] == cls]
    regions(ax, g.r_target_med, g.r_achieved_med, UNIT[cls], (0, 1), (-0.05, 1))
    ax.plot(g.r_target_med, g.r_achieved_med, ls="none", marker=MARK[cls], ms=0.7, mfc=UNIT[cls], mec="none", alpha=0.6, zorder=3, rasterized=True)
ax.plot([0, 1], [0, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=4)
ax.set_xlim(0, 1)
ax.set_ylim(-0.05, 1)
ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_xlabel("signature r")
ax.set_ylabel("decoy r")
key(pg, 125.2, 10.3, [(COEX, "co-expressed", "s")], sw=1.4)
key(pg, 125.2, 12.7, [(RAND, "random-gene", "D")], sw=1.4)
pg.src("c", SET1[["decoy", "target", "class", "r_target_med", "r_achieved_med"]])

# ---------------------------------------------------------------- d  ten decoy sets
pg.panel(149.5, 3, "d", "Ten decoy sets")
ax = pg.ax(158, 9, 21, 26)
out = []
for i, e in enumerate(EP):
    b = B[e][B[e].stat == "fold"]
    r = row(F[e])
    point(ax, i + np.linspace(-0.26, 0.26, len(b)), b.fdr.to_numpy(), EST, ms=2.3)
    ax.plot([i - 0.32, i + 0.32], [r.fdr] * 2, color=INK, lw=0.8, zorder=5, solid_capstyle="butt")
    ax.text(i, b.fdr.max() + 2.2, f"{r.fdr:.1f}", ha="center", va="bottom", fontsize=5)
    assert abs(b.fdr.min() - r.fdr_set_min) < 1e-6 and abs(b.fdr.max() - r.fdr_set_max) < 1e-6
    out += [(e, int(s), f) for s, f in zip(b.set, b.fdr)] + [(e, "pooled", r.fdr)]
ax.set_xticks(range(3), EP, fontsize=5.5)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.55, 2.55)
ax.set_ylim(0, 80)
ax.set_yticks([0, 25, 50, 75])
ax.set_ylabel("estimated FDR (%)")
pg.src("d", pd.DataFrame(out, columns=["endpoint", "set", "fdr_pct"]))

# ---------------------------------------------------------------- e, f, g  cancer types, shared rows
TORD = UT.groupby("cohort").fold_med.median().sort_values(ascending=False).index.tolist()
NT = UT.drop_duplicates("cohort").set_index("cohort").n
Y0, PT = 49.4, 2.72
HT = PT * len(TORD)


def ymm(i):
    return Y0 + PT * (i + 0.5)


pg.panel(1, 43, "e", "Null median folded C")
ax = pg.ax(9.4, Y0, 19.6, HT)
for i, co in enumerate(TORD):
    v = UT[UT.cohort == co].fold_med.to_numpy()
    q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
    point(ax, v, i + swarm(v, 0.0015, 0.065, 0.42), NULLC, ms=1.5, mew=0)
    ci(ax, q[0], q[4], i, MUTED, lw=0.5, zorder=5)
    ci(ax, q[1], q[3], i, INK, lw=1.2, zorder=5)
    ax.plot([q[2], q[2]], [i - 0.42, i + 0.42], color=INK, lw=0.7, zorder=6, solid_capstyle="butt")
    pg.txt(8.8, ymm(i), co, ha="right")
    pg.txt(35.2, ymm(i), f"{NT[co]:,}", ha="right")
pg.txt(35.2, Y0 - 1.5, "n", ha="right", weight="bold")
zero(ax, 0.5)
ax.set_ylim(len(TORD) - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlim(0.495, 0.655)
ax.set_xticks([0.5, 0.55, 0.6, 0.65], ["0.50", "0.55", "0.60", "0.65"])
ax.set_xlabel("null median folded C")
pg.src("e", UT[["cohort", "n", "target", "fold_med"]])

pg.panel(37.5, 43, "f", "Above the null")
ax = pg.ax(38.6, Y0, 18, HT)
above = UO.assign(above=UO.fold_obs > UO.fold_p95)
A = above.assign(unit_class=above["class"].where(above.kind == "decoy", "signature")).groupby(["cohort", "unit_class"]).above.agg(["sum", "size"])
out = []
for i, co in enumerate(TORD):
    s, c, rr = (A.loc[(co, k)] for k in ("signature", "coex", "rand"))
    ax.plot([0, 60], [i, i], color=EMPTY, lw=0.3, zorder=1)
    point(ax, 100 * rr["sum"] / rr["size"], i, RAND, "D", ms=2.0, zorder=3)
    point(ax, 100 * c["sum"] / c["size"], i, COEX, "s", ms=2.2, zorder=4)
    point(ax, 100 * s["sum"] / s["size"], i, SIG, "o", ms=2.6, zorder=5)
    pg.txt(60.8, ymm(i), f"{int(s['sum'])}", ha="right")
    out.append((co, int(s["sum"]), int(s["size"]), 100 * s["sum"] / s["size"], 100 * c["sum"] / c["size"], int(c["size"]), 100 * rr["sum"] / rr["size"], int(rr["size"])))
pg.txt(60.8, Y0 - 1.5, "of 52", ha="right", weight="bold")
zero(ax, 5)
ax.set_ylim(len(TORD) - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlim(-2, 60)
ax.set_xticks([0, 25, 50])
ax.set_xlabel("above null 95th\npercentile (%)", linespacing=1.0)
for k, (c, l_, m_) in enumerate(((SIG, "signatures", "o"), (COEX, "co-expressed", "s"), (RAND, "random-gene", "D"))):
    key(pg, 45.0, ymm(15 + k), [(c, l_, m_)], sw=1.4)
pg.src("f", pd.DataFrame(out, columns=["cohort", "signatures_above", "signatures", "signatures_pct", "coex_decoys_pct", "coex_decoy_tests", "rand_decoys_pct", "rand_decoy_tests"]))

pg.panel(62.5, 43, "g", "Shortfall by type")
ax = pg.ax(63, Y0, 21.5, HT)
RP["gap"] = RP.r_target - RP.r_decoy_med
out = []
for i, co in enumerate(TORD):
    q = np.nanpercentile(RP[RP.cohort == co].gap, [5, 25, 50, 75, 95])
    ci(ax, q[0], q[4], i, tint(COEX, 0.5), lw=0.6)
    ci(ax, q[1], q[3], i, COEX, lw=1.4)
    point(ax, q[2], i, "white", "o", ms=1.9, zorder=5, mew=0.4)
    ax.plot(q[2], i, marker="o", ms=1.9, mfc="white", mec=COEX, mew=0.5, zorder=6)
    out.append((co, *q))
zero(ax, 0)
ax.set_ylim(len(TORD) - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([0, 0.2, 0.4], ["0", "0.2", "0.4"])
ax.set_xlabel("signature − decoy r")
pg.src("g", pd.DataFrame(out, columns=["cohort", "q05", "q25", "median", "q75", "q95"]))

# ---------------------------------------------------------------- h  permuted survival
pg.panel(87.5, 43, "h", "Permuted survival")
ax = pg.ax(95.5, 48.5, 18, 21.5)
g = GP[GP.stat == "fold"]
gs = GN[GN.stat == "fold"].iloc[0]
hist(ax, g.R, np.arange(-0.5, 18.5))
ax.axvline(g.R.mean(), color=INK, lw=0.7, zorder=4)
ax.axvline(g.E.mean(), color=EST, lw=0.7, ls=(0, (2.2, 1.6)), zorder=5)
ax.set_xlim(-0.8, 17.8)
ax.set_xticks([0, 5, 10, 15])
ax.set_ylim(0, 56)
ax.set_yticks([0, 25, 50])
ax.set_xlabel("false calls per library")
ax.set_ylabel("permutations")
ax.text(17.8, 56, f"ratio {gs.ratio:.2f}\n({gs.ratio_lo:.2f}–{gs.ratio_hi:.2f})", ha="right", va="top", fontsize=5, linespacing=1.1)
key(pg, 100.6, 56.4, [(INK, f"mean {half_up(g.R.mean())}", "line")], sw=1.4)
pg.line(100.6, 58.8, 102.0, 58.8, color=EST, lw=0.8, ls=(0, (2.2, 1.0)))
pg.txt(102.6, 58.8, f"estimate {half_up(g.E.mean())}")
pg.src("h", g.assign(mean_realized=g.R.mean(), mean_estimate=g.E.mean(), ratio=gs.ratio, ratio_lo=gs.ratio_lo, ratio_hi=gs.ratio_hi))

# ---------------------------------------------------------------- i  permuted survival, per signature
pg.panel(115.5, 43, "i", "Call rate, permuted")
ax = pg.ax(123.5, 48.5, 18.5, 21.5)
gt = GT[GT.stat == "fold"].set_index("target")
OSC = CL["OS"][CL["OS"].kind == "target"].set_index("target").call
ax.plot([0, 9], [0, 9], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
for c_ in (False, True):
    m_ = OSC.reindex(gt.index).astype(bool) == c_
    point(ax, 100 * gt.decoy_rate[m_], 100 * gt.rate[m_], SIG, ms=2.6, filled=c_, mew=None if c_ else 0.5, zorder=5 if c_ else 4)
ax.set_xlim(0, 9)
ax.set_ylim(0, 9)
ax.set_xticks([0, 3, 6, 9])
ax.set_yticks([0, 3, 6, 9])
ax.set_xlabel("decoys (%)")
ax.set_ylabel("signature (%)")
ax.text(8.9, 0.25, f"r = {np.corrcoef(gt.decoy_rate, gt.rate)[0, 1]:.2f}, n = {len(gt)}", ha="right", va="bottom", fontsize=5)
gt = gt.assign(called=OSC.reindex(gt.index)).reset_index()
pg.src("i", gt[["target", "rate", "decoy_rate", "called"]])

# ---------------------------------------------------------------- j  response against survival
pg.panel(144.5, 43, "j", "Response vs survival")
ax = pg.ax(153, 48.5, 26, 21.5)
ost = CL["OS"][CL["OS"].kind == "target"].set_index("target")
rsp = U0.set_index("unit")
J = pd.DataFrame({"resp": rsp.med_delta, "os": ost.med_delta, "rc": rsp.call.astype(str) == "True", "oc": ost.call}).dropna()
J["family"] = J.index.map(FAMILY)
for s, r in J.iterrows():
    c = FAM[r.family]
    if r.rc and r.oc:
        point(ax, r.resp, r.os, c, ms=3.0, zorder=6)
    elif r.rc or r.oc:
        point(ax, r.resp, r.os, c, ms=2.8, filled=False, zorder=5, mew=0.6)
    else:
        ax.plot(r.resp, r.os, marker="o", ms=1.8, mfc=tint(c, 0.45), mec="none", ls="none", zorder=3)
zero(ax, 0)
zero(ax, 0, axis="y")
ax.set_xlabel("median ΔAUC, response")
ax.set_ylabel("median ΔC, survival")
ax.set_xlim(-0.012, 0.112)
ax.set_ylim(-0.0135, 0.061)
ax.set_xticks([0, 0.05, 0.1], ["0", "0.05", "0.1"])
ax.set_yticks([-0.01, 0, 0.01, 0.02, 0.03], ["−0.01", "0", "0.01", "0.02", "0.03"])
both, one, none = int((J.rc & J.oc).sum()), int((J.rc ^ J.oc).sum()), int((~J.rc & ~J.oc).sum())
pg.src("j", J.reset_index(names="signature").rename(columns={"rc": "called_response", "oc": "called_os"}))
key(pg, 154.2, 49.6, [(INK, f"both ({both})", "o")], sw=1.4)
key(pg, 154.2, 51.8, [(INK, f"one ({one})", "open-o")], sw=1.4)
pg.c.plot(154.9, 54.0, "o", ms=1.8, mfc=tint(INK, 0.45), mec="none", zorder=3)
pg.txt(156.2, 54.0, f"neither ({none})")
for k, f_ in enumerate(("Inflamed/IFN", "Confounder", "Other")):
    fam_mark(pg, 167.2, 49.6 + 2.2 * k, f_, w=1.2, h=1.4)
    pg.txt(168.9, 49.6 + 2.2 * k, FLAB[f_])

# ---------------------------------------------------------------- k  the example across cancer types
pg.panel(87.5, 80.5, "k", "CIN70 (Carter) and its 32 co-expressed decoys, every cancer type")
ax = pg.ax(95.5, 85.5, 83.5, 19.5)
Ud = UO[(UO.target == EXS) & (((UO.kind == "decoy") & (UO["class"] == "coex") & (UO.set == 1)) | (UO.kind == "target"))]
tg = Ud[Ud.kind == "target"].drop_duplicates("cohort").set_index("cohort").dC.sort_values(ascending=False)
for i, co in enumerate(tg.index):
    dv = Ud[(Ud.cohort == co) & (Ud.kind == "decoy")].dC.to_numpy()
    point(ax, i + swarm(dv, 0.0065, 0.09, 0.38), dv, COEX, "s", ms=1.5, zorder=3, mew=0.15)
    ax.plot([i - 0.32, i + 0.32], [np.median(dv)] * 2, color=INK, lw=0.7, zorder=4, solid_capstyle="butt")
    point(ax, i, tg[co], SIG, ms=3.4, zorder=5)
zero(ax, 0, axis="y")
ax.set_xticks(range(len(tg)), tg.index, rotation=90, fontsize=5)
ax.tick_params(axis="x", length=0, pad=1)
ax.set_xlim(-0.6, len(tg) - 0.4)
ax.set_ylim(-0.1, 0.22)
ax.set_yticks([-0.1, 0, 0.1, 0.2], ["−0.1", "0", "0.1", "0.2"])
ax.set_ylabel("ΔC")
key(pg, 136.5, 87, [(SIG, "CIN70 (Carter)", "o"), (COEX, "32 decoys", "s"), (INK, "decoy median", "line")], sw=1.4, gap=1.5)
pg.src("k", Ud[["cohort", "kind", "unit", "dC"]].assign(unit=lambda d: d.unit.where(d.kind == "decoy", EXS)))

# ---------------------------------------------------------------- l, m, n  every signature, shared rows
T0 = CL["OS"][CL["OS"].kind == "target"].set_index("target")
ORD = T0.med_delta.sort_values(ascending=False).index.tolist()
Y1, PS = 127.2, 1.95
HS = PS * len(ORD)


def ysig(i):
    return Y1 + PS * (i + 0.5)


pg.panel(1, 118, "l", "Signature and 320 co-expressed decoys")
ax = pg.ax(25.5, Y1, 31.5, HS)
dq = CL["OS"][(CL["OS"].kind == "decoy") & (CL["OS"]["class"] == "coex")].groupby("target").med_delta.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
out = []
for i, s in enumerate(ORD):
    q = dq.loc[s].to_numpy()
    ci(ax, q[0], q[4], i, tint(COEX, 0.45), lw=0.6, zorder=2)
    ci(ax, q[1], q[3], i, COEX, lw=1.4, zorder=3)
    ax.plot([q[2], q[2]], [i - 0.36, i + 0.36], color=INK, lw=0.5, zorder=4, solid_capstyle="butt")
    point(ax, T0.loc[s, "med_delta"], i, SIG, ms=2.8, filled=bool(T0.loc[s, "call"]), zorder=5)
    pg.txt(23.4, ysig(i), lab(s), ha="right")
    fam_mark(pg, 23.8, ysig(i), FAMILY[s], w=1.0, h=PS - 0.15)
    out.append((s, FAMILY[s], T0.loc[s, "med_delta"], bool(T0.loc[s, "call"]), *q))
zero(ax, 0)
ax.set_ylim(len(ORD) - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlim(-0.022, 0.05)
ax.set_xticks([-0.02, 0, 0.02, 0.04], ["−0.02", "0", "0.02", "0.04"])
ax.set_xlabel("median ΔC, overall survival")
pg.src("l", pd.DataFrame(out, columns=["signature", "family", "median_dC", "called", "d05", "d25", "d50", "d75", "d95"]))
key(pg, 46.4, ysig(30.5), [(SIG, "called", "o")], sw=1.4)
key(pg, 46.4, ysig(31.8), [(SIG, "not called", "open-o")], sw=1.4)
pg.txt(46.4, ysig(33.1), "decoys")
pg.line(47.1, ysig(34.4) - 0.7, 47.1, ysig(34.4) + 0.7, color=INK, lw=0.5)
pg.txt(48.4, ysig(34.4), "median")
pg.line(46.4, ysig(35.7), 47.8, ysig(35.7), color=COEX, lw=1.4)
pg.txt(48.4, ysig(35.7), "25–75%")
pg.line(46.4, ysig(37), 47.8, ysig(37), color=tint(COEX, 0.45), lw=0.6)
pg.txt(48.4, ysig(37), "5–95%")
for k, f_ in enumerate(("Inflamed/IFN", "Confounder", "Other")):
    fam_mark(pg, 9.0, 231.2 + 2.4 * k, f_, w=1.2, h=1.5)
    pg.txt(10.8, 231.2 + 2.4 * k, FLAB[f_])

# m  median gain and calls per endpoint
pg.panel(58.5, 118, "m", "Endpoints")
xm, cwm = 59.2, 2.8
LM = 0.03
out = []
for j, e in enumerate(EP):
    t = CL[e][CL[e].kind == "target"].set_index("target")
    pg.txt(xm + cwm * j + cwm / 2, Y1 - 0.8, e, ha="center", va="bottom", rotation=90)
    for i, s in enumerate(ORD):
        v, c = t.loc[s, "med_delta"], bool(t.loc[s, "call"])
        tile(pg, xm + cwm * j, Y1 + PS * i, cwm, PS, TEALOR((v + LM) / (2 * LM)))
        if c:
            pg.c.plot(xm + cwm * j + (cwm - 0.15) / 2, Y1 + PS * i + (PS - 0.15) / 2, "o", ms=1.7, mfc=SIG, mec="white", mew=0.3, zorder=4)
        out.append((s, e, v, c))
pg.src("m", pd.DataFrame(out, columns=["signature", "endpoint", "median_dC", "called"]))
pg.cbar(59.2, 231.0, 11, 1.4, TEALOR, -LM, LM, [-0.03, 0, 0.03], "median ΔC", labels=["−0.03", "0", "≥ 0.03"])

# n  gain per signature and cancer type, same rows
pg.panel(75, 118, "n", "Gain over the null, by type")
HM = UT.pivot_table(index="cohort", columns="target", values="dC")
tord = HM.median(axis=1).sort_values(ascending=False).index
xn, cwn = 70.2, 2.4
HL = 0.1
for j, co in enumerate(tord):
    pg.txt(xn + cwn * j + cwn / 2, Y1 - 0.8, co, ha="center", va="bottom", rotation=90)
    for i, s in enumerate(ORD):
        tile(pg, xn + cwn * j, Y1 + PS * i, cwn, PS, TEALOR((np.clip(HM.loc[co, s], -HL, HL) + HL) / (2 * HL)))
pg.cbar(90, 231.0, 22, 1.4, TEALOR, -HL, HL, [-0.1, 0, 0.1], "ΔC", labels=["−0.1", "0", "≥ 0.1"])
pg.src("n", UT[["cohort", "target", "fold_obs", "fold_med", "dC"]])

# ---------------------------------------------------------------- o, p  Kaplan-Meier pair
XMAX, TICKS = 10, [0, 2.5, 5, 7.5, 10]
KMY, KMH = 123, 20
KM_OUT = {}


def km_pair(x, w, col, c, letter, title, show_y):
    pg.panel(x - (7.4 if show_y else 3.0), 118, letter, title)
    a = pg.ax(x, KMY, w, KMH)
    v = KMD[col].to_numpy()
    hi = v > np.median(v)
    t, e = KMD.years.to_numpy(), KMD.event.to_numpy()
    lo_c = tint(c, 0.45)
    for m, cc in ((hi, c), (~hi, lo_c)):
        xs, ys_ = km_curve(t[m], e[m])
        a.step(xs, ys_, where="post", color=cc, lw=0.8, zorder=3)
        ct = t[m][e[m] == 0]
        a.plot(ct, ys_[np.clip(np.searchsorted(xs, ct, side="right") - 1, 0, len(ys_) - 1)], ls="none", marker="|", ms=1.4, mew=0.3, color=cc, zorder=4)
    p = logrank_p(t, e, hi)
    a.set_xlim(0, XMAX)
    a.set_ylim(0, 1.04)
    a.set_xticks(TICKS, ["0", "2.5", "5", "7.5", "10"])
    a.set_yticks([0, 0.5, 1], ["0", "0.5", "1"] if show_y else [])
    a.set_xlabel("years")
    if show_y:
        a.set_ylabel("overall survival")
    xl = x + 0.8
    yl = KMY + KMH - 8.0
    pg.line(xl, yl, xl + 1.8, yl, color=c, lw=0.8)
    pg.txt(xl + 2.4, yl, "above median")
    pg.line(xl, yl + 2.2, xl + 1.8, yl + 2.2, color=lo_c, lw=0.8)
    pg.txt(xl + 2.4, yl + 2.2, "below")
    # two figures down to 0.001; below, no superscript (p_small: 0.0004, or < 0.0001)
    if p >= 0.001:
        pg.txt(xl, yl + 4.5, f"log-rank P = {p:.2f}")
    else:
        op, s_ = p_small(p)
        pg.txt(xl, yl + 4.5, f"log-rank P {op} {s_}")
    rows_y = (KMY + KMH + 9.3, KMY + KMH + 11.6)
    for yy, m, cc, nm in ((rows_y[0], hi, c, "above"), (rows_y[1], ~hi, INK, "below")):
        for tt in TICKS:
            pg.txt(x + w * tt / XMAX, yy, str(int((t[m] >= tt).sum())), ha="center", color=cc)
        if show_y:
            pg.txt(x - 2.6, yy, nm, ha="right")
    if show_y:
        pg.txt(x - 2.6, KMY + KMH + 7.0, "at risk", ha="right", weight="bold")
    KM_OUT[letter] = p
    return p


km_pair(128.4, 23, "target", SIG, "o", "CIN70 (Carter), KIRP", True)
km_pair(155, 23, "decoy", COEX, "p", "Median decoy, KIRP", False)
km_src = KMD[["time", "event", "target", "decoy"]].sort_values("time").reset_index(drop=True)
pg.src("o", km_src[["time", "event", "target"]].assign(log_rank_P=KM_OUT["o"]).rename(columns={"time": "days", "target": "score"}))
pg.src("p", km_src[["time", "event", "decoy"]].assign(log_rank_P=KM_OUT["p"]).rename(columns={"time": "days", "decoy": "score"}))

# ---------------------------------------------------------------- q  decoys called per signature
pg.panel(121, 160, "q", "Decoys called, per signature")
ax = pg.ax(129, 165.5, 50, 26.5)
out = []
for i, e in enumerate(EP):
    d = CL[e][CL[e].kind == "decoy"].groupby(["target", "class"]).call.mean().unstack()
    for off, cls in ((-0.21, "coex"), (0.21, "rand")):
        v = 100 * d[cls]
        strip2(ax, i + off, v, UNIT[cls], m=MARK[cls], ms=1.9 if cls == "coex" else 1.7, bin_w=1.25 if cls == "coex" else 1.25, step=0.042, limit=0.19, med_half=0.19, mew=0.15)
        ax.text(i + off, 38, f"{row(F[e], cls).decoy_rate:.1f}", ha="center", va="top", fontsize=5)
        out += [(e, cls, s, x_) for s, x_ in v.items()]
ax.set_xticks(range(3), EP, fontsize=5.5)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.55, 2.55)
ax.set_ylim(-1, 38)
ax.set_yticks([0, 10, 20, 30])
ax.set_ylabel("decoys called (%)")
key(pg, 152.6, 163.0, [(COEX, "co-expressed", "s"), (RAND, "random-gene", "D")], sw=1.4, gap=1.2)
pg.src("q", pd.DataFrame(out, columns=["endpoint", "class", "signature", "decoys_called_pct"]))

# ---------------------------------------------------------------- r  rank among decoys
pg.panel(121, 198.5, "r", "Rank among decoys")
ax = pg.ax(129, 204, 19.5, 25.5)


def rank_q(calls):
    """q as decoy_fdr() computes it for pi0: (1 + decoys at or above the target's signed -log10 P) / (1 + decoys)"""
    st = lambda x: np.sign(x.med_delta) * -np.log10(np.maximum(x.p, 1e-300))  # noqa: E731
    t = calls[calls.kind == "target"]
    d = calls[(calls.kind == "decoy") & (calls["class"] == "coex")]
    sd = st(d).groupby(d.target)
    return pd.Series({s: (1 + (sd.get_group(s) >= v - 1e-12).sum()) / (1 + sd.get_group(s).size) for s, v in zip(t.target, st(t))})


QS = {e: rank_q(CL[e]) for e in EP}
for e in EP:
    assert abs(2 * (QS[e] > 0.5).mean() - row(F[e]).pi0) < 1e-9, e
hist(ax, QS["OS"], np.linspace(0, 1, 11))
ax.axhline(row(F["OS"]).pi0 * 52 / 10, color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=4)
ax.set_xlim(0, 1)
ax.set_ylim(0, 24)
ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_yticks([0, 10, 20])
ax.set_xlabel("q, overall survival")
ax.text(0.98, 23.6 - 2.9 * 3, "n = 52", ha="right", va="top", fontsize=5)
ax.set_ylabel("signatures")
for k, e in enumerate(EP):
    ax.text(0.98, 23.6 - 2.9 * k, f"π0 {e} {row(F[e]).pi0:.2f}", ha="right", va="top", fontsize=5)
pg.src("r", pd.DataFrame(QS).reset_index(names="signature").rename(columns={e: f"q_{e}" for e in EP}))

# ---------------------------------------------------------------- s  called, by family
pg.panel(150.5, 198.5, "s", "Called, by family")
x0, y0, cw, ch = 164.2, 207.0, 3.7, 5.6
HD = ("resp.", "OS", "DSS", "PFI")
FORD = ("Inflamed/IFN", "Confounder", "Other")
TAB = {"resp.": J.groupby("family").rc.sum().to_dict()}
for e in EP:
    t = CL[e][CL[e].kind == "target"]
    TAB[e] = t[t.call].target.map(FAMILY).value_counts().to_dict()
NF = pd.Series(FAMILY).value_counts()
out = []
for j, h in enumerate(HD):
    pg.txt(x0 + cw * j + cw / 2, y0 - 1.3, h, ha="center", weight="bold")
for i, f_ in enumerate(FORD + ("all",)):
    nf = NF[f_] if f_ != "all" else int(NF.sum())
    if f_ != "all":
        fam_mark(pg, 150.6, y0 + ch * i + ch / 2 - 0.08, f_, w=1.1, h=ch - 0.15)
    pg.txt(152.2, y0 + ch * i + ch / 2 - 1.0, FLAB.get(f_, f_))
    pg.txt(152.2, y0 + ch * i + ch / 2 + 0.9, f"n = {nf}")
    for j, h in enumerate(HD):
        n_ = int(TAB[h].get(f_, 0)) if f_ != "all" else int(sum(TAB[h].values()))
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(n_ / nf), str(n_))
        out.append((f_, int(nf), h, n_))
pg.src("s", pd.DataFrame(out, columns=["family", "signatures", "endpoint", "called"]))
pg.cbar(164.2, 231.0, 12.6, 1.4, BLUE, 0, 100, [0, 50, 100], None, labels=["0", "50", "100"])
pg.txt(151.0, 231.7, "family called (%)")

pg.save("Figure5")
print("called", {e: (int(row(F[e]).called), round(row(F[e]).fdr, 1), round(row(F[e]).fdr_pi0, 1)) for e in EP}, "both/one/neither", both, one, none, "km P", KM_OUT)
