"""Figure 4 redrawn at the density and in the look of the AML multi-omic figures (pixl_nc.py).
A target-decoy false discovery rate for 52 immune signatures and response to checkpoint blockade.
Tables as in fig4_v10.py: calib_v3_fdr.csv, calib_v3_units.csv, calib_v3/ICB_*_v3.csv, calib_v3_fdr_diff_close.csv, calib_v3_rcurve.csv
(first decoy set); v12_pkg/response_fdr.csv, response_fdr_by_set.csv, response_calls.csv.gz, response_decoy_info.csv.gz (ten sets);
v13_fig/threshold_curve.csv, convergence_sets.csv (written by fig4_v10.py). New panels are descriptive reads of those tables
(validation/NOTE_figures_v17_descriptive_only.md). The per-test intervals of panel h repeat the seeded bootstrap of fig4_v10.py,
same seed, same order, same draws; nothing is resampled afresh."""
import glob
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt
from matplotlib import font_manager as fm
from matplotlib.textpath import TextPath

I_ = D + "/icb/"
P_ = I_ + "v12_pkg/"
F = pd.read_csv(I_ + "calib_v3_fdr.csv")
U = pd.read_csv(I_ + "calib_v3_units.csv")
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(I_ + "calib_v3/ICB_*_v3.csv"))])
SC = pd.read_csv(I_ + "calib_v3_fdr_diff_close.csv")
RCV = pd.read_csv(I_ + "calib_v3_rcurve.csv")
RF = pd.read_csv(P_ + "response_fdr.csv")
RB = pd.read_csv(P_ + "response_fdr_by_set.csv")
RC = pd.read_csv(P_ + "response_calls.csv.gz")
RI = pd.read_csv(P_ + "response_decoy_info.csv.gz", usecols=["decoy", "target", "class", "set", "r_target_med", "r_achieved_med"])
TH = pd.read_csv(I_ + "v13_fig/threshold_curve.csv")
CONV = pd.read_csv(I_ + "v13_fig/convergence_sets.csv")
FAMS = pd.read_csv(I_ + "v5_signature_clusters.csv").set_index("signature").family

# The two statistics. Black and blue as in the old sheet's dark and blue; checked with their tints.
FOLD, DIRN = INK, "#3B8ED0"
STAT = {"fold": FOLD, "dir": DIRN}
LIGHT = {"fold": 0.40, "dir": 0.30}
check_palette("statistic", {"fold": FOLD, "dir": DIRN, "fold false": tint(FOLD, LIGHT["fold"]), "dir false": tint(DIRN, LIGHT["dir"])})

SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)",
         "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)", "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}
GENE_NAMED = ()   # signature names are labels, roman on every sheet
FP = fm.FontProperties(family="Helvetica", size=5)


def lab(s):
    if s in SHORT:
        return SHORT[s]
    a, b = s.rsplit("_", 1)
    return f"{a.replace('_', ' ')} ({b})"


def width_mm(s, size=5.0):
    fp = fm.FontProperties(family="Helvetica", size=size)
    return TextPath((0, 0), s, prop=fp).get_extents().width * 25.4 / 72


def name_right(pg, x, y, s):
    """a signature name right-aligned at x, the gene-named part in italics"""
    t = lab(s)
    head = t.split(" (")[0]
    if head in GENE_NAMED:
        tail = "(" + t.split(" (")[1]
        pg.txt(x, y, tail, ha="right")
        pg.txt(x - width_mm(tail) - 0.55, y, head, ha="right", style="italic")
    else:
        pg.txt(x, y, t, ha="right")


def fmt(v):
    return f"{v:.0f}" if round(v, 1) >= 10 else f"{v:.1f}"


def rng_txt(a, b):
    return f"{fmt(a)}–{fmt(b)}"


pg = Page(240)

# ================================================================== a  the estimate, by decoy class (old f)
pg.panel(1, 3, "a", "Estimated FDR, by decoy class")
ax = pg.ax(21.5, 9.2, 44.5, 38.3)
CLS = [("coex", "co-expressed", "s"), ("coex closest", "closest-matched", "s"), ("rand", "random-gene", "D")]
ST = [("fold", "folded"), ("dir", "directed")]
units_y = []
yy = 0.0
for c_i in range(3):
    for s_i in range(2):
        units_y.append((c_i, s_i, yy))
        yy += 2 + (0.35 if s_i == 0 else 0)
    yy += 0.75 - 0.0
span = yy - 0.75
ax.set_ylim(span - 0.45, -0.55)
y2mm = lambda v: 9.2 + 38.3 * (v + 0.55) / (span - 0.45 + 0.55)  # noqa: E731
cols = {"called": 69.0, "false": 76.6, "fdr": 84.0, "sig": 92.0, "clu": 101.2, "pi0": 108.0}
for k, h in (("called", "called"), ("false", "exp. false"), ("fdr", "FDR"), ("sig", "signatures"), ("clu", "clusters"), ("pi0", "π0")):
    pg.txt(cols[k], 7.6, h, ha="center", weight="bold")
pg.txt((cols["sig"] + cols["clu"]) / 2, 5.3, "95% interval, over", ha="center", weight="bold")
out = []
for c_i, s_i, y in units_y:
    dc, cname, mk = CLS[c_i]
    st, sname = ST[s_i]
    rc = RF[(RF.stat == st) & (RF.decoys == dc) & (RF.boot == "clusters")].iloc[0]
    rt = RF[(RF.stat == st) & (RF.decoys == dc) & (RF.boot == "targets")].iloc[0]
    # plain estimate, then the pi0-corrected estimate on the line below
    ci(ax, rc.fdr_lo, rc.fdr_hi, y, tint(EST, 0.45), lw=0.5, zorder=2)
    ci(ax, rt.fdr_lo, rt.fdr_hi, y, EST, lw=1.2, zorder=3)
    point(ax, rc.fdr, y, EST, mk, ms=3.6, zorder=5)
    ci(ax, rc.fdr_pi0_lo, rc.fdr_pi0_hi, y + 1, tint(PI0, 0.45), lw=0.5, zorder=2)
    ci(ax, rt.fdr_pi0_lo, rt.fdr_pi0_hi, y + 1, PI0, lw=1.2, zorder=3)
    point(ax, rc.fdr_pi0, y + 1, PI0, mk, ms=3.2, filled=False, zorder=5)
    yt, yb = y2mm(y), y2mm(y + 1)
    pg.txt(cols["called"], yt, f"{int(rc.called)}", ha="center")
    pg.txt(cols["false"], yt, f"{rc.exp_false:.1f}", ha="center")
    pg.txt(cols["fdr"], yt, f"{fmt(rc.fdr)}%", ha="center")
    pg.txt(cols["sig"], yt, rng_txt(rt.fdr_lo, rt.fdr_hi), ha="center")
    pg.txt(cols["clu"], yt, rng_txt(rc.fdr_lo, rc.fdr_hi), ha="center")
    pg.txt(cols["fdr"], yb, f"{fmt(rc.fdr_pi0)}%", ha="center")
    pg.txt(cols["sig"], yb, rng_txt(rt.fdr_pi0_lo, rt.fdr_pi0_hi), ha="center")
    pg.txt(cols["clu"], yb, rng_txt(rc.fdr_pi0_lo, rc.fdr_pi0_hi), ha="center")
    pg.txt(cols["pi0"], yb, f"{rc.pi0:.2f}", ha="center")
    pg.txt(20.9, (yt + yb) / 2, sname, ha="right")
    if s_i == 0:
        pg.txt(1.2, (yb + y2mm(y + 2.35)) / 2, cname, ha="left")
    out.append(dict(decoys=cname, statistic=sname, m=int(rc.m), called=int(rc.called), exp_false=rc.exp_false, fdr=rc.fdr, fdr_lo_signatures=rt.fdr_lo,
                    fdr_hi_signatures=rt.fdr_hi, fdr_lo_clusters=rc.fdr_lo, fdr_hi_clusters=rc.fdr_hi, pi0=rc.pi0, fdr_pi0=rc.fdr_pi0,
                    fdr_pi0_lo_signatures=rt.fdr_pi0_lo, fdr_pi0_hi_signatures=rt.fdr_pi0_hi, fdr_pi0_lo_clusters=rc.fdr_pi0_lo, fdr_pi0_hi_clusters=rc.fdr_pi0_hi))
zero(ax, 100)
ax.set_xscale("log"); ax.set_xlim(0.25, 320); ax.set_xticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"]); ax.minorticks_off()
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("estimated FDR (%)")
key(pg, 68, 49.8, [(EST, "estimate", "s"), (PI0, "π0-corrected", "open-s")], gap=2.0)
key(pg, 68, 52.6, [(INK, "signatures", "line")], gap=1.4)
pg.line(83.3, 52.6, 85.1, 52.6, color=tint(INK, 0.45), lw=0.5)
pg.txt(85.7, 52.6, "clusters")
pg.src("a", pd.DataFrame(out))

# ================================================================== b, c  calls with their composition (old b, d)
BX, BW = 133.5, 23.5
XMAX = 31
BCOL = {"called": 160.2, "true": 165.8, "false": 171.0, "fdr": 176.6}


def stack(ax, y, called, false, st):
    true_ = called - false
    ax.barh(y, true_, height=0.5, left=0, color=STAT[st], lw=0, zorder=3)
    ax.barh(y, false, height=0.5, left=true_, color=tint(STAT[st], LIGHT[st]), lw=0, zorder=3)


def row_values(y_mm, called, false, fdr):
    pg.txt(BCOL["called"], y_mm, f"{int(called)}", ha="center")
    pg.txt(BCOL["true"], y_mm, f"{called - false:.1f}", ha="center")
    pg.txt(BCOL["false"], y_mm, f"{false:.1f}", ha="center")
    pg.txt(BCOL["fdr"], y_mm, f"{fdr:.0f}%", ha="center")


pg.panel(112, 3, "b", "Signatures called of 52, by null and statistic")
for k, h in (("called", "called"), ("true", "true"), ("false", "false"), ("fdr", "FDR")):
    pg.txt(BCOL[k], 7.6, h, ha="center", weight="bold")
fa = F[(F.pool == "A") & (F.decoys == "coex")].set_index(["null", "stat"])
BY, BH = 9.2, 22.6
ax = pg.ax(BX, BY, BW, BH)
ax.set_ylim(11.65, -0.55)
y2b = lambda v: BY + BH * (v + 0.55) / 12.2  # noqa: E731
GAP = {"N0": 0, "N1": 0.6, "N2": 0.9, "N1nat": 1.5, "N2nat": 1.8}
NAMES = {"N0": "N0", "N1": "N1", "N2": "N2", "N1nat": "N1", "N2nat": "N2"}
order = ("N0", "N1", "N2", "N1nat", "N2nat")
out = []
for i, nl in enumerate(order):
    for j, (st, let) in enumerate((("fold", "F"), ("dir", "D"))):
        y = 2 * i + j + GAP[nl]
        r = fa.loc[(nl, st)]
        stack(ax, y, r.called, r.exp_false, st)
        pg.txt(BX - 0.6, y2b(y), let, ha="right")
        row_values(y2b(y), r.called, r.exp_false, r.fdr)
        out.append(dict(null=nl, statistic=st, called=int(r.called), exp_true=r.exp_true, exp_false=r.exp_false, fdr=r.fdr))
    pg.txt(BX - 3.1, y2b(2 * i + 0.5 + GAP[nl]), NAMES[nl], ha="right", weight="bold")
for grp, labt in ((("N0",), "size, engine"), (("N1", "N2"), "enriched pool"), (("N1nat", "N2nat"), "natural pool")):
    ys = [2 * order.index(n) + j + GAP[n] for n in grp for j in (0, 1)]
    pg.line(BX - 6.7, y2b(min(ys)) - 0.8, BX - 6.7, y2b(max(ys)) + 0.8, lw=0.4)
    pg.txt(BX - 7.3, y2b(np.mean(ys)), labt, ha="right")
ax.set_xlim(0, XMAX); ax.set_xticks([0, 10, 20, 30], []); ax.set_yticks([]); ax.spines["left"].set_position(("data", 0))
kend = key(pg, 115.6, 6.9, [(FOLD, "folded (F)", None), (DIRN, "directed (D)", None)], gap=1.1, sw=1.5)
pg.box(kend, 6.3, 1.5, 1.2, tint(FOLD, LIGHT["fold"]), z=2)
pg.box(kend + 1.7, 6.3, 1.5, 1.2, tint(DIRN, LIGHT["dir"]), z=2)
pg.txt(kend + 3.8, 6.9, "expected false")
pg.src("b", pd.DataFrame(out))

pg.panel(112, 34.6, "c", "Calls by internal r, closest-matched decoys, N0")
CY, CH = 38.3, 9.2
ax = pg.ax(BX, CY, BW, CH)
ax.set_ylim(3.95, -0.55)
y2c = lambda v: CY + CH * (v + 0.55) / 4.5  # noqa: E731
out = []
for gi, (sub, labt) in enumerate((("targets r <= 0.5", "r ≤ 0.5 (30)"), ("targets r > 0.5", "r > 0.5 (17 of 22)"))):
    q = SC[(SC.subset == sub) & (SC.alt == "N0-dir")]
    et = q[q.quantity == "exp_true"].iloc[0]
    fd = q[q.quantity == "fdr"].iloc[0]
    for si, (called, true_, fdr_, st, let) in enumerate(((et.ref_called, et.ref_value, fd.ref_value, "fold", "F"), (et.alt_called, et.alt_value, fd.alt_value, "dir", "D"))):
        y = gi * 2 + si + gi * 0.4
        stack(ax, y, called, called - true_, st)
        pg.txt(BX - 0.6, y2c(y), let, ha="right")
        row_values(y2c(y), called, called - true_, fdr_)
        out.append(dict(subset=labt, statistic=st, called=int(called), exp_true=true_, exp_false=called - true_, fdr=fdr_))
    pg.txt(BX - 3.1, y2c(gi * 2.4 + 0.5), labt, ha="right")
ax.set_xlim(0, XMAX); ax.set_xticks([0, 10, 20, 30]); ax.set_yticks([]); ax.spines["left"].set_position(("data", 0))
ax.set_xlabel("signatures called", labelpad=1)
pg.src("c", pd.DataFrame(out))

# ================================================================== d  each signature among its 320 co-expressed decoys (old h)
c = RC[RC.stat == "fold"]
t = c[c.kind == "target"].set_index("target")
dcx = c[(c.kind == "decoy") & (c["class"] == "coex")]
QD = dcx.groupby("target").med_delta.quantile([.05, .25, .5, .75, .95]).unstack()
ORDER = t.med_delta.sort_values(ascending=False).index.tolist()
RH = 1.9
Y0 = 62.6
DX, DW = 25.2, 42.8
pg.panel(1, 56, "d", "Signatures against their 320 co-expressed decoys")
ax = pg.ax(DX, Y0, DW, RH * 52)
out = []
for i, s in enumerate(ORDER):
    q = QD.loc[s].to_numpy()
    ax.plot([q[0], q[4]], [i, i], color=tint(COEX, 0.45), lw=0.6, solid_capstyle="butt", zorder=2)
    ax.plot([q[1], q[3]], [i, i], color=COEX, lw=1.5, solid_capstyle="butt", zorder=3)
    ax.plot([q[2], q[2]], [i - 0.36, i + 0.36], color="white", lw=0.55, solid_capstyle="butt", zorder=3.5)
    point(ax, t.loc[s, "med_delta"], i, SIG, ms=3.0, filled=bool(t.loc[s, "call"]), mew=0.55 if not t.loc[s, "call"] else 0.25, zorder=5)
    ym = Y0 + RH * (i + 0.5)
    name_right(pg, DX - 1.7, ym, s)
    pg.box(DX - 1.4, Y0 + RH * i + 0.08, 1.0, RH - 0.16, FAM[FAMS[s]], z=2)
    out.append(dict(signature=lab(s), family=FAMS[s], median_delta_auc=t.loc[s, "med_delta"], p=t.loc[s, "p"], called=bool(t.loc[s, "call"]),
                    decoy_p05=q[0], decoy_p25=q[1], decoy_median=q[2], decoy_p75=q[3], decoy_p95=q[4]))
zero(ax, 0)
ax.set_ylim(51.5, -0.5); ax.set_xlim(-0.034, 0.113); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([-0.025, 0, 0.025, 0.05, 0.075, 0.1], ["−0.025", "0", "0.025", "0.05", "0.075", "0.1"])
ax.set_xlabel("median ΔAUC, folded")
kx = 50.6
for k_y, it in zip(np.arange(133.4, 160, 2.45),
                   ((SIG, "called", "o"), (SIG, "not called", "open-o"), (COEX, "decoys 25–75%", "line"), (tint(COEX, 0.45), "decoys 5–95%", "line"),
                    (INFL, "inflamed/IFN", None), (CONF, "confounder", None), (OTH, "other", None))):
    key(pg, kx, k_y, [it], sw=1.8)
pg.src("d", pd.DataFrame(out))

# ---- q per signature, as fig4_v10.py computes it
QS = {}
PS = {}
NDEC = {}
for st in ("fold", "dir"):
    c_ = RC[RC.stat == st]
    t_ = c_[c_.kind == "target"].set_index("target")
    for cl in ("coex", "rand"):
        d_ = c_[(c_.kind == "decoy") & (c_["class"] == cl)]
        if st == "dir":
            d_ = d_[d_.dirn == d_.target.map(t_.dirn)]
        PS[(st, cl)] = d_.groupby("target").call.mean()
        NDEC[(st, cl)] = d_.groupby("target").size()
        if cl == "coex":
            s_t = np.sign(t_.med_delta) * -np.log10(t_.p.clip(lower=1e-300))
            s_d = np.sign(d_.med_delta) * -np.log10(d_.p.clip(lower=1e-300))
            QS[st] = pd.Series({tg: (1 + (s_d[d_.target == tg] >= s_t[tg] - 1e-12).sum()) / (1 + (d_.target == tg).sum()) for tg in t_.index})
PI0V = {st: RF[(RF.stat == st) & (RF.decoys == "coex") & (RF.boot == "clusters")].iloc[0].pi0 for st in ("fold", "dir")}

# ================================================================== e  aligned: the decoys called per signature and rank q (new)
EX, ECW = 70.2, 4.8
pg.panel(69, 56, "e", "Decoys called, rank q")
ECOLS = [("fold", "coex"), ("fold", "rand"), ("dir", "coex"), ("dir", "rand")]
for j, h in enumerate(("coex.", "rand.", "coex.", "rand.", "F", "D")):
    pg.txt(EX + ECW * j + ECW / 2, 61.3, h, ha="center")
for j0, h in ((0, "folded, %"), (2, "directed, %"), (4, "q")):
    pg.txt(EX + ECW * j0 + ECW, 59.0, h, ha="center", weight="bold")
    pg.line(EX + ECW * j0 + 0.3, 60.2, EX + ECW * j0 + 2 * ECW - 0.3, 60.2, lw=0.3)
out = []
for i, s in enumerate(ORDER):
    y = Y0 + RH * i
    row = dict(signature=lab(s))
    for j, (st, cl) in enumerate(ECOLS):
        v = 100 * PS[(st, cl)][s]
        tile(pg, EX + ECW * j, y, ECW, RH, PURPLE(min(v / 36, 1.0)), fmt(v))
        row[f"{st}_{cl}_called_pct"] = v
        row[f"{st}_{cl}_decoys"] = int(NDEC[(st, cl)][s])
    for j, st in enumerate(("fold", "dir")):
        v = QS[st][s]
        tile(pg, EX + ECW * (4 + j), y, ECW, RH, BLUE(v), f"{v:.2f}")
        row[f"q_{st}"] = v
    out.append(row)
ysum = Y0 + RH * 52 + 1.0
for j, (st, cl) in enumerate(ECOLS):
    pg.txt(EX + ECW * j + ECW / 2, ysum, f"{PS[(st, cl)].sum():.1f}", ha="center", weight="bold")
for j, st in enumerate(("fold", "dir")):
    pg.txt(EX + ECW * (4 + j) + ECW / 2, ysum, f"{PI0V[st]:.2f}", ha="center", weight="bold")
pg.txt(EX + 2 * ECW, ysum + 2.1, "expected false", ha="center")
pg.txt(EX + 5 * ECW, ysum + 2.1, "π0", ha="center")
pg.src("e", pd.DataFrame(out))

# ================================================================== f  aligned: calls under each null and statistic, first decoy set analysis (new)
FX, FCW, FGAP = 100.5, 2.28, 0.4
pg.panel(100, 56, "f", "Calls by null")
TC = U[U.kind == "target"].set_index(["target", "null", "stat"])
fx = {}
for i, nl in enumerate(order):
    for j, (st, let) in enumerate((("fold", "F"), ("dir", "D"))):
        fx[(nl, st)] = FX + FCW * (2 * i + j) + FGAP * i + FCW / 2
        pg.txt(fx[(nl, st)], 61.3, let, ha="center")
    pg.txt(fx[(nl, "fold")] + FCW / 2, 59.0, NAMES[nl], ha="center", weight="bold")
    pg.line(fx[(nl, "fold")] - FCW / 2 + 0.15, 60.2, fx[(nl, "dir")] + FCW / 2 - 0.15, 60.2, lw=0.3)
out = []
for i, s in enumerate(ORDER):
    y = Y0 + RH * (i + 0.5)
    for k, xx in fx.items():
        r = TC.loc[(s, k[0], k[1])]
        if bool(r.call):
            pg.c.plot(xx, y, "o", ms=2.1, mfc=SIG, mec="white", mew=0.2, zorder=4)
        else:
            pg.c.plot(xx, y, "o", ms=1.5, mfc="white", mec=SIG, mew=0.3, zorder=4)
        out.append(dict(signature=lab(s), null=k[0], statistic=k[1], called=bool(r.call), median_delta=r.med_delta, p=r.p))
for k, xx in fx.items():
    pg.txt(xx, ysum, f"{int(U[(U.kind == 'target') & (U['null'] == k[0]) & (U.stat == k[1])].call.sum())}", ha="center", weight="bold")
for g0, g1, labt in (("N1", "N2", "enriched"), ("N1nat", "N2nat", "natural")):
    pg.txt((fx[(g0, "fold")] + fx[(g1, "dir")]) / 2, ysum + 2.1, labt, ha="center")
pg.src("f", pd.DataFrame(out))

# ================================================================== g  ten decoy sets (old e)
pg.panel(127.5, 56, "g", "Ten decoy sets")
ax = pg.ax(134.8, 62, 17.5, 29)
F0 = F[(F["null"] == "N0") & (F.pool == "A") & (F.decoys == "coex")].set_index("stat")
out = []
for gi, st in enumerate(("fold", "dir")):
    b = RB[RB.stat == st].fdr.to_numpy()
    pooled = RF[(RF.stat == st) & (RF.decoys == "coex") & (RF.boot == "clusters")].iloc[0]
    xs = gi + np.linspace(-0.2, 0.2, len(b))
    point(ax, xs, b, EST, ms=2.5)
    ax.plot([gi - 0.27, gi + 0.27], [pooled.fdr] * 2, color=INK, lw=0.8, zorder=5, solid_capstyle="butt")
    point(ax, gi + 0.37, F0.loc[st, "fdr"], INK, ms=2.8, filled=False)
    ax.text(gi, b.max() + 1.0, f"{pooled.fdr:.0f}%", ha="center", va="bottom", fontsize=5)
    out += [dict(statistic=st, source=f"set {k + 1}", fdr=v) for k, v in enumerate(b)]
    out += [dict(statistic=st, source="pooled", fdr=pooled.fdr), dict(statistic=st, source="first set, calibration", fdr=F0.loc[st, "fdr"])]
ax.set_xticks([0, 1], ["folded", "directed"]); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.5, 1.55); ax.set_ylim(0, 24); ax.set_yticks([0, 5, 10, 15, 20]); ax.set_ylabel("estimated FDR (%)")
key(pg, 142.8, 63.2, [(EST, "each set", "o")], sw=1.5)
pg.line(142.8, 65.7, 144.3, 65.7, lw=0.8)
pg.txt(144.9, 65.7, "pooled")
pg.c.plot(143.55, 68.2, "o", ms=2.8, mfc="white", mec=INK, mew=0.6, zorder=3)
pg.txt(144.9, 68.2, "first set")
pg.src("g", pd.DataFrame(out))

# ================================================================== h  per test and per module (old c)
RNG = np.random.default_rng(20261019)


def boot_rate(df, col="call", nb=4000):
    per = df.groupby("target")[col].mean().to_numpy(float)
    bs = [per[RNG.integers(0, len(per), len(per))].mean() for _ in range(nb)]
    return 100 * per.mean(), 100 * np.percentile(bs, 2.5), 100 * np.percentile(bs, 97.5)


pg.panel(153.5, 56, "h", "Per test, per module")
ax = pg.ax(160, 62, 18.8, 29)
xd = X[(X["null"] == "N0") & (X.kind == "decoy") & (X.pool == "A")].copy()
xd["clear"] = (xd.fold_obs > xd.fold_p95).astype(float)
out = []
for gi, (g, nom) in enumerate((("test", 5.0), ("module", 2.5))):
    ax.plot([gi - 0.45, gi + 0.45], [nom, nom], color=NOM, lw=0.7, ls=(0, (1.0, 1.5)), zorder=1)
    ax.text(gi, nom - 0.2, f"{nom:g}%", fontsize=5, ha="center", va="top")
    for si, (s, cc, mk) in enumerate((("rand", RAND, "D"), ("coex", COEX, "s"))):
        v = boot_rate(xd[xd.structure == s], "clear") if g == "test" else boot_rate(
            U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == s) & (U["null"] == "N0") & (U.stat == "fold")])
        x = gi + (si - 0.5) * 0.48
        ci(ax, v[1], v[2], x, cc, lw=0.7, horizontal=False)
        point(ax, x, v[0], cc, mk, ms=3.0 if mk == "D" else 3.4)
        ax.text(x, v[2] + 0.25, f"{v[0]:.1f}", ha="center", va="bottom", fontsize=5)
        out.append(dict(level=g, decoys=s, called_pct=v[0], lo=v[1], hi=v[2], nominal=nom))
ax.set_xticks([0, 1], ["each\ntest", "each\nmodule"], fontsize=5.5, linespacing=1.0); ax.tick_params(axis="x", length=0)
ax.set_ylabel("decoys called (%)"); ax.set_ylim(0, 13); ax.set_yticks([0, 4, 8, 12]); ax.set_xlim(-0.55, 1.55)
key(pg, 160.8, 63.2, [(COEX, "co-expressed", "s")], sw=1.5)
key(pg, 160.8, 65.7, [(RAND, "random-gene", "D")], sw=1.5)
pg.src("h", pd.DataFrame(out))

# ================================================================== i  every decoy set, values (new)
pg.panel(127.5, 98.2, "i", "Each decoy set")
IX, ICW, IY, ICH = 134.5, 5.55, 105.4, 2.3
ICOLS = [("decoy_rate", "rate", lambda v: BLUE(v / 12), lambda v: f"{v:.1f}"), ("fdr", "FDR", lambda v: BLUE(v / 24), lambda v: f"{v:.1f}"),
         ("pi0", "π0", lambda v: PURPLE(v), lambda v: f"{v:.2f}"), ("fdr_pi0", "corr.", lambda v: PURPLE(v / 12), lambda v: f"{v:.1f}")]
for g, (st, name) in enumerate(ST):
    x_g = IX + ICW * 4 * g
    pg.txt(x_g + 2 * ICW, 101.8, name, ha="center", weight="bold")
    pg.line(x_g + 0.3, 103.0, x_g + 4 * ICW - 0.5, 103.0, lw=0.3)
    for j, (_, h, _f, _g) in enumerate(ICOLS):
        pg.txt(x_g + ICW * j + ICW / 2, 104.2, h, ha="center")
for k in range(1, 11):
    pg.txt(IX - 0.7, IY + ICH * (k - 1) + ICH / 2, str(k), ha="right")
    for g, (st, _) in enumerate(ST):
        r = RB[(RB.set == k) & (RB.stat == st)].iloc[0]
        for j, (c_, _h, cm, f_) in enumerate(ICOLS):
            tile(pg, IX + ICW * (4 * g + j), IY + ICH * (k - 1), ICW, ICH, cm(r[c_]), f_(r[c_]))
pg.txt(IX - 0.7, IY + ICH * 10 + ICH / 2 + 0.2, "pooled", ha="right", weight="bold")
for g, (st, _) in enumerate(ST):
    r = RF[(RF.stat == st) & (RF.decoys == "coex") & (RF.boot == "clusters")].iloc[0]
    for j, (c_, _h, cm, f_) in enumerate(ICOLS):
        pg.txt(IX + ICW * (4 * g + j) + ICW / 2, IY + ICH * 10 + ICH / 2 + 0.2, f_(r[c_]), ha="center", weight="bold")
pg.src("i", pd.concat([RB, RF[(RF.decoys == "coex") & (RF.boot == "clusters")][["stat", "called", "decoy_rate", "exp_false", "fdr", "pi0", "fdr_pi0"]].assign(set="pooled")]))

# ================================================================== j, k  rank among own decoys (old g, and the directed q)
for lt, xa, st, name in (("j", 134.5, "fold", "folded"), ("k", 160.4, "dir", "directed")):
    pg.panel(xa - 7, 134.8, lt, f"q, {name}")
    ax = pg.ax(xa, 140.2, 18.5, 22)
    bins = np.linspace(0, 1, 11)
    hist(ax, QS[st], bins)
    ax.axhline(PI0V[st] * 52 / 10, color=PI0, lw=0.7, zorder=4)
    ax.text(0.98, 35.5, f"π0 = {PI0V[st]:.2f}", ha="right", va="top", fontsize=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 36); ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, 10, 20, 30])
    ax.set_xlabel("rank q"); ax.set_ylabel("signatures")
    pg.src(lt, pd.DataFrame({"signature": [lab(s) for s in QS[st].index], "q": QS[st].values}))

# ================================================================== l  matching (old a)
pg.panel(1, 171.6, "l", "Internal r, signature and decoys")
ax = pg.ax(9.3, 177.4, 27, 21.2)
x0 = X[X["null"] == "N0"]
tr = x0[x0.kind == "target"].set_index(["target", "cohort"]).r_int
out = []


def regions(ax, x, y, color, levels=(0.5, 0.75, 0.9, 0.97), n=120, bw=0.18):
    """highest-density regions drawn as hairlines, as fig4_v10.py computes them"""
    from scipy.stats import gaussian_kde
    k = gaussian_kde(np.vstack([x, y]), bw_method=bw)
    gx, gy = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    dz = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    s = np.sort(dz.ravel())[::-1]
    cs = np.cumsum(s) / s.sum()
    th = sorted(s[np.searchsorted(cs, lv)] for lv in levels)
    ax.contour(gx, gy, dz, levels=th, colors=[color], linewidths=0.4, zorder=3)


for s_, cc, mk in (("rand", RAND, "D"), ("coex", COEX, "s")):
    dd = x0[(x0.kind == "decoy") & (x0.pool == "A") & (x0.structure == s_)].groupby(["target", "cohort"]).r_int.median()
    j_ = pd.concat([tr, dd], axis=1, keys=["t", "d"]).dropna()
    ax.plot(j_.t, j_.d, mk, ms=0.9, mfc=cc, mec="none", ls="none", zorder=2, alpha=0.8)
    regions(ax, j_.t.to_numpy(), j_.d.to_numpy(), cc)
    out.append(j_.reset_index().assign(decoys=s_))
ax.plot([0, 1], [0, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_xlabel("signature internal r"); ax.set_ylabel("decoys' median r")
key(pg, 10.2, 178.4, [(COEX, "co-expressed", "s")], sw=1.5)
key(pg, 10.2, 180.9, [(RAND, "random-gene", "D")], sw=1.5)
pg.src("l", pd.concat(out).rename(columns={"t": "signature_r", "d": "decoy_median_r"}))

# ================================================================== m  how close the decoys sit, all ten sets (new)
pg.panel(41, 171.6, "m", "Decoy r minus signature r")
ax = pg.ax(49.4, 177.4, 27, 21.2)
RI["dev"] = RI.r_achieved_med - RI.r_target_med
bins = np.arange(-0.75, 0.4001, 0.05)
dcoex = RI[RI["class"] == "coex"].dev
drand = RI[RI["class"] == "rand"].dev
hist(ax, dcoex, bins, weights=100 * np.ones(len(dcoex)) / len(dcoex))
hr, _ = np.histogram(drand, bins)
ax.step(np.r_[bins[:-1], bins[-1]], np.r_[100 * hr / len(drand), 0], where="post", color=RAND, lw=0.7, zorder=4)
for v in (-0.05, 0.05):
    ax.plot([v, v], [0, 40], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=5)
ax.text(0, 40.5, "±0.05", ha="center", va="bottom", fontsize=5)
ax.set_xlim(-0.75, 0.4); ax.set_xticks([-0.5, 0, 0.25], ["−0.5", "0", "0.25"]); ax.set_ylim(0, 60); ax.set_yticks([0, 20, 40])
ax.set_xlabel("decoy r − signature r"); ax.set_ylabel("decoys (%)")
nc, nr = int((dcoex.abs() <= 0.05).sum()), int((drand.abs() <= 0.05).sum())
pg.box(50.2, 177.8, 1.5, 1.2, NULLC, z=3)
pg.txt(52.3, 178.4, f"co-expressed, {nc:,} of {len(dcoex):,}")
pg.line(50.2, 180.9, 51.7, 180.9, color=RAND, lw=0.7)
pg.txt(52.3, 180.9, f"random-gene, {nr:,} of {len(drand):,}")
pg.src("m", pd.DataFrame({"bin_left": np.round(bins[:-1], 2), "coex_pct": np.histogram(dcoex, bins)[0] * 100 / len(dcoex), "rand_pct": 100 * hr / len(drand)}))

# ================================================================== n  each cohort, tests above the null (new)
pg.panel(81, 171.6, "n", "Tests above null P95, by cohort")
ax = pg.ax(93.6, 177.4, 22, 21.2)
byc = xd.groupby(["cohort", "structure"]).clear.mean().unstack() * 100
byc["n"] = xd.groupby("cohort").n.first()
byc = byc.sort_values("n")
for i, (co, r) in enumerate(byc.iterrows()):
    point(ax, r.rand, i, RAND, "D", ms=2.1)
    point(ax, r.coex, i, COEX, "s", ms=2.3)
    ax.text(13.6, i, f"{int(r.n)}", ha="right", va="center", fontsize=5, clip_on=False)
ax.text(13.6, len(byc) - 0.25, "n", ha="right", va="bottom", fontsize=5, weight="bold", clip_on=False)
ax.axvline(5, color=NOM, lw=0.7, ls=(0, (1.0, 1.5)), zorder=1)
ax.text(5, len(byc) - 0.25, "5%", ha="center", va="bottom", fontsize=5, clip_on=False)
ax.set_yticks(range(len(byc)), [c_.replace("ICB_", "").replace("_", " ") for c_ in byc.index], fontsize=5)
ax.set_ylim(-0.6, len(byc) - 0.4); ax.set_xlim(0, 11.5); ax.set_xticks([0, 5, 10]); ax.set_xlabel("decoy tests above (%)")
pg.src("n", byc.reset_index().rename(columns={"coex": "coexpressed_pct", "rand": "random_gene_pct", "n": "patients"}))

# ================================================================== o  decoy call rate by own internal r (old k, with the directed rows)
pg.panel(126.5, 171.6, "o", "Decoys called, by own internal r")
OX, OCW, OY, OCH = 137.5, 6.9, 180.3, 2.45
bins_ = list(dict.fromkeys(RCV.rbin))
RMAX = RCV[RCV["null"].isin(["N0", "N1nat", "N2nat"])].call_rate.max()
hi_ = 12.0
assert RMAX <= hi_


def rate_colour(v):
    return PUOR(0.5 + (v - 2.5) / (2 * (hi_ - 2.5)))


pg.txt(OX + 3 * OCW, 175.5, "decoy internal r", ha="center")
for j, b in enumerate(bins_):
    pg.txt(OX + OCW * j + OCW / 2, 177.5, b.replace("(", "").replace("]", "").replace(", ", "–"), ha="center")
    pg.txt(OX + OCW * j + OCW / 2, 179.4, f"{int(RCV[(RCV['null'] == 'N0') & (RCV.stat == 'fold') & (RCV.rbin == b)].decoys.iloc[0])}", ha="center", color=MUTED)
pg.txt(OX - 0.7, 179.4, "n", ha="right", color=MUTED)
for g, (st, name) in enumerate(ST):
    for i, (nl, labt) in enumerate((("N0", "N0"), ("N1nat", "N1"), ("N2nat", "N2"))):
        y = OY + OCH * (3 * g + i) + (0.5 if g else 0)
        pg.txt(OX - 0.7, y + OCH / 2, labt, ha="right")
        for j, b in enumerate(bins_):
            v = RCV[(RCV["null"] == nl) & (RCV.stat == st) & (RCV.rbin == b)].call_rate.iloc[0]
            tile(pg, OX + OCW * j, y, OCW, OCH, rate_colour(v), f"{v:.1f}")
    pg.txt(126.8, OY + OCH * (3 * g + 1.5) + (0.5 if g else 0), name, ha="left")
a_ = pg.ax(OX, 196.9, 22, 1.4)
vv = np.linspace(0, hi_, 200)
a_.imshow(np.array([[rate_colour(v) for v in vv]]), aspect="auto", extent=(0, hi_, 0, 1), rasterized=True)
a_.set_yticks([]); a_.set_xticks([0, 2.5, 6, 12], ["0", "2.5%", "6", "12"], fontsize=5); a_.tick_params(length=1.5, pad=1, width=0.4)
for s_ in a_.spines.values():
    s_.set_visible(False)
pg.txt(OX + 24.6, 197.6, "decoys called")
pg.txt(OX + 24.6, 200.6, "N1, N2 natural pool")
pg.src("o", RCV[RCV["null"].isin(["N0", "N1nat", "N2nat"])])

# ================================================================== p, q  the threshold (old i, split)
for lt, xa, title, col, ylab, lim, ticks in (("p", 9.3, "Calls, by P threshold", "called", "signatures called", 40, [0, 20, 40]),
                                             ("q", 54.3, "FDR, by P threshold", "fdr", "estimated FDR (%)", 40, [0, 20, 40])):
    pg.panel(xa - 8.3, 208, lt, title)
    ax = pg.ax(xa, 213.6, 31, 18)
    for st in ("fold", "dir"):
        ax.plot(TH.threshold, TH[f"{col}_{st}"], color=STAT[st], lw=0.8)
    ax.set_xscale("log"); ax.set_xlim(0.001, 0.2); ax.minorticks_off(); zero(ax, 0.05)
    ax.set_xticks([0.001, 0.01, 0.05, 0.2], ["0.001", "0.01", "0.05", "0.2"]); ax.set_ylim(0, lim); ax.set_yticks(ticks)
    ax.set_xlabel("P threshold"); ax.set_ylabel(ylab)
    key(pg, xa + 1.0, 214.6, [(FOLD, "folded", "line")], sw=1.8)
    key(pg, xa + 1.0, 217.1, [(DIRN, "directed", "line")], sw=1.8)
    pg.src(lt, TH[["threshold", f"{col}_fold", f"{col}_dir"]])

# ================================================================== r  QQ (old j)
pg.panel(92, 208, "r", "P values against uniform")
ax = pg.ax(100.2, 213.6, 24, 18)
c_ = RC[RC.stat == "fold"]
out = []
for msk, cc, mk, s_, lab_, ms in (((c_.kind == "decoy") & (c_["class"] == "rand"), RAND, "D", False, "random-gene", 0.95),
                                  ((c_.kind == "decoy") & (c_["class"] == "coex"), COEX, "s", False, "co-expressed", 0.95),
                                  ((c_.kind == "target"), SIG, "o", True, "signatures", 2.2)):
    pv = np.sort(c_[msk].p.to_numpy())
    n_ = len(pv)
    ex = (np.arange(1, n_ + 1) - 0.5) / n_
    ax.plot(-np.log10(ex), -np.log10(pv), mk, ms=ms, mfc=cc, mec="white" if s_ else "none", mew=0.2 if s_ else 0, ls="none", zorder=4 if s_ else 3 if cc == COEX else 2, rasterized=n_ > 2000)
    out.append(pd.DataFrame({"unit": lab_, "expected_log10": -np.log10(ex), "observed_log10": -np.log10(pv)}))
ax.plot([0, 4.5], [0, 4.5], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.set_xlim(0, 4.5); ax.set_ylim(0, 4.5); ax.set_xticks([0, 2, 4]); ax.set_yticks([0, 2, 4])
ax.set_xlabel("expected −log10 P"); ax.set_ylabel("observed −log10 P")
for yk, it in zip((214.6, 217.1, 219.6), ((SIG, "signatures", "o"), (COEX, "co-expressed", "s"), (RAND, "random-gene", "D"))):
    key(pg, 125.2, yk, [it], sw=1.5)
pg.src("r", pd.concat(out))

# ================================================================== s  pooling decoy sets (old l)
pg.panel(139.5, 208, "s", "Pooling decoy sets")
ax = pg.ax(147.8, 213.6, 31, 18)
for st in ("fold", "dir"):
    ax.plot(CONV.sets_pooled, CONV[f"fdr_{st}"], color=STAT[st], lw=0.6, zorder=3)
    point(ax, CONV.sets_pooled, CONV[f"fdr_{st}"], STAT[st], ms=2.4)
ax.set_xticks([1, 5, 10]); ax.set_xlim(0.5, 10.5); ax.set_ylim(0, 25); ax.set_yticks([0, 10, 20])
ax.set_xlabel("decoy sets pooled"); ax.set_ylabel("estimated FDR (%)")
key(pg, 161, 214.6, [(FOLD, "folded", "o"), (DIRN, "directed", "o")], sw=1.5, gap=1.5)
pg.src("s", CONV)

pg.save("Fig4")
