"""Figure 7 redrawn at the density and in the look of the AML multi-omic figures (pixl_nc.py).
The breast-cancer prognostic signatures of Venet et al. in 13 cohorts (folded C-index, first principal component unless marked,
ten decoy sets). Tables: icb/v12_pkg/breast_* (applications/v12_rerun.R, v12_breast_origin.R), icb/v11_pkg/breast_pc1_strength.csv
(v11_breast_pc1_strength.R, post hoc), icb/v13_fig/km_breast.csv (scripts/figprep_km.R), breast/cohorts/manifest.csv,
validation/venet_provenance*.csv. New panels (j, m, n, o, s) and the random-gene rows of c are reads of those tables, descriptive
(validation/NOTE_figures_v17_descriptive_only.md). Nothing here is resampled or tested anew; the log-rank P values are the two the old
sheet printed, recomputed by the same formula as cns.logrank_p and checked against the old sheet."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt  # noqa: F401
from matplotlib.colors import to_hex

P = D + "/icb/v12_pkg/"
V = D + "/signature-decoy-fdr/validation/"
F = {k: pd.read_csv(P + f"breast_{k}_fdr.csv") for k in ("pc1", "mean", "abs_pc1")}
OR = pd.read_csv(P + "breast_origin_fdr.csv")
B = {k: pd.read_csv(P + f"breast_{k}_fdr_by_set.csv") for k in ("pc1", "mean", "abs_pc1")}
CL = pd.read_csv(P + "breast_pc1_calls.csv.gz")
CO = pd.read_csv(P + "breast_pc1_origin_calls.csv.gz")
UA = pd.read_csv(P + "breast_pc1_units.csv.gz")
T = UA[UA.kind == "target"].copy()
T["dC"] = T.fold_obs - T.fold_med
GN = pd.read_csv(P + "breast_pc1_globalnull_summary.csv").iloc[0]
GP = pd.read_csv(P + "breast_pc1_globalnull_per_perm.csv")
GT = pd.read_csv(P + "breast_pc1_globalnull_per_target.csv")
PIL = pd.read_csv(P + "breast_pc1_pi0_lambda.csv")
INFO = pd.read_csv(P + "breast_decoy_info.csv.gz", usecols=["decoy", "target", "class", "set", "k", "r_target_med", "r_achieved_med"])
ST = pd.read_csv(D + "/icb/v11_pkg/breast_pc1_strength.csv")
EXC = pd.read_csv(P + "breast_origin_excluded_tests.csv")
MAN = pd.read_csv(D + "/breast/cohorts/manifest.csv")
MAN = MAN[MAN.included.astype(str) == "True"].set_index("cohort")
PV = pd.read_csv(V + "venet_provenance.csv").set_index("signature")
BL = pd.read_csv(V + "venet_provenance_borderline.csv")
KM = pd.read_csv(D + "/icb/v13_fig/km_breast.csv")
COH = ["NKI", "Loi_GSE6532", "TCGA_BRCA", "GSE7390", "GSE11121", "GSE1456", "GSE4922", "GSE20685", "GSE21653", "GSE42568",
       "GSE17705", "GSE96058", "METABRIC"]
CLAB = {"Loi_GSE6532": "Loi", "TCGA_BRCA": "TCGA", "GSE96058": "SCAN-B"}
EXS = "SOTIRIOU-GGI"


def cl_(c):
    return CLAB.get(c, c)


def row(df, decoys="coex", boot="clusters"):
    return df[(df.stat == "fold") & (df.decoys == decoys) & (df.boot == boot)].iloc[0]


def orow(score, rule, boot="clusters"):
    return OR[(OR.score == score) & (OR["rule"] == rule) & (OR.boot == boot)].iloc[0]


def ref(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=3)


def km_curve(time, event):
    t, e = np.asarray(time, float), np.asarray(event, int)
    s, xs, ys = 1.0, [0.0], [1.0]
    for u in np.unique(t[e == 1]):
        s *= 1 - ((t == u) & (e == 1)).sum() / (t >= u).sum()
        xs.append(u)
        ys.append(s)
    return np.array(xs), np.array(ys)


def logrank_p(time, event, group):
    """two-group log-rank P, as cns.logrank_p; it reproduces the two values on the old sheet"""
    from scipy.stats import chi2
    t, e, g = np.asarray(time, float), np.asarray(event, int), np.asarray(group, bool)
    oe, v = 0.0, 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u
        n, n1 = at.sum(), (at & g).sum()
        d, d1 = ((t == u) & (e == 1)).sum(), ((t == u) & (e == 1) & g).sum()
        if n > 1:
            oe += d1 - d * n1 / n
            v += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    return float(chi2.sf(oe ** 2 / v, 1))


def swarm(v, dv, dx, cap):
    """offsets that stack points of nearly equal value side by side, centre out; a stack wider than cap wraps back inside it"""
    v = np.asarray(v, float)
    idx = np.floor((v - v.min()) / dv).astype(int)
    off = np.zeros(len(v))
    for b in np.unique(idx):
        k = np.where(idx == b)[0]
        k = k[np.argsort(v[k])]
        off[k] = np.array([0] + [q_ for i_ in range(1, len(k)) for q_ in (i_, -i_)])[: len(k)] * dx
    return (off + cap) % (2 * cap) - cap


def width(s, size=5.0, weight="normal"):
    return TextPath((0, 0), s, prop=fm.FontProperties(family="Helvetica", size=size, weight=weight)).get_extents().width * 25.4 / 72


pg = Page(240)

# ---------------------------------------------------------------- a  the null bank by cohort; rows shared with b
Y1, RH1, N1 = 9.0, 2.5, 13
order = T.groupby("cohort").fold_med.median().sort_values(ascending=False).index
pg.panel(1, 3, "a", "Null median folded C")
ax = pg.ax(12, Y1, 28, RH1 * N1)
out = []
for i, co in enumerate(order):
    v = T[T.cohort == co].fold_med.to_numpy()
    q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
    point(ax, v, i + swarm(v, 0.003, 0.2, 0.42), NULLC, ms=1.5, mew=0)
    ci(ax, q[0], q[4], i, MUTED, lw=0.5, zorder=5)
    ci(ax, q[1], q[3], i, INK, lw=1.2, zorder=5)
    ax.plot([q[2], q[2]], [i - 0.36, i + 0.36], color=INK, lw=0.8, zorder=6, solid_capstyle="butt")
    pg.txt(47.5, Y1 + RH1 * (i + 0.5), f"{np.median(v):.3f}", ha="right")
    out.append((cl_(co), *q))
zero(ax, 0.5)
ax.set_ylim(N1 - 0.5, -0.5)
ax.set_yticks(range(N1), [cl_(c) for c in order], fontsize=5)
ax.set_xlim(0.495, 0.66)
ax.set_xticks([0.5, 0.55, 0.6, 0.65], ["0.50", "0.55", "0.60", "0.65"])
ax.set_xlabel("null median folded C")
pg.txt(47.5, 7.3, "median", ha="right", weight="bold")
for i in range(1, N1):
    pg.line(41, Y1 + RH1 * i, 94.5, Y1 + RH1 * i, color="#e8e8e8", lw=0.25)
pg.src("a", T[["cohort", "target", "fold_med"]].assign(cohort=T.cohort.map(cl_)))

# ---------------------------------------------------------------- b  patients and events, same rows
pg.panel(50, 3, "b", "Patients and events, 13 cohorts")
ax = pg.ax(52, Y1, 18, RH1 * N1)
mo = MAN.loc[list(order)]
for i, (co, r) in enumerate(mo.iterrows()):
    point(ax, r.n_samples, i, INK, ms=2.7)
    point(ax, r.n_events, i, OTHER, ms=2.7)
    yy = Y1 + RH1 * (i + 0.5)
    pg.txt(79, yy, f"{int(r.n_samples):,}", ha="right")
    pg.txt(85.5, yy, f"{int(r.n_events):,}", ha="right")
    pg.txt(87, yy, r.endpoint)
ax.set_xscale("log")
ax.set_xlim(20, 5000)
ax.set_xticks([30, 300, 3000], ["30", "300", "3,000"])
ax.minorticks_off()
ax.set_ylim(N1 - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("count (log scale)")
key(pg, 51.5, 7.3, [(INK, "patients", "o"), (OTHER, "events", "o")], gap=1.2)
pg.txt(79, 7.3, "patients", ha="right", weight="bold")
pg.txt(85.5, 7.3, "events", ha="right", weight="bold")
pg.txt(87, 7.3, "endpoint", weight="bold")
ytot = Y1 + RH1 * N1
pg.line(72.5, ytot, 94.5, ytot, color=INK, lw=0.3)
pg.txt(79, ytot + 1.6, f"{int(mo.n_samples.sum()):,}", ha="right")
pg.txt(85.5, ytot + 1.6, f"{int(mo.n_events.sum()):,}", ha="right")
pg.txt(87, ytot + 1.6, "all")
pg.src("b", mo[["n_samples", "n_events", "endpoint"]].reset_index().assign(cohort=lambda d_: d_.cohort.map(cl_)))

# ---------------------------------------------------------------- c  the six old estimates, every value printed, and the random-gene decoys
pg.panel(97, 3, "c", "Calls and the estimate")
rows = [("PC1, every cohort", row(F["pc1"]), "coex"), ("PC1, development out", orow("pc1", "D and C"), "coex"),
        ("PC1, + borderline out", orow("pc1", "D, C and every borderline"), "coex"), ("PC1, closest decoys", row(F["pc1"], "coex closest"), "coex"),
        ("PC1, abs. r decoys (post hoc)", row(F["abs_pc1"]), "coex"), ("PC1, random decoys", row(F["pc1"], "rand"), "rand"),
        ("mean z", row(F["mean"]), "coex"), ("mean z, random decoys", row(F["mean"], "rand"), "rand")]
XMAX = 70
HC = RH1 * N1
ax = pg.ax(122.5, Y1, 18.5, HC)
ysc = [k + (0.55 if k >= 6 else 0) for k in range(len(rows))]
out = []
for y, (lab, r, u) in zip(ysc, rows):
    ci(ax, r.fdr_lo, min(r.fdr_hi, XMAX), y, EST, lw=0.7)
    if r.fdr_hi > XMAX:
        ax.plot(XMAX, y, marker=">", ms=2.2, mfc=EST, mec="none", clip_on=False, zorder=4)
    point(ax, r.fdr_pi0, y, PI0, ms=2.9, filled=False, mew=0.5)
    point(ax, r.fdr, y, EST, ms=3.3)
    yy = Y1 + HC * (y + 0.6) / (max(ysc) + 1.2)
    pg.txt(147, yy, f"{int(r.called)}", ha="right")
    pg.txt(163.2, yy, f"{r.fdr:.1f} ({r.fdr_lo:.1f}–{r.fdr_hi:.1f})", ha="right")
    pg.txt(179.6, yy, f"{r.fdr_pi0:.1f} ({r.fdr_pi0_lo:.1f}–{r.fdr_pi0_hi:.1f})", ha="right")
    out.append((lab, u, int(r.called), r.exp_false, r.fdr, r.fdr_lo, r.fdr_hi, r.pi0, r.fdr_pi0, r.fdr_pi0_lo, r.fdr_pi0_hi))
ax.set_ylim(max(ysc) + 0.6, -0.6)
ax.set_yticks(ysc, [x[0] for x in rows], fontsize=5)
ax.set_xlim(0, XMAX)
ax.set_xticks([0, 20, 40, 60])
ax.set_xlabel("estimated FDR (%)")
pg.txt(147, 7.3, "called", ha="right", weight="bold")
pg.txt(163.2, 7.3, "estimate", ha="right", weight="bold")
pg.txt(179.6, 7.3, "π0-corrected", ha="right", weight="bold")
key(pg, 110.5, 7.3, [(EST, "estimate", "o"), (PI0, "π0-corrected", "open-o")], sw=1.4, gap=1.4)
pg.src("c", pd.DataFrame(out, columns=["analysis", "decoys", "called", "exp_false", "fdr", "fdr_lo", "fdr_hi", "pi0", "fdr_pi0", "fdr_pi0_lo", "fdr_pi0_hi"]))

# ---------------------------------------------------------------- d  development cohorts out
Y2, H2 = 55.5, 25
pg.panel(1, 50.5, "d", "Development cohorts out")
a_ = CL[(CL.stat == "fold") & (CL.kind == "target")].set_index("target")
b_ = CO[(CO.stat == "fold") & (CO.kind == "target")].set_index("target")
aff = sorted(EXC.target.unique(), key=lambda s: -a_.loc[s, "med_delta"])
ax = pg.ax(22.5, Y2, 18.5, H2)
RHD = H2 / len(aff)
out = []
for i, s in enumerate(aff):
    m0, m1 = a_.loc[s, "med_delta"], b_.loc[s, "med_delta"]
    ax.plot([m0, m1], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, m0, i, INK, ms=2.8, filled=bool(a_.loc[s, "call"]))
    point(ax, m1, i, SIG, ms=2.8, filled=bool(b_.loc[s, "call"]), zorder=5)
    yy = Y2 + RHD * (i + 0.5)
    pg.txt(47.5, yy, minus(f"{m0:.3f}"), ha="right")
    pg.txt(53.5, yy, minus(f"{m1:.3f}"), ha="right")
    out.append((s, int((EXC.target == s).sum()), m0, bool(a_.loc[s, "call"]), m1, bool(b_.loc[s, "call"])))
zero(ax, 0)
ax.set_ylim(len(aff) - 0.5, -0.5)
ax.set_yticks(range(len(aff)), [f"{s} ({(EXC.target == s).sum()})" for s in aff], fontsize=5)
ax.set_xlim(-0.02, 0.095)
ax.set_xticks([0, 0.05], ["0.00", "0.05"])
ax.set_xlabel("median ΔC")
pg.txt(47.5, Y2 - 1.5, "all", ha="right", weight="bold")
pg.txt(53.5, Y2 - 1.5, "out", ha="right", weight="bold")
key(pg, 29, Y2 + RHD * 8.5, [(INK, "all cohorts", "o")])
key(pg, 29, Y2 + RHD * 9.8, [(SIG, "cohorts out", "o")])
pg.src("d", pd.DataFrame(out, columns=["signature", "tests_left_out", "median_dC_all", "called_all", "median_dC_out", "called_out"]))

# ---------------------------------------------------------------- e, f  variance on PC1, one panel per decoy class
LIMV = 0.8
for lt, xa, col, unit, m, name, yl in (("e", 64.5, "pc1_coex_med", "coex", "s", "PC1 variance, co-expr.", "co-expressed (median)"),
                                       ("f", 96.5, "pc1_rand_med", "rand", "D", "PC1 variance, random", "random-gene (median)")):
    pg.panel(xa - 7.5, 50.5, lt, name)
    ax = pg.ax(xa, Y2, 21.5, H2)
    ax.hexbin(ST.pc1_target, ST[col], gridsize=15, cmap=GREYHEX, mincnt=1, linewidths=0, rasterized=True, extent=(0, LIMV, 0, LIMV))
    ref(ax, 0, LIMV)
    point(ax, ST.pc1_target.median(), ST[col].median(), UNIT[unit], m, ms=3.2, zorder=6)
    ax.set_xlim(0, LIMV)
    ax.set_ylim(0, LIMV)
    ax.set_xticks([0, 0.4, 0.8], ["0", "0.4", "0.8"])
    ax.set_yticks([0, 0.4, 0.8], ["0", "0.4", "0.8"])
    ax.set_xlabel("signature")
    ax.set_ylabel(yl)
    ax.text(0.03, 0.78, f"r = {np.corrcoef(ST.pc1_target, ST[col])[0, 1]:.2f}\nn = {len(ST)}", va="top", fontsize=5, linespacing=1.15)
    pg.src(lt, ST[["cohort", "target", "k", "pc1_target", col]].assign(cohort=ST.cohort.map(cl_)))

# ---------------------------------------------------------------- g  permuted outcome
pg.panel(120.5, 50.5, "g", "Permuted outcome")
ax = pg.ax(127.5, Y2, 21, H2)
XG = 13
hist(ax, GP.R, np.arange(-0.5, XG + 0.5, 1))
ax.axvline(GP.R.mean(), color=INK, lw=0.7, zorder=4)
ax.axvline(GP.E.mean(), color=EST, lw=0.7, ls=(0, (2.2, 1.6)), zorder=5)
ax.set_xlim(-0.6, XG)
ax.set_ylim(0, 110)
ax.set_xticks([0, 5, 10])
ax.set_yticks([0, 50, 100])
ax.set_xlabel("false calls per library")
ax.set_ylabel("permutations")
pg.txt(148.5, Y2 + 0.4, f"ratio {GN.ratio:.2f}\n({GN.ratio_lo:.2f}–{GN.ratio_hi:.2f})", ha="right", va="top", linespacing=1.15)
key(pg, 134.5, Y2 + 7.6, [(INK, f"mean {GP.R.mean():.2f}", "line")])
key(pg, 134.5, Y2 + 10.2, [(EST, f"estimate {GP.E.mean():.2f}", "line")])
beyond = sorted(GP.R[GP.R >= XG].astype(int))
pg.txt(148.5, Y2 + H2 - 3.0, f"› {', '.join(map(str, beyond))}", ha="right")
pg.src("g", GP[["perm", "R", "E"]].assign(realized_mean=GN.realized, estimated_mean=GN.estimated, ratio=GN.ratio, ratio_lo=GN.ratio_lo, ratio_hi=GN.ratio_hi))

# ---------------------------------------------------------------- h  ten decoy sets
pg.panel(152, 50.5, "h", "Ten decoy sets")
ax = pg.ax(160, Y2, 19.5, H2)
out = []
for i, (k_, lab) in enumerate((("pc1", "PC1"), ("abs_pc1", "abs. r"), ("mean", "mean z"))):
    d_ = B[k_][B[k_].stat == "fold"]
    point(ax, i + np.linspace(-0.26, 0.26, len(d_)), d_.fdr, EST, ms=2.3)
    assert abs(d_.fdr.min() - row(F[k_]).fdr_set_min) < 1e-6 and abs(d_.fdr.max() - row(F[k_]).fdr_set_max) < 1e-6
    pool = row(F[k_]).fdr
    ax.plot([i - 0.34, i + 0.34], [pool, pool], color=INK, lw=0.8, zorder=5, solid_capstyle="butt")
    ax.text(i, d_.fdr.max() + 2.2, f"{pool:.1f}", ha="center", va="bottom", fontsize=5)
    out += [(lab, int(r_.set), r_.fdr, pool) for r_ in d_.itertuples()]
ax.set_xticks(range(3), ["PC1", "abs. r", "mean z"], fontsize=5)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.55, 2.55)
ax.set_ylim(0, 55)
ax.set_yticks([0, 20, 40])
ax.set_ylabel("estimated FDR (%)")
pg.src("h", pd.DataFrame(out, columns=["score", "set", "fdr", "pooled"]))

# ---------------------------------------------------------------- i  forest; rows shared with j, k and l
Y3, RH3 = 103.0, 2.0
c_ = CL[CL.stat == "fold"]
tt = c_[c_.kind == "target"].set_index("target")
dcx = c_[(c_.kind == "decoy") & (c_["class"] == "coex")]
qd = dcx.groupby("target").med_delta.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
sorder = list(tt.med_delta.sort_values(ascending=False).index)
NS = len(sorder)
assert NS == 45
dev = set(EXC.target)
pg.panel(1, 89, "i", "Signature and 320 co-expressed decoys")
ax = pg.ax(20.5, Y3, 34, RH3 * NS)
for i, s in enumerate(sorder):
    q = qd.loc[s].to_numpy()
    ax.plot([q[0], q[4]], [i, i], color=tint(COEX, 0.45), lw=0.6, solid_capstyle="butt", zorder=2)
    ax.plot([q[1], q[3]], [i, i], color=COEX, lw=1.4, solid_capstyle="butt", zorder=3)
    ax.plot([q[2], q[2]], [i - 0.36, i + 0.36], color=INK, lw=0.5, solid_capstyle="butt", zorder=4)
    point(ax, tt.loc[s, "med_delta"], i, SIG, ms=2.8, filled=bool(tt.loc[s, "call"]), zorder=5)
    pg.txt(20, Y3 + RH3 * (i + 0.5), s, ha="right", **({"weight": "bold"} if s in dev else {}))
zero(ax, 0)
lo_ = min(qd[0.05].min(), tt.med_delta.min())
assert lo_ > -0.055 and tt.med_delta.max() < 0.108
ax.set_xlim(-0.055, 0.108)
ax.set_ylim(NS - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([-0.05, 0, 0.05, 0.1], ["−0.05", "0.00", "0.05", "0.10"])
ax.set_xlabel("median ΔC")
key(pg, 20.5, 93.6, [(SIG, "called", "o"), (SIG, "not called", "open-o")])
pg.line(20.5, 96.7, 22.3, 96.7, color=tint(COEX, 0.45), lw=0.6)
pg.txt(22.9, 96.7, "decoys 5–95%")
pg.line(20.5, 99.8, 22.3, 99.8, color=COEX, lw=1.4)
pg.txt(22.9, 99.8, "25–75%")
pg.line(32.5, 98.9, 32.5, 100.7, color=INK, lw=0.5)
pg.txt(33.5, 99.8, "median")
pg.src("i", qd.rename(columns=lambda q_: f"decoys_q{int(q_ * 100):02d}").assign(med_delta=tt.med_delta, called=tt.call, developed_here=[s in dev for s in qd.index]).loc[sorder].reset_index())

# ---------------------------------------------------------------- j  share of each signature's decoys called, same rows
pg.panel(54.6, 89, "j", "Decoys called")
ps = dcx.groupby("target").call.mean()
fp = row(F["pc1"])
assert abs(ps.sum() - fp.exp_false) < 1e-3 and abs(ps.mean() * 100 - fp.decoy_rate) < 1e-3
ax = pg.ax(57, Y3, 14.5, RH3 * NS)
point(ax, 100 * ps.loc[sorder].to_numpy(), range(NS), COEX, "s", ms=2.1)
zero(ax, 100 * ps.mean())
ax.set_ylim(NS - 0.5, -0.5)
ax.set_yticks([])
ax.set_xlim(-1, 34)
ax.set_xticks([0, 15, 30])
ax.set_xlabel("of 320 (%)")
pg.txt(59, 98.4, f"mean {100 * ps.mean():.1f}")
pg.txt(59, 100.7, f"sum {ps.sum():.2f}")
pg.src("j", pd.DataFrame({"signature": sorder, "share_decoys_called": ps.loc[sorder].to_numpy(), "n_decoys": 320}))

# ---------------------------------------------------------------- k  provenance, same rows
CW = 2.7
XK, XL = 75.6, 112.8
PROV = {"D": to_hex(PURPLE(0.95)), "C": to_hex(PURPLE(0.62)), "V": to_hex(PURPLE(0.32))}
check_palette("provenance", {**PROV, "blank": EMPTY})
pg.panel(XK - 0.5, 89, "k", "Provenance")
blk = {(r_.target, r_.cohort): r_.tier for r_ in BL.itertuples()}
out = []
for j, co in enumerate(COH):
    pg.txt(XK + CW * j + CW / 2, Y3 - 0.8, cl_(co), ha="center", va="bottom", rotation=90)
for i, s in enumerate(sorder):
    for j, co in enumerate(COH):
        code = PV.loc[s, co]
        code = code.strip() if isinstance(code, str) else ""
        x_, y_ = XK + CW * j, Y3 + RH3 * i
        tile(pg, x_, y_, CW, RH3, PROV.get(code, EMPTY), code or None)
        if (s, co) in blk:
            pg.c.plot(x_ + CW - 0.55, y_ + 0.42, "o", ms=1.0, mfc=INK if code != "C" else "white", mec="none", zorder=4)
        out.append((s, cl_(co), code, blk.get((s, co), "")))
pk = key(pg, XK, 195.2, [(PROV["D"], "developed (D)"), (PROV["C"], "shaped (C)")], gap=1.6)
pk = key(pg, XK, 198.1, [(PROV["V"], "tested only (V)")], gap=1.6)
pg.c.plot(pk + 0.9, 198.1, "o", ms=1.4, mfc=INK, mec="none", zorder=4)
pg.txt(pk + 2.4, 198.1, "borderline")
pg.src("k", pd.DataFrame(out, columns=["signature", "cohort", "code", "borderline_tier"]))

# ---------------------------------------------------------------- l  gain per signature and cohort, same rows
HL = 0.1
pg.panel(XL - 0.5, 89, "l", "Gain over the null")
HM = T.pivot_table(index="target", columns="cohort", values="dC").reindex(index=sorder, columns=COH)
for j, co in enumerate(COH):
    pg.txt(XL + CW * j + CW / 2, Y3 - 0.8, cl_(co), ha="center", va="bottom", rotation=90)
for i, s in enumerate(sorder):
    for j, co in enumerate(COH):
        tile(pg, XL + CW * j, Y3 + RH3 * i, CW, RH3, TEALOR(float(np.clip((HM.loc[s, co] + HL) / (2 * HL), 0, 1))))
pg.cbar(XL + 1, 194.0, 24, 1.3, TEALOR, -HL, HL, [-0.1, 0, 0.1], "ΔC", labels=["−0.1", "0", "≥ 0.1"])
pg.src("l", T[["target", "cohort", "fold_obs", "fold_med", "dC"]].assign(cohort=T.cohort.map(cl_)))

# ---------------------------------------------------------------- m  each signature's rank among its own decoys, the q behind pi0
sstat = lambda d_: np.sign(d_.med_delta) * -np.log10(np.maximum(d_.p, 1e-300))  # noqa: E731
qq = pd.Series({s: (1 + (sstat(g).to_numpy() >= sstat(tt.loc[s]) - 1e-12).sum()) / (1 + len(g)) for s, g in dcx.groupby("target")})
for r_ in PIL.itertuples():
    assert abs((qq > r_[1]).sum() / (NS * (1 - r_[1])) - r_.pi0) < 1e-3  # the pi0 in the table at each lambda
XR = 158.6
pg.panel(150.5, 89, "m", "Rank among decoys")
ax = pg.ax(XR, 94.5, 20, 20.5)
hist(ax, qq, np.linspace(0, 1, 11))
pi05 = float(PIL[PIL["lambda"] == 0.5].pi0.iloc[0])
zero(ax, pi05 * NS * 0.1, axis="y")
ax.set_xlim(0, 1)
ax.set_ylim(0, 37)
ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_yticks([0, 10, 20, 30])
ax.set_xlabel("q")
ax.set_ylabel("signatures")
ax.text(0.98, 36, f"n = {NS}\nπ0 = {pi05:.2f}", ha="right", va="top", fontsize=5, linespacing=1.15)
pg.src("m", pd.DataFrame({"signature": qq.index, "q": qq.to_numpy(), "called": tt.call.reindex(qq.index).to_numpy()}))

# ---------------------------------------------------------------- n  decoy co-expression (set 1)
pg.panel(150.5, 125, "n", "Decoy co-expression")
ax = pg.ax(XR, 130.5, 20, 20.5)
I1 = INFO[INFO.set == 1]
ax.plot(I1[I1["class"] == "coex"].r_target_med, I1[I1["class"] == "coex"].r_achieved_med, "s", ms=0.9, mfc=COEX, mec="none", ls="none", rasterized=True, zorder=3)
ax.plot(I1[I1["class"] == "rand"].r_target_med, I1[I1["class"] == "rand"].r_achieved_med, "D", ms=1.0, mfc=RAND, mec="none", ls="none", rasterized=True, zorder=4)
ref(ax, 0, 0.62)
ax.set_xlim(-0.03, 0.62)
ax.set_ylim(-0.03, 0.62)
ax.set_xticks([0, 0.3, 0.6], ["0", "0.3", "0.6"])
ax.set_yticks([0, 0.3, 0.6], ["0", "0.3", "0.6"])
ax.set_xlabel("signature r")
ax.set_ylabel("decoy r")
key(pg, XR + 0.5, 132.2, [(COEX, "co-expressed", "s")])
key(pg, XR + 0.5, 134.8, [(RAND, "random-gene", "D")])
pg.src("n", I1.drop(columns=["set"]))

# ---------------------------------------------------------------- o  call rate per signature under the permutations of g
pg.panel(150.5, 161, "o", "Call rate, permuted")
ax = pg.ax(XR, 166.5, 20, 20.5)
G_ = GT[GT.stat == "fold"].set_index("target")
for v_, f_ in ((True, True), (False, False)):
    s_ = tt.call.reindex(G_.index) == v_
    point(ax, 100 * G_.decoy_rate[s_], 100 * G_.rate[s_], SIG, ms=2.4, filled=f_, mew=None if f_ else 0.5)
ref(ax, 0, 6.5)
ax.set_xlim(0, 6.5)
ax.set_ylim(0, 6.5)
ax.set_xticks([0, 3, 6])
ax.set_yticks([0, 3, 6])
ax.set_xlabel("decoys (%)")
ax.set_ylabel("signature (%)")
ax.text(6.4, 0.25, f"r = {np.corrcoef(G_.rate, G_.decoy_rate)[0, 1]:.2f}, n = {len(G_)}", ha="right", va="bottom", fontsize=5)
pg.src("o", G_.reset_index()[["target", "rate", "decoy_rate"]].assign(called=tt.call.reindex(G_.index).to_numpy()))


# ---------------------------------------------------------------- p, q  Kaplan-Meier pair with numbers at risk
def km_panel(lt, x_, title, col, color, show_y):
    pg.panel(x_ - (9.5 if show_y else 2.5), 201.5, lt, title)
    w_, y_, h_ = 28.0, 207.0, 17.0
    ax_ = pg.ax(x_, y_, w_, h_)
    v_ = KM[col].to_numpy()
    hi_ = v_ > np.median(v_)
    light = tint(color, 0.5)
    for m_, cc in ((hi_, color), (~hi_, light)):
        xs, ys = km_curve(KM.time[m_], KM.event[m_])
        ax_.step(xs, ys, where="post", color=cc, lw=0.8, zorder=3)
        ct = KM.time[m_ & (KM.event == 0)].to_numpy()
        ax_.plot(ct, ys[np.clip(np.searchsorted(xs, ct, side="right") - 1, 0, len(ys) - 1)], "|", ms=1.5, mew=0.2, color=cc, ls="none", rasterized=True, zorder=4)
    p_ = logrank_p(KM.time, KM.event, hi_)
    ax_.set_xlim(0, 300)
    ax_.set_ylim(0, 1.04)
    ticks = [0, 100, 200, 300]
    ax_.set_xticks(ticks)
    ax_.set_yticks([0, 0.5, 1], ["0", "0.5", "1"] if show_y else [])
    ax_.set_xlabel("months")
    if show_y:
        ax_.set_ylabel("relapse-free")
    key(pg, x_ + 1.5, y_ + h_ - 6.4, [(color, "above median", "line")])
    key(pg, x_ + 1.5, y_ + h_ - 4.1, [(light, "below", "line")])
    # no superscript: two figures down to 0.001, below that p_small (0.0004, or < 0.0001)
    op, s_ = p_small(p_) if p_ < 0.001 else ("=", f"{p_:.2g}" if p_ < 0.01 else f"{p_:.2f}")
    pg.txt(x_ + 1.5, y_ + h_ - 1.6, f"log-rank P {op} {s_}")
    for jj, (mk, nm) in enumerate(((hi_, "above"), (~hi_, "below"))):
        yy = 232.6 + 2.6 * jj
        for tk in ticks:
            pg.txt(x_ + w_ * tk / 300, yy, str(int((KM.time[mk] >= tk).sum())), ha="center", color=color if jj == 0 else INK)
        if show_y:
            pg.txt(x_ - 2.2, yy, nm, ha="right")
    if show_y:
        pg.txt(x_ - 2.2, 229.9, "at risk", ha="right", weight="bold")
    pg.src(lt, KM.assign(above_median=hi_)[["time", "event", "above_median"]].sort_values("time").reset_index(drop=True).assign(log_rank_P=p_))
    return p_, [[int((KM.time[mk] >= tk).sum()) for tk in ticks] for mk in (hi_, ~hi_)]


kmv = [km_panel("p", 10.5, f"{EXS}, METABRIC", "target", SIG, True), km_panel("q", 43.0, "Median decoy, METABRIC", "decoy", COEX, False)]
assert kmv[0][1] == [[987, 441, 142, 4], [988, 567, 156, 6]] and kmv[1][1] == [[987, 473, 154, 8], [988, 535, 144, 2]]
assert f"{kmv[0][0]:.0e}" == "2e-14" and f"{kmv[1][0]:.0e}" == "1e-04"

# ---------------------------------------------------------------- r  the example across cohorts
pg.panel(75, 201.5, "r", f"{EXS} and its 32 co-expressed decoys")
Ud = UA[(UA.target == EXS) & (((UA.kind == "decoy") & (UA["class"] == "coex") & (UA.set == 1)) | (UA.kind == "target"))].copy()
Ud["dC"] = Ud.fold_obs - Ud.fold_med
oc = Ud[Ud.kind == "target"].set_index("cohort").dC.reindex(COH).sort_values(ascending=False).index
ax = pg.ax(82.5, 207, 57, 19.5)
for i, co in enumerate(oc):
    dv = Ud[(Ud.cohort == co) & (Ud.kind == "decoy")].dC.to_numpy()
    assert len(dv) == 32
    point(ax, i + swarm(dv, 0.0075, 0.13, 0.44), dv, COEX, "s", ms=1.5, mew=0.15, zorder=3)
    ax.plot([i - 0.4, i + 0.4], [np.median(dv)] * 2, color=INK, lw=0.7, zorder=5, solid_capstyle="butt")
    point(ax, i, Ud[(Ud.cohort == co) & (Ud.kind == "target")].dC.iloc[0], SIG, ms=3.6, zorder=6)
zero(ax, 0, axis="y")
ax.set_xticks(range(len(oc)), [cl_(c) for c in oc], rotation=90, fontsize=5)
ax.tick_params(axis="x", length=0, pad=1)
ax.set_xlim(-0.6, len(oc) - 0.4)
assert Ud.dC.min() > -0.13 and Ud.dC.max() < 0.14
ax.set_ylim(-0.13, 0.17)
ax.set_yticks([-0.1, 0, 0.1], ["−0.1", "0.0", "0.1"])
ax.set_ylabel("ΔC")
key(pg, 97, 208.2, [(SIG, EXS, "o"), (COEX, "32 decoys", "s"), (INK, "decoy median", "line")], sw=1.4, gap=1.5)
pg.src("r", Ud[["cohort", "unit", "kind", "dC"]].assign(cohort=Ud.cohort.map(cl_)))

# ---------------------------------------------------------------- s  above the null, by cohort, in the row order of a
pg.panel(143, 201.5, "s", "Above the null")
RHS = 1.85
YS = 207.0
ax = pg.ax(155, YS, 19.5, RHS * N1)
ab = UA[UA.fold_obs.notna()].assign(above=lambda d_: d_.fold_obs > d_.fold_p95, unit_class=lambda d_: d_["class"].where(d_.kind == "decoy", "signature"))
A = ab.groupby(["cohort", "unit_class"]).above.agg(["sum", "size"])
out = []
for i, co in enumerate(order):
    S_, C_, R_ = (A.loc[(co, k_)] for k_ in ("signature", "coex", "rand"))
    assert S_["size"] == NS
    ax.plot([0, 100], [i, i], color=EMPTY, lw=0.3, zorder=1)
    point(ax, 100 * R_["sum"] / R_["size"], i, RAND, "D", ms=2.0, zorder=3)
    point(ax, 100 * C_["sum"] / C_["size"], i, COEX, "s", ms=2.2, zorder=4)
    point(ax, 100 * S_["sum"] / S_["size"], i, SIG, "o", ms=2.6, zorder=5)
    pg.txt(154.4, YS + RHS * (i + 0.5), cl_(co), ha="right")
    pg.txt(179.6, YS + RHS * (i + 0.5), f"{int(S_['sum'])}", ha="right")
    out.append((cl_(co), int(S_["sum"]), int(S_["size"]), 100 * S_["sum"] / S_["size"], int(C_["sum"]), int(C_["size"]), 100 * C_["sum"] / C_["size"], int(R_["sum"]), int(R_["size"]), 100 * R_["sum"] / R_["size"]))
zero(ax, 5)
ax.set_ylim(N1 - 0.5, -0.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
XA = 90
ax.set_xlim(-2, XA)
ax.set_xticks([0, 40, 80])
ax.set_xlabel("above null 95th\npercentile (%)", linespacing=1.0, labelpad=1)
pg.txt(179.6, YS - 1.4, "of 45", ha="right", weight="bold")
key(pg, 143.2, YS - 1.4, [(SIG, "signatures", "o"), (COEX, "co-expr.", "s"), (RAND, "random", "D")], sw=1.4, gap=1.0)
SD = pd.DataFrame(out, columns=["cohort", "signatures_above", "signatures", "signatures_pct", "coex_above", "coex", "coex_pct", "rand_above", "rand", "rand_pct"])
assert SD[["signatures_pct", "coex_pct", "rand_pct"]].max().max() < XA
pg.src("s", SD)

pg.save("Fig7")
print("c:", [(r[0], int(r[1].called), round(float(r[1].fdr), 2), round(float(r[1].fdr_lo), 2), round(float(r[1].fdr_hi), 2), round(float(r[1].fdr_pi0), 2)) for r in rows])
print(SD.round(1).to_string())
print("sum p_s", round(ps.sum(), 3), "| q>0.5", int((qq > 0.5).sum()), "| dC>0.1:", int((T.dC > HL).sum()), "max", round(T.dC.max(), 3), "| KM P", [f"{k_[0]:.1e}" for k_ in kmv])
print("perms beyond axis", beyond, "| patients", int(mo.n_samples.sum()), "events", int(mo.n_events.sum()), "| flips", [s for s in a_.index if a_.call[s] != b_.call[s]])
