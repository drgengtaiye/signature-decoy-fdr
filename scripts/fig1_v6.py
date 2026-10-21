"""Figure 1 (v6) - one null for both settings, and what it shows.

a  the 52 signatures by size and scoring engine, with the sizes of the matched null
b  why the null has to be matched: its 95th percentile moves with size and engine
c  the matched null against a fixed 48-gene null, signature by signature (all 52: the fixed null does not
   depend on the signature, so each cohort's fixed-null median and 95th percentile apply to every signature)
d  response effect against prognostic effect, in null s.d., with family density regions and medians (hero)
e  the shift between settings by family (raincloud)
f  each signature's two effects, paired
g  the shift against gene-set size, with the size-adjusted family fit
h  the family contrast with each checkpoint cohort or cancer type removed
Statistics are in the legend; each is printed here so the legend can be checked.
"""
import numpy as np, pandas as pd, warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr, gaussian_kde
import statsmodels.formula.api as smf
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D
R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med; R["clear"] = R.AUC > R.rand_p95; R["z"] = (R.AUC - R.rand_med) / R.rand_sd
O = pd.read_csv(f"{D}/icb/native_response.csv"); O = O[O.cohort != "ICB_VanDenEnde"]
FIX = O.groupby("cohort")[["rand_med", "rand_p95"]].first()          # one fixed null per cohort
assert (O.groupby("cohort")[["rand_med", "rand_p95"]].nunique() == 1).all().all()
PP = pd.read_csv(f"{D}/icb/native_prognosis_matched.csv"); PP = PP[PP.endpoint == "OS"].copy()
PP["z"] = (PP.C - PP.rand_med) / PP.rand_sd
X = pd.read_csv(f"{D}/icb/v3_dissociation.csv")
X = X.merge(pd.read_csv(f"{D}/icb/v4_standardised_effects.csv")[["signature", "z_resp", "z_prog", "z_shift"]], on="signature")
V5 = pd.read_csv(f"{D}/icb/v5_dependence_aware.csv")
PERM = V5[V5.statistic == "family contrast: rank-biserial r"].iloc[0]
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]; ORD3 = ["Other", "Confounder", "Inflamed/IFN"]
famd = X.set_index("signature").family

C = K.Canvas(202)

# ============================================================ a  the library
SIZ = (R.drop_duplicates("signature")[["signature", "n_genes", "engine"]].assign(family=lambda d: d.signature.map(famd)))
LADDER = sorted(R.null_size.unique()); MEDSZ = float(SIZ.n_genes.median())
ax = C.ax(27, 9, 118, 24)
rng = np.random.default_rng(7)
for i, (e, name) in enumerate([("gsva", "GSVA"), ("wmean", "weighted mean"), ("ssgsea", "ssGSEA")]):
    g = SIZ[SIZ.engine == e]
    ax.scatter(g.n_genes, i + rng.uniform(-.24, .24, len(g)), s=13, c=[fs.FAM[f] for f in g.family],
               lw=.35, ec="white", zorder=4, clip_on=False)
    ax.text(2.55, i, f"{name} ({len(g)})", ha="right", va="center", fontsize=K.T_LAB)
ax.axvline(MEDSZ, color=fs.INK, lw=.6, ls=(0, (3, 2)), zorder=2)
ax.text(MEDSZ * 1.07, -.72, f"median {MEDSZ:g} genes", fontsize=K.T_SMALL, va="center")
for r_ in LADDER: ax.plot([r_, r_], [2.62, 2.86], color=fs.SLATE, lw=.8, clip_on=False)
ax.text(2.55, 2.74, "null sizes", ha="right", va="center", fontsize=K.T_SMALL, color=fs.SLATE)
ax.set_xscale("log"); ax.set_xlim(2.7, 560); ax.set_ylim(2.95, -1.0); ax.minorticks_off()
ax.set_xticks([3, 10, 30, 100, 300]); ax.set_xticklabels(["3", "10", "30", "100", "300"])
ax.set_yticks([]); K.frame(ax, left=False); ax.set_xlabel("genes in signature")
K.key(ax, [(fs.FAM[f], fs.FAMLAB[f]) for f in FAMORD], x=1.035, y=.86, dy=.19, s=16)
C.letter(1.5, 3, "a")

# ============================================================ b  the null moves
ax = C.ax(14, 52, 42, 34)
ECOL = {"gsva": fs.SLATE, "ssgsea": fs.TEAL, "wmean": fs.AMBER}
Q = R.groupby(["engine", "null_size"]).agg(p95=("rand_p95", "median"), n=("rand_p95", "size")).reset_index()
for e in ["gsva", "ssgsea", "wmean"]:
    g = Q[(Q.engine == e) & (Q.n >= 10)].sort_values("null_size")
    ax.plot(g.null_size, g.p95, color=ECOL[e], lw=1.1, zorder=3)
    ax.scatter(g.null_size, g.p95, s=np.clip(g.n / 4, 5, 22), c=ECOL[e], lw=.35, ec="white", zorder=4)
xk = .02
for e, nm in [("gsva", "GSVA"), ("ssgsea", "ssGSEA"), ("wmean", "weighted mean")]:
    ax.text(xk, .985, nm, transform=ax.transAxes, color=ECOL[e], fontsize=K.T_SMALL, va="top")
    xk += {"gsva": .21, "ssgsea": .27, "wmean": 0}[e]
old95 = R.cohort.map(FIX.rand_p95).median()          # over the same 623 tests
ax.axhline(old95, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
ax.text(62, old95 + .003, "fixed 48-gene null", va="bottom", fontsize=K.T_SMALL)
ax.set_xscale("log"); ax.set_xlim(3.2, 560); ax.minorticks_off()
ax.set_xticks([4, 10, 40, 100, 400]); ax.set_xticklabels(["4", "10", "40", "100", "400"])
ax.set_ylim(.59, .75); ax.set_yticks([.6, .65, .7, .75])
ax.set_xlabel("random gene-set size"); ax.set_ylabel("null 95th percentile (AUC)"); K.frame(ax)
C.letter(1.5, 46, "b")

# ============================================================ c  fixed against matched
ax = C.ax(18, 104, 36, 36)
m = R.assign(delta_old=R.AUC - R.cohort.map(FIX.rand_med), rand_p95_old=R.cohort.map(FIX.rand_p95))
s = m.groupby("signature")[["delta", "delta_old"]].median(); s["family"] = s.index.map(famd)
for f in ORD3:
    g = s[s.family == f]; ax.scatter(g.delta_old, g.delta, s=12, c=fs.FAM[f], lw=.3, ec="white", zorder=3)
lo_, hi_ = -.025, .125
ax.plot([lo_, hi_], [lo_, hi_], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.axhline(0, color=fs.LIGHT, lw=.6); ax.axvline(0, color=fs.LIGHT, lw=.6)
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_xticks([0, .05, .1]); ax.set_yticks([0, .05, .1])
ax.set_xlabel("median ΔAUC, fixed null"); ax.set_ylabel("median ΔAUC, matched null"); K.frame(ax)
rr = spearmanr(s.delta_old, s.delta)[0]
C.letter(1.5, 98, "c")

# ============================================================ d  the dissociation (hero)
xl, yl = (-.85, 1.65), (-.55, 2.65)
ax = C.ax(84, 62, 66, 72); mt = C.ax(84, 52.5, 66, 8.5); mr = C.ax(151.5, 62, 8.5, 72)
for f in ORD3:
    g = X[X.family == f]
    ax.scatter(g.z_prog, g.z_resp, s=22 if f != "Other" else 15, c=fs.FAM[f], lw=.45, ec="white",
               zorder=4 if f != "Other" else 3, alpha=1 if f != "Other" else .75)
    gx = np.linspace(*xl, 200); gy = np.linspace(*yl, 200)
    kx = gaussian_kde(g.z_prog)(gx); ky = gaussian_kde(g.z_resp)(gy); a_ = .40 if f != "Other" else .20
    mt.fill_between(gx, kx / kx.max(), color=fs.FAM[f], alpha=a_, lw=0); mt.plot(gx, kx / kx.max(), color=fs.FAM[f], lw=.6)
    mr.fill_betweenx(gy, ky / ky.max(), color=fs.FAM[f], alpha=a_, lw=0); mr.plot(ky / ky.max(), gy, color=fs.FAM[f], lw=.6)
mt.set_xlim(xl); mt.set_ylim(0, 1.1); mr.set_ylim(yl); mr.set_xlim(0, 1.1); fs.clean(mt); fs.clean(mr)
# joint density: 50% and 80% highest-density regions of the two contrasted families, and each family's median with a
# bootstrap 95% interval (2,000 resamples of signatures)
GX, GY = np.meshgrid(np.linspace(*xl, 160), np.linspace(*yl, 160)); HDR = {}
for f in ["Inflamed/IFN", "Confounder"]:
    g = X[X.family == f]; kd = gaussian_kde(np.vstack([g.z_prog, g.z_resp]))
    Zd = kd(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape); zs = np.sort(Zd.ravel())[::-1]; cm = np.cumsum(zs) / zs.sum()
    lv = [zs[np.searchsorted(cm, q)] for q in (.8, .5)]; HDR[f] = lv
    ax.contourf(GX, GY, Zd, levels=[lv[0], lv[1], Zd.max() * 1.01], colors=[fs.FAM[f]], alpha=.10, zorder=1)
    ax.contour(GX, GY, Zd, levels=lv, colors=[fs.FAM[f]], linewidths=[.4, .7], alpha=.8, zorder=2)
    bs = np.array([np.median(g[["z_prog", "z_resp"]].to_numpy()[rng.integers(0, len(g), len(g))], axis=0) for _ in range(2000)])
    mx, my = g.z_prog.median(), g.z_resp.median(); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    HDR[f] = (mx, my, lo.round(2).tolist(), hi.round(2).tolist())
    ax.errorbar(mx, my, xerr=[[mx - lo[0]], [hi[0] - mx]], yerr=[[my - lo[1]], [hi[1] - my]], fmt="D", ms=4.2, mfc="white",
                mec=fs.FAM[f], mew=1.0, ecolor=fs.FAM[f], elinewidth=.9, capsize=0, zorder=7)
ax.plot([-1, 2.7], [-1, 2.7], color=fs.INK, lw=.5, ls=":", zorder=2)
ax.text(1.60, 1.72, "y = x", ha="right", va="bottom", fontsize=K.T_SMALL, color=fs.SLATE)
ax.axhline(0, color=fs.LIGHT, lw=.7, zorder=1); ax.axvline(0, color=fs.LIGHT, lw=.7, zorder=1)
ax.set_xlim(xl); ax.set_ylim(yl); ax.set_xticks([-.5, 0, .5, 1, 1.5]); ax.set_yticks([0, 1, 2])
ax.set_xlabel("prognostic effect (null s.d.)"); ax.set_ylabel("response effect (null s.d.)"); K.frame(ax)
xs_ = X.set_index("signature")
LBL = {"Chemokine_Messina": (-.80, 2.50, "left"), "APM_Thompson": (.25, 2.50, "left"),
       "CD39-CD8Tcell_Chow": (1.62, 2.18, "right"), "proliferation_Pabla": (1.62, .98, "right"),
       "CIN70_Carter": (1.62, -.32, "right")}
TWO = {"proliferation_Pabla": "proliferation\nPabla"}
for sg, (tx, ty, ha) in LBL.items():
    r = xs_.loc[sg]
    ax.annotate(TWO.get(sg, K.lab(sg)), (r.z_prog, r.z_resp), xytext=(tx, ty), textcoords="data",
                fontsize=K.T_SMALL, ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", lw=.4, color=fs.SLATE, shrinkA=1, shrinkB=2.2), zorder=6)
C.letter(70, 46, "d")

# ============================================================ row 4
RY, RH = 160, 34
# e  shift by family
ax = C.ax(14, RY, 28, RH)
for i, f in enumerate(FAMORD):   # raincloud: half-density, interquartile box, points
    v = X.loc[X.family == f, "z_shift"].to_numpy(); gy = np.linspace(v.min() - .35, v.max() + .35, 200)
    kd = gaussian_kde(v)(gy); kd = kd / kd.max() * .34
    ax.fill_betweenx(gy, i - .06 - kd, i - .06, color=fs.FAM[f], alpha=.35, lw=0); ax.plot(i - .06 - kd, gy, color=fs.FAM[f], lw=.5)
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    ax.add_patch(plt.Rectangle((i - .03, q1), .06, q3 - q1, fc="white", ec=fs.INK, lw=.5, zorder=5))
    ax.plot([i - .03, i + .03], [q2, q2], color=fs.INK, lw=1.0, zorder=6)
    ax.scatter(i + .16 + K.swarm(v, .12, 16) * .9, v, s=7, c=fs.FAM[f], lw=.25, ec="white", zorder=3)
ax.axhline(0, color=fs.LIGHT, lw=.7, zorder=1)
ax.plot([0, 0, 2, 2], [2.55, 2.68, 2.68, 2.55], color=fs.INK, lw=.55)
ax.text(1, 2.76, f"$P$ = {PERM.P:.3f}", ha="center", va="bottom", fontsize=K.T_NOTE)
ax.set_xticks(range(3)); ax.set_xticklabels(["inflamed/\nIFN", "other", "confounder\ncontrol"], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-.5, 2.45); ax.set_ylim(-1.35, 3.1)
ax.set_ylabel("shift (null s.d.)"); K.frame(ax)
C.letter(1.5, RY - 6, "e")
# f  paired effects
ax = C.ax(58, RY, 28, RH)
for f in ORD3:
    g = X[X.family == f]; a_ = .85 if f != "Other" else .30
    for _, r in g.iterrows():
        ax.plot([0, 1], [r.z_prog, r.z_resp], color=fs.FAM[f], lw=.45, alpha=a_, zorder=3 if f != "Other" else 2)
    ax.scatter(np.zeros(len(g)), g.z_prog, s=5, c=fs.FAM[f], lw=0, zorder=4, alpha=a_)
    ax.scatter(np.ones(len(g)), g.z_resp, s=5, c=fs.FAM[f], lw=0, zorder=4, alpha=a_)
for f in ["Inflamed/IFN", "Confounder"]:
    g = X[X.family == f]
    ax.plot([1.08, 1.25], [g.z_resp.median()] * 2, color=fs.FAM[f], lw=1.8, solid_capstyle="butt")
    ax.plot([-.25, -.08], [g.z_prog.median()] * 2, color=fs.FAM[f], lw=1.8, solid_capstyle="butt")
ax.axhline(0, color=fs.LIGHT, lw=.6, zorder=1)
ax.set_xticks([0, 1]); ax.set_xticklabels(["prognosis", "response"], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_xlim(-.32, 1.32); ax.set_ylim(-.8, 2.6); ax.set_ylabel("effect (null s.d.)"); K.frame(ax)
C.letter(45, RY - 6, "f")
# g  shift and gene-set size
ax = C.ax(102, RY, 30, RH)
for f in ORD3:
    g = X[X.family == f]; ax.scatter(g.n_genes, g.z_shift, s=11, c=fs.FAM[f], lw=.3, ec="white", zorder=3)
ax.set_xscale("log"); ax.axhline(0, color=fs.LIGHT, lw=.6); ax.minorticks_off()
ax.set_xticks([3, 10, 30, 100, 300]); ax.set_xticklabels(["3", "10", "30", "100", "300"])
ax.set_xlabel("genes in signature"); ax.set_ylabel("shift (null s.d.)"); ax.set_ylim(-1.3, 2.6); K.frame(ax)
rs, ps = spearmanr(X.n_genes, X.z_shift)
mfit = smf.ols("z_shift ~ lg + infl", data=X.assign(lg=np.log10(X.n_genes), infl=(X.family == "Inflamed/IFN").astype(int))).fit()
gl = np.linspace(np.log10(3), np.log10(430), 60)
for inf_, col in [(1, fs.FAM["Inflamed/IFN"]), (0, fs.SLATE)]:     # size-adjusted fit, 95% confidence band
    pr = mfit.get_prediction(pd.DataFrame({"lg": gl, "infl": inf_})).summary_frame(alpha=.05)
    ax.fill_between(10 ** gl, pr.mean_ci_lower, pr.mean_ci_upper, color=col, alpha=.13, lw=0, zorder=1)
    ax.plot(10 ** gl, pr["mean"], color=col, lw=.9, zorder=2)
ax.text(.98, .97, "inflamed/IFN", transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL, color=fs.FAM["Inflamed/IFN"])
ax.text(.98, .88, "all others", transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL, color=fs.SLATE)
C.letter(89, RY - 6, "g")
# h  leave-one-out
ax = C.ax(147, RY, 31, RH)
def contrast(Rsub, Psub):
    dr = Rsub.groupby("signature").z.median(); dp = Psub.groupby("signature").z.median()
    sh = (dr - dp).reindex(X.signature); fam = X.set_index("signature").family
    a_, b_ = sh[fam == "Inflamed/IFN"], sh[fam == "Confounder"]
    return a_.median() - b_.median(), mannwhitneyu(a_, b_).pvalue
full, _ = contrast(R, PP)
rows = [("ICB", c, *contrast(R[R.cohort != c], PP)) for c in sorted(R.cohort.unique())]
rows += [("TCGA", c, *contrast(R, PP[PP.cancer != c])) for c in sorted(PP.cancer.unique())]
L = pd.DataFrame(rows, columns=["set", "left_out", "contrast", "p"])
for i, st in enumerate(["ICB", "TCGA"]):
    g = L[L.set == st]
    col = [fs.TUMCOL[fs.ICB_TUM[c]] for c in g.left_out] if st == "ICB" else fs.SLATE
    ax.scatter(np.full(len(g), i) + K.swarm(g.contrast.values, .26, 10), g.contrast, s=10, c=col, lw=.3, ec="white", zorder=3)
ax.axhline(full, color=fs.BLUE, lw=.9, ls=(0, (3, 2)))
ax.text(.5, full + .03, "all data", ha="center", va="bottom", fontsize=K.T_SMALL, color=fs.BLUE)
ax.axhline(0, color=fs.LIGHT, lw=.6)
ax.set_xticks([0, 1]); ax.set_xticklabels(["cohort\nremoved", "cancer type\nremoved"], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-.6, 1.6); ax.set_ylim(0, full * 1.45)
ax.set_ylabel("inflamed − confounder\nshift (null s.d.)"); K.frame(ax)
C.letter(133, RY - 6, "h")

C.save("Figure1")
vi = X.loc[X.family == "Inflamed/IFN", "z_shift"]; vc = X.loc[X.family == "Confounder", "z_shift"]
u = mannwhitneyu(vi, vc); rb = 2 * u.statistic / (len(vi) * len(vc)) - 1
rho, prho = spearmanr(X.z_prog, X.z_resp); pw = wilcoxon(X.z_prog, X.z_resp).pvalue
print(f"Figure1 v6 | b fixed 95th {old95:.3f} | h leave-one-out {L.contrast.min():.2f}-{L.contrast.max():.2f} (all {full:.2f}) max P {L.p.max():.1e}")
print(f"Figure1 v6 | c rho {rr:.2f} n {len(s)} fixed {100*(m.AUC > m.rand_p95_old).mean():.1f}% matched {100*m.clear.mean():.1f}%"
      f" | d rho {rho:+.2f} P {prho:.2f} paired {pw:.2g} | e rb {rb:.2f} perm {PERM.P:.3f}"
      f" | g rho {rs:+.2f} P {ps:.2f} size-adj {mfit.pvalues['infl']:.1e} | h full {full:.2f} min {L.contrast.min():.2f} maxP {L.p.max():.1e}")
print('d medians (prog, resp) [lo] [hi]:', {k: (round(v[0], 2), round(v[1], 2), v[2], v[3]) for k, v in HDR.items()}, '| g infl coef', round(mfit.params['infl'], 2), mfit.conf_int().loc['infl'].round(2).tolist())
