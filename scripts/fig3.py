"""Figure 3 (GB v10) - immune signatures predict response to checkpoint blockade (restyled from fig3_v6.py; same data).

a  every signature in every checkpoint cohort against its matched null (hero)
b  how often each family clears the null
c  the number of significant signatures against a permutation null that keeps their correlation
d  absolute discrimination: published signatures against their matched random sets
e  per cohort, how often inflamed/IFN and confounder signatures clear the null
f  median gain by tumor type within each family (descriptive)
g  each signature's median gain against its signed-rank P over cohorts (volcano)
h  every test's gain by family (raincloud)
i  the share of signatures clearing the null against cohort size
j  each signature's median gain against its number of genes
"""
from scipy.stats import binom
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.stats import wilcoxon
from statsmodels.stats.proportion import proportion_confint
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import cns10 as K
import figstyle as fs
K.apply()
D = K.D

R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv")
R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med
R["clear"] = R.AUC > R.rand_p95
X = pd.read_csv(f"{D}/icb/v3_dissociation.csv").set_index("signature")
R["family"] = R.signature.map(X.family)
famd = X.family
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
TUMS = ["melanoma", "kidney", "lung", "bladder", "gastric"]
sigp = R.groupby("signature").delta.apply(lambda v: wilcoxon(v).pvalue)
med = R.groupby("signature").delta.median()
nsig = int((sigp < .05).sum())

C = K.Canvas(226)

# ============================================================ a  the matrix
M = R.pivot_table(index="signature", columns="cohort", values="delta")
rows = []
for f in FAMORD:
    rows += sorted([s for s in M.index if famd[s] == f], key=lambda s: -med[s])
M = M.loc[rows]
ncoh = R.groupby("cohort").n.first()
cols = sorted(M.columns, key=lambda c: (TUMS.index(fs.ICB_TUM[c]), -ncoh[c]))
M = M[cols]
nr, nc = M.shape
MX, MY, MW, MH = 38, 30, 63, 108
LIM = .15
ax = C.ax(MX, MY, MW, MH)
ax.imshow(M.values, aspect="auto", cmap=K.DIV, norm=Normalize(-LIM, LIM), interpolation="nearest")
for i in range(nr + 1): ax.axhline(i - .5, color="white", lw=.3)
for j in range(nc + 1): ax.axvline(j - .5, color="white", lw=.3)
CL = R.pivot_table(index="signature", columns="cohort", values="clear").loc[rows][cols]
yy_, xx_ = np.where(CL.values == 1)
ax.scatter(xx_, yy_, s=1.6, color=K.INK, lw=0, zorder=4)
ax.set_yticks(range(nr)); ax.set_yticklabels([K.lab(s) for s in M.index], fontsize=K.T_SMALL)
for t, s in zip(ax.get_yticklabels(), M.index):
    t.set_color(K.INK if famd[s] != "Other" else K.SLATE)
ax.set_xticks(range(nc))
ax.set_xticklabels([f"{c.replace('ICB_', '').replace('_', ' ')} ({ncoh[c]})" for c in cols],
                   rotation=90, fontsize=K.T_SMALL)
ax.tick_params(length=0, pad=1.5)
for sp in ax.spines.values(): sp.set_visible(False)
# family strip with names
st = C.ax(4.5, MY, 1.6, MH); st.set_xlim(0, 1); st.set_ylim(nr - .5, -.5); st.axis("off")
start = 0
for f in FAMORD:
    n = int(sum(famd[s] == f for s in M.index))
    st.add_patch(Rectangle((0, start - .5), 1, n, fc=K.FAM[f], ec="white", lw=.5))
    st.text(-.9, start + n / 2 - .5, K.FAMLAB[f].replace(" control", ""), rotation=90,
            ha="center", va="center", fontsize=K.T_SMALL, color=K.FAM[f], fontweight="bold")
    start += n
# tumor groups, named above the columns
tt = C.ax(MX, MY - 3.0, MW, 2.4); tt.set_xlim(-.5, nc - .5); tt.set_ylim(0, 1); tt.axis("off")
start = 0
for tname in TUMS:
    n_ = sum(fs.ICB_TUM[c] == tname for c in cols)
    if not n_: continue
    tt.plot([start - .4, start + n_ - .6], [0.15, 0.15], color=K.INK2, lw=0.6, solid_capstyle="butt")
    tt.text(start + (n_ - 1) / 2, 0.45, tname if n_ > 1 else tname[:3] + ".", ha="center", va="bottom", fontsize=K.T_SMALL, color=K.INK)
    start += n_
# top: share of signatures clearing the null, per cohort
top = C.ax(MX, MY - 13.5, MW, 8.2)
pc = [100 * R[R.cohort == c].clear.mean() for c in cols]
TUMC = {"melanoma": K.BLUE, "kidney": K.ORANGE, "lung": K.AQUA, "bladder": K.SLATE, "gastric": "#b8862b"}
top.bar(range(nc), pc, width=.72, color=[TUMC[fs.ICB_TUM[c]] for c in cols], lw=0)
top.set_xlim(-.5, nc - .5); top.set_ylim(0, 60); top.set_yticks([0, 50]); top.set_xticks([])
top.set_ylabel("% clear", fontsize=K.T_SMALL, labelpad=1.5)
top.tick_params(axis="y", labelsize=K.T_SMALL)
K.frame(top)
# right: median gain per signature
rb = C.ax(MX + MW + 2, MY, 17, MH)
rb.barh(range(nr), med[M.index].values, height=.74, color=[K.FAM[famd[s]] for s in M.index], lw=0)
for i, s in enumerate(M.index):
    if sigp[s] < .05:
        rb.text(max(med[s], 0) + .004, i, "*", fontsize=K.T_NOTE, va="center", ha="left")
rb.set_ylim(nr - .5, -.5); rb.set_yticks([]); rb.axvline(0, color=K.INK, lw=.45)
rb.set_xlim(-.015, .13); rb.set_xticks([0, .1]); rb.set_xticklabels(["0", "0.1"])
rb.xaxis.tick_top(); rb.xaxis.set_label_position("top")
rb.set_xlabel("median ΔAUC", fontsize=K.T_SMALL, labelpad=2)
rb.tick_params(axis="x", labelsize=K.T_SMALL)
rb.spines["left"].set_visible(False); rb.spines["bottom"].set_visible(False)
rb.spines["top"].set_visible(True); rb.spines["right"].set_visible(False)
# colour bar and tumor key, in the corner above the row labels
cax = C.ax(9, 9.5, 22, 1.8)
cb = C.fig.colorbar(plt.cm.ScalarMappable(Normalize(-LIM, LIM), K.DIV), cax=cax,
                    orientation="horizontal", ticks=[-.15, 0, .15])
cb.outline.set_linewidth(.3); cax.tick_params(labelsize=K.T_SMALL, length=1.5, width=.4, pad=1)
cax.set_title("ΔAUC vs matched null", fontsize=K.T_SMALL, pad=2)
C.text(9, 15.3, "•  clears null 95th pct", size=K.T_SMALL, color=K.INK)
C.letter(1.5, 3, "a")

# ============================================================ b  clearance by family
ax = C.ax(141, 14, 36, 34)
PSIG = R.groupby(["family", "signature"]).clear.mean().mul(100).reset_index()
for i, f in enumerate(FAMORD):
    g = R[R.family == f]; k, n = int(g.clear.sum()), len(g)
    lo, hi = proportion_confint(k, n, method="wilson")
    v = PSIG.loc[PSIG.family == f, "clear"].to_numpy()          # each signature's own clearance rate
    ax.add_patch(Rectangle((i - .38, 100 * lo), .76, 100 * (hi - lo), fc=K.tint(K.FAM[f], .82), ec="none", zorder=1))
    ax.plot([i - .38, i + .38], [100 * k / n] * 2, color=K.FAM[f], lw=1.4, zorder=2, solid_capstyle="butt")
    ax.scatter(i + K.swarm(v, .32, 14), v, s=7, c=K.FAM[f], lw=.3, ec="white", zorder=3)
K.ref(ax, 5)
ax.text(2.52, 4.2, "chance (5%)", ha="right", va="top", fontsize=K.T_SMALL, color=K.INK2)
ax.set_xticks(range(3)); ax.set_xticklabels(["inflamed/\nIFN", "other", "confounder\ncontrol"],
                                             fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0)
ax.set_ylim(-2, 56); ax.set_xlim(-.55, 2.55); ax.set_yticks([0, 10, 20, 30, 40, 50])
ax.set_ylabel("tests clearing the null (%)")
K.frame(ax)
C.letter(128, 8, "b")

# ============================================================ c  count vs permutation null
PERM = pd.read_csv(f"{D}/icb/v4_permutation_inference.csv")
ph = PERM[(PERM.analysis == "12 cohorts (headline)") &
          (PERM.statistic == "signatures with P<0.05")].iloc[0]
ax = C.ax(141, 67, 36, 22)
DR = pd.read_csv(f"{D}/icb/v4_permutation_null_draws.csv.gz")
dv = DR[(DR.analysis == "12 cohorts (headline)") & (DR.statistic == "signatures with P<0.05")].value.to_numpy()
kk = np.arange(0, 31); fr = np.array([(dv == k_).mean() for k_ in kk]) * 100
frp = np.where(fr > 0, fr, 1e-3)
xs_ = np.repeat(np.r_[kk - .5, kk[-1] + .5], 2)[1:-1]; ys_ = np.repeat(frp, 2)
ax.fill_between(xs_, 1e-3, ys_, where=xs_ <= nsig + .5, color=K.tint(K.SLATE, .45), lw=0, zorder=2, step=None)
ax.fill_between(xs_, 1e-3, ys_, where=xs_ >= nsig - .5, color=K.tint(K.BLUE, .2), lw=0, zorder=2)
ax.plot(xs_, ys_, color=K.SLATE, lw=.5, zorder=2.5)
pmf = binom.pmf(kk, 52, .05) * 100                                   # what 52 independent tests would give
ax.plot(kk[pmf > 1e-2], pmf[pmf > 1e-2], color=K.ORANGE, lw=.8, zorder=3)
ax.axvline(nsig, color=K.BLUE, lw=1.1, zorder=4)
ax.text(nsig + .6, 30, f"observed {nsig}\n$P$ = {ph.perm_P:.3f}", ha="left", va="top", fontsize=K.T_SMALL, color=K.INK)
ax.scatter([.02], [1.1], s=10, c=K.SLATE, marker="s", lw=0, transform=ax.transAxes, clip_on=False)
ax.text(.06, 1.1, "permutation null", transform=ax.transAxes, fontsize=K.T_SMALL, va="center")
ax.plot([.52, .6], [1.1, 1.1], color=K.ORANGE, lw=.8, transform=ax.transAxes, clip_on=False)
ax.text(.63, 1.1, "binomial", transform=ax.transAxes, fontsize=K.T_SMALL, va="center")
ax.set_yscale("log"); ax.set_ylim(.01, 60); ax.set_yticks([.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"]); ax.minorticks_off()
ax.set_xlim(-1, 30); ax.set_xticks([0, 10, 20, 30])
ax.set_xlabel("signatures with $P$ < 0.05 (of 52)"); ax.set_ylabel("permutations (%)")
K.frame(ax)
C.letter(128, 59, "c")

# ============================================================ d  published vs random
ax = C.ax(141, 106, 36, 32)
bins = np.arange(.30, .96, .025)
h1, _ = np.histogram(R.AUC, bins); h0, _ = np.histogram(R.rand_med, bins)
xs_ = np.repeat(bins, 2)[1:-1]
ax.fill_between(xs_, 0, np.repeat(h1, 2), color=K.tint(K.BLUE, .35), lw=0); ax.plot(xs_, np.repeat(h1, 2), color=K.BLUE, lw=.6)
ax.fill_between(xs_, 0, -np.repeat(h0, 2), color=K.tint(K.SLATE, .6), lw=0); ax.plot(xs_, -np.repeat(h0, 2), color=K.SLATE, lw=.6)
ax.axhline(0, color=K.INK, lw=.45)
ma, mr = R.AUC.median(), R.rand_med.median()
ax.plot([ma] * 2, [0, h1.max() * 1.04], color=K.INK, lw=.7); ax.plot([mr] * 2, [0, -h0.max() * 1.02], color=K.INK, lw=.7)
ax.text(.95, h1.max() * .95, f"published\nmedian {ma:.3f}", ha="right", va="top", fontsize=K.T_SMALL, color=K.BLUE, linespacing=1.0)
ax.text(.95, -h0.max() * .55, f"matched random\nmedian {mr:.3f}", ha="right", va="top", fontsize=K.T_SMALL, color=K.SLATE, linespacing=1.0)
ax.set_xlim(.33, .96); ax.set_ylim(-h0.max() * 1.08, h1.max() * 1.12)
yt = [-200, -100, 0, 100]
ax.set_yticks(yt); ax.set_yticklabels([str(abs(v)) for v in yt])
ax.set_xticks([.4, .6, .8])
ax.set_xlabel("AUC for response"); ax.set_ylabel("signature × cohort tests")
K.frame(ax)
C.letter(128, 99, "d")


# ============================================================ e  per cohort
Cc = []
for c in cols:
    g = R[R.cohort == c]
    Cc.append(dict(cohort=c, tum=fs.ICB_TUM[c], infl=100 * g[g.family == "Inflamed/IFN"].clear.mean(),
                   conf=100 * g[g.family == "Confounder"].clear.mean()))
Cc = pd.DataFrame(Cc).sort_values("infl").reset_index(drop=True)
ax = C.ax(30, 164, 34, 28)
for i, r in Cc.iterrows():
    ax.plot([r.conf, r.infl], [i, i], color=K.LGRAY, lw=1.0, zorder=1)
    ax.scatter(r.conf, i, s=12, c=K.ORANGE, lw=.3, ec="white", zorder=3)
    ax.scatter(r.infl, i, s=12, c=K.BLUE, lw=.3, ec="white", zorder=4)
ax.set_yticks(range(len(Cc))); ax.set_yticklabels([f"{K.COHORT[c]} ({t[:4] if t != 'melanoma' else 'mel'})" for c, t in zip(Cc.cohort, Cc.tum)], fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0); ax.set_ylim(-.7, len(Cc) - .3)
K.ref(ax, 5, axis="v")
ax.set_xlim(-4, 104); ax.set_xticks([0, 50, 100]); ax.set_xlabel("signatures clearing the null (%)")
K.key(ax, [(K.BLUE, "inflamed/IFN")], x=.02, y=1.09, s=12)
K.key(ax, [(K.ORANGE, "confounder")], x=.55, y=1.09, s=12)
K.frame(ax, left=False)
C.letter(1.5, 156, "e")
wins = int((Cc.infl > Cc.conf).sum())

# ============================================================ f  by tumor type
ax = C.ax(82, 164, 54, 28)
TT = [t for t in TUMS]
R["tum"] = R.cohort.map(fs.ICB_TUM)
for i, t in enumerate(TT):
    for j, f in enumerate(FAMORD):
        v = R[(R.tum == t) & (R.family == f)].groupby("signature").delta.median()
        xx = i + (j - 1) * .22
        ax.scatter(xx + K.swarm(v.to_numpy(), .09, 12), v, s=3.5, c=K.tint(K.FAM[f], .25), lw=0, zorder=2)
        K.band(ax, xx + .085, v=v, color=K.FAM[f], width=.5, dot_s=6, z=4)
ax.axhline(0, color=K.INK, lw=.5)
for i in range(len(TT)):
    ax.axvspan(i - .5, i + .5, color=K.GRID if i % 2 == 0 else "white", lw=0, zorder=0)
ax.set_xticks(range(len(TT)))
ax.set_xticklabels([f"{t}\n({R[R.tum == t].cohort.nunique()})" for t in TT], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-.5, len(TT) - .5)
ax.set_ylim(-.16, .30); ax.set_yticks([-.1, 0, .1, .2, .3])
ax.set_ylabel("median ΔAUC per signature")
K.frame(ax)
K.key(ax, [(K.FAM[f], K.FAMLAB[f]) for f in FAMORD][:1], x=.02, y=1.09, s=12)
K.key(ax, [(K.FAM["Other"], "other")], x=.30, y=1.09, s=12)
K.key(ax, [(K.FAM["Confounder"], "confounder control")], x=.50, y=1.09, s=12)
C.letter(71, 156, "f")

# ============================================================ g  volcano
ax = C.ax(152, 164, 26, 28); C.letter(142, 156, "g")
lp = -np.log10(sigp)
for f in ("Other", "Confounder", "Inflamed/IFN"):
    s_ = [x_ for x_ in sigp.index if famd[x_] == f]
    K.dots(ax, med[s_], lp[s_], K.FAM[f], s=4 + 2.2 * np.log2(X.n_genes[s_]), alpha=0.9)
ax.axhspan(0, -np.log10(.05), color=K.GRID, lw=0, zorder=0); ax.axhline(-np.log10(.05), color=K.INK2, lw=0.5, ls=(0, (2, 1.5)), zorder=1); K.ref(ax, 0, axis="v")
ax.text(0.125, -np.log10(.05) + 0.08, "P = 0.05", ha="right", va="bottom", fontsize=K.T_SMALL, color=K.INK2)
top4 = lp.sort_values(ascending=False).index[:4]
for k, s_ in enumerate(top4):
    ax.annotate(K.lab(s_), (med[s_], lp[s_]), xytext=(-0.026, lp.max() * 1.06 - k * 0.32), fontsize=4.3, ha="left",
                va="center", color=K.INK, arrowprops=dict(arrowstyle="-", lw=0.25, color=K.MUTED, shrinkA=0, shrinkB=1.5))
ax.set_xlim(-0.03, 0.13); ax.set_xticks([0, 0.05, 0.1]); ax.set_ylim(0, lp.max() * 1.12)
ax.set_xlabel("median ΔAUC"); ax.set_ylabel("−log10 P (signed rank)")
K.frame(ax)
C.head(147, 156.2, "Effect and evidence")

# ============================================================ h  raincloud by family
ax = C.ax(30, 210, 50, 22); C.letter(1.5, 202, "h"); C.head(6.5, 202.2, "Every test, by family")
for i, f in enumerate(FAMORD):
    v = R[R.family == f].delta.to_numpy()
    K.violin(ax, v, i, K.FAM[f], width=0.42, side="right", alpha=0.55)
    ax.scatter(i - 0.12 - np.random.default_rng(i).uniform(0, 0.22, len(v)), v, s=1.4, color=K.FAM[f], alpha=0.6, lw=0, zorder=2)
K.ref(ax, 0)
ax.set_xticks(range(3)); ax.set_xticklabels([f"{K.FAMLAB[f]} ({(R.family == f).sum()})" for f in FAMORD], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-0.5, 2.8)
ax.set_ylim(-0.32, 0.42); ax.set_yticks([-0.2, 0, 0.2, 0.4]); ax.set_ylabel("ΔAUC")
K.frame(ax)

# ============================================================ i  share clearing vs cohort size
ax = C.ax(100, 210, 32, 22); C.letter(90, 202, "i"); C.head(95, 202.2, "Clearing and cohort size")
nco = R.groupby("cohort").n.first(); shr = R.groupby("cohort").clear.mean() * 100
K.dots(ax, nco[cols], shr[cols], K.INK, s=10)
OFF = {"Gide": (1.08, 0, "left"), "Puch": (1.08, 0, "left"), "Jung": (0.92, 0, "right"), "Kim": (1.08, 0, "left"),
       "Limagne": (1.08, 0, "left"), "Mariathasan": (1.0, 6.5, "center"), "Liu": (1.08, 0, "left"), "Braun": (1.08, 1.5, "left"),
       "Hugo": (0.92, 0, "right"), "Miao": (0.92, 1.0, "right"), "Riaz": (1.08, 3.2, "left"), "Van Allen": (1.08, -1.2, "left")}
for c in cols:
    fx, dy, ha = OFF[K.COHORT[c]]
    ax.text(nco[c] * fx, shr[c] + dy, K.COHORT[c], fontsize=4.3, va="center", ha=ha, color=K.INK2)
K.ref(ax, 5)
ax.set_xscale("log"); ax.set_xticks([30, 60, 120, 240]); ax.set_xticklabels(["30", "60", "120", "240"]); ax.minorticks_off()
ax.set_xlim(19, 320); ax.set_ylim(-2, 72); ax.set_yticks([0, 25, 50])
ax.set_xlabel("patients"); ax.set_ylabel("% of 52 clearing")
K.frame(ax)

# ============================================================ j  gain vs size
ax = C.ax(150, 210, 28, 22); C.letter(140, 202, "j"); C.head(145, 202.2, "Gain and signature size")
ng_ = X.n_genes
for f in ("Other", "Confounder", "Inflamed/IFN"):
    s_ = [x_ for x_ in med.index if famd[x_] == f]
    K.dots(ax, ng_[s_], med[s_], K.FAM[f], s=7, alpha=0.9)
K.ref(ax, 0)
ax.set_xscale("log"); ax.set_xticks([3, 10, 30, 100, 300]); ax.set_xticklabels(["3", "10", "30", "100", "300"]); ax.minorticks_off()
ax.set_xlabel("genes"); ax.set_ylabel("median ΔAUC"); ax.set_ylim(-0.03, 0.13); ax.set_yticks([0, 0.05, 0.1])
K.frame(ax)

C.save("Figure3_GB10")

k_ = {f: (int(R[R.family == f].clear.sum()), len(R[R.family == f])) for f in FAMORD}
print(f"Figure3_GB10 | tests {len(R)} | {nsig}/52 at P<0.05 | perm null median {ph.perm_null_median:.0f} "
      f"p95 {ph.perm_null_p95:.0f} max {ph.perm_null_max:.0f} P {ph.perm_P:.4f}")
for f in FAMORD:
    k, n = k_[f]; lo, hi = proportion_confint(k, n, method="wilson")
    print(f"  {f}: {k}/{n} = {100*k/n:.0f}% ({100*lo:.0f}-{100*hi:.0f}%); family median gain "
          f"{med[[s for s in med.index if famd[s]==f]].median():+.3f}")
brk = sigp[sigp < .05].index.map(famd).value_counts().to_dict()
print(f"  e: inflamed > confounder in {wins}/12 cohorts; f medians", {t: {f: round(R[(R.tum == t) & (R.family == f)].groupby("signature").delta.median().median(), 3) for f in FAMORD} for t in TT})
print(f"  significant by family {brk}; AUC median {ma:.3f} vs random {mr:.3f}; best signature "
      f"{R.groupby('signature').AUC.median().max():.3f}")
