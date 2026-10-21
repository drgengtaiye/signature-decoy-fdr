"""Figure 3 (v6) - immune signatures predict response to checkpoint blockade.

a  every signature in every checkpoint cohort against its matched null (hero)
b  how often each family clears the null
c  the number of significant signatures against a permutation null that keeps their correlation
d  absolute discrimination: published signatures against their matched random sets
e  per cohort, how often inflamed/IFN and confounder signatures clear the null
f  median gain by tumor type within each family (descriptive)
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
import v4kit as K
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

C = K.Canvas(218)

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
MX, MY, MW, MH = 38, 30, 63, 120
LIM = .15
ax = C.ax(MX, MY, MW, MH)
ax.imshow(M.values, aspect="auto", cmap=fs.DIV, norm=Normalize(-LIM, LIM), interpolation="nearest")
for i in range(nr + 1): ax.axhline(i - .5, color="white", lw=.3)
for j in range(nc + 1): ax.axvline(j - .5, color="white", lw=.3)
ax.set_yticks(range(nr)); ax.set_yticklabels([K.lab(s) for s in M.index], fontsize=K.T_SMALL)
for t, s in zip(ax.get_yticklabels(), M.index):
    t.set_color(fs.INK if famd[s] != "Other" else fs.SLATE)
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
    st.add_patch(Rectangle((0, start - .5), 1, n, fc=fs.FAM[f], ec="white", lw=.5))
    st.text(-.9, start + n / 2 - .5, fs.FAMLAB[f].replace(" control", ""), rotation=90,
            ha="center", va="center", fontsize=K.T_SMALL, color=fs.FAM[f], fontweight="bold")
    start += n
# tumor strip
tt = C.ax(MX, MY - 2.6, MW, 1.9); tt.set_xlim(-.5, nc - .5); tt.set_ylim(0, 1); tt.axis("off")
for j, c in enumerate(cols):
    tt.add_patch(Rectangle((j - .5, 0), 1, 1, fc=fs.TUMCOL[fs.ICB_TUM[c]], ec="white", lw=.6))
# top: share of signatures clearing the null, per cohort
top = C.ax(MX, MY - 12.5, MW, 8.2)
pc = [100 * R[R.cohort == c].clear.mean() for c in cols]
top.bar(range(nc), pc, width=.72, color=fs.SLATE, lw=0)
top.set_xlim(-.5, nc - .5); top.set_ylim(0, 60); top.set_yticks([0, 50]); top.set_xticks([])
top.set_ylabel("% clear", fontsize=K.T_SMALL, labelpad=1.5)
top.tick_params(axis="y", labelsize=K.T_SMALL)
K.frame(top)
# right: median gain per signature
rb = C.ax(MX + MW + 2, MY, 17, MH)
rb.barh(range(nr), med[M.index].values, height=.74, color=[fs.FAM[famd[s]] for s in M.index], lw=0)
for i, s in enumerate(M.index):
    if sigp[s] < .05:
        rb.text(max(med[s], 0) + .004, i, "*", fontsize=K.T_NOTE, va="center", ha="left")
rb.set_ylim(nr - .5, -.5); rb.set_yticks([]); rb.axvline(0, color=fs.INK, lw=.45)
rb.set_xlim(-.015, .13); rb.set_xticks([0, .1]); rb.set_xticklabels(["0", "0.1"])
rb.xaxis.tick_top(); rb.xaxis.set_label_position("top")
rb.set_xlabel("median ΔAUC", fontsize=K.T_SMALL, labelpad=2)
rb.tick_params(axis="x", labelsize=K.T_SMALL)
rb.spines["left"].set_visible(False); rb.spines["bottom"].set_visible(False)
rb.spines["top"].set_visible(True); rb.spines["right"].set_visible(False)
# colour bar and tumor key, in the corner above the row labels
cax = C.ax(9, 9.5, 22, 1.8)
cb = C.fig.colorbar(plt.cm.ScalarMappable(Normalize(-LIM, LIM), fs.DIV), cax=cax,
                    orientation="horizontal", ticks=[-.15, 0, .15])
cb.outline.set_linewidth(.3); cax.tick_params(labelsize=K.T_SMALL, length=1.5, width=.4, pad=1)
cax.set_title("ΔAUC vs matched null", fontsize=K.T_SMALL, pad=2)
kx = C.ax(MX + MW + 3, 3.5, 18, 13); kx.axis("off"); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
for i, t in enumerate(TUMS):
    yy = .92 - i * .21
    kx.add_patch(Rectangle((0, yy - .07), .14, .14, fc=fs.TUMCOL[t], ec="none"))
    kx.text(.22, yy, t, fontsize=K.T_SMALL, va="center")
C.letter(1.5, 3, "a")

# ============================================================ b  clearance by family
ax = C.ax(141, 14, 36, 40)
PSIG = R.groupby(["family", "signature"]).clear.mean().mul(100).reset_index()
for i, f in enumerate(FAMORD):
    g = R[R.family == f]; k, n = int(g.clear.sum()), len(g)
    lo, hi = proportion_confint(k, n, method="wilson")
    v = PSIG.loc[PSIG.family == f, "clear"].to_numpy()          # each signature's own clearance rate
    ax.scatter(i - .14 + K.swarm(v, .2, 14) * .8, v, s=6, c=fs.FAM[f], lw=.25, ec="white", alpha=.75, zorder=3)
    ax.plot([i + .2, i + .2], [100 * lo, 100 * hi], color=fs.FAM[f], lw=1.3, solid_capstyle="butt", zorder=4)
    ax.scatter([i + .2], [100 * k / n], s=26, c="white", ec=fs.FAM[f], lw=1.2, zorder=5)
ax.axhline(5, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.text(2.52, 4.2, "chance (5%)", ha="right", va="top", fontsize=K.T_SMALL, color=fs.SLATE)
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
ax = C.ax(141, 78, 36, 24)
DR = pd.read_csv(f"{D}/icb/v4_permutation_null_draws.csv.gz")
dv = DR[(DR.analysis == "12 cohorts (headline)") & (DR.statistic == "signatures with P<0.05")].value.to_numpy()
kk = np.arange(0, 31); fr = np.array([(dv == k_).mean() for k_ in kk]) * 100
ax.bar(kk, np.where(fr > 0, fr, np.nan), 1.0, color=[fs.SLATE if k_ < nsig else fs.BLUE for k_ in kk], lw=0, alpha=.85, zorder=2)
pmf = binom.pmf(kk, 52, .05) * 100                                   # what 52 independent tests would give
ax.plot(kk[pmf > 1e-2], pmf[pmf > 1e-2], color=fs.VERM, lw=.8, marker="o", ms=1.8, zorder=3)
ax.axvline(nsig, color=fs.BLUE, lw=1.1, zorder=4)
ax.text(nsig + .6, 30, f"observed {nsig}\n$P$ = {ph.perm_P:.3f}", ha="left", va="top", fontsize=K.T_SMALL, color=fs.BLUE)
ax.scatter([.02], [1.1], s=10, c=fs.SLATE, marker="s", lw=0, transform=ax.transAxes, clip_on=False)
ax.text(.06, 1.1, "permutation null", transform=ax.transAxes, fontsize=K.T_SMALL, va="center")
ax.plot([.52, .6], [1.1, 1.1], color=fs.VERM, lw=.8, transform=ax.transAxes, clip_on=False)
ax.text(.63, 1.1, "binomial", transform=ax.transAxes, fontsize=K.T_SMALL, va="center")
ax.set_yscale("log"); ax.set_ylim(.01, 60); ax.set_yticks([.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"]); ax.minorticks_off()
ax.set_xlim(-1, 30); ax.set_xticks([0, 10, 20, 30])
ax.set_xlabel("signatures with $P$ < 0.05 (of 52)"); ax.set_ylabel("permutations (%)")
K.frame(ax)
C.letter(128, 70, "c")

# ============================================================ d  published vs random
ax = C.ax(141, 122, 36, 38)
bins = np.arange(.30, .96, .025)
h1, _ = np.histogram(R.AUC, bins); h0, _ = np.histogram(R.rand_med, bins)
ctr = (bins[:-1] + bins[1:]) / 2
ax.bar(ctr, h1, width=.022, color=fs.BLUE, lw=0)
ax.bar(ctr, -h0, width=.022, color=fs.GREY, lw=0)
ax.axhline(0, color=fs.INK, lw=.45)
ma, mr = R.AUC.median(), R.rand_med.median()
ax.plot([ma] * 2, [0, h1.max() * 1.04], color=fs.INK, lw=.6, ls=(0, (2, 1.5)))
ax.plot([mr] * 2, [0, -h0.max() * 1.02], color=fs.INK, lw=.6, ls=(0, (2, 1.5)))
ax.text(.95, h1.max() * .9, "published", ha="right", va="top", fontsize=K.T_SMALL, color=fs.BLUE)
ax.text(.95, -h0.max() * .75, "matched\nrandom", ha="right", va="top", fontsize=K.T_SMALL,
        color=fs.SLATE)
ax.set_xlim(.33, .96); ax.set_ylim(-h0.max() * 1.08, h1.max() * 1.12)
yt = [-200, -100, 0, 100]
ax.set_yticks(yt); ax.set_yticklabels([str(abs(v)) for v in yt])
ax.set_xticks([.4, .6, .8])
ax.set_xlabel("AUC for response"); ax.set_ylabel("signature × cohort tests")
K.frame(ax)
C.letter(128, 114, "d")


# ============================================================ e  per cohort
Cc = []
for c in cols:
    g = R[R.cohort == c]
    Cc.append(dict(cohort=c, tum=fs.ICB_TUM[c], infl=100 * g[g.family == "Inflamed/IFN"].clear.mean(),
                   conf=100 * g[g.family == "Confounder"].clear.mean()))
Cc = pd.DataFrame(Cc).sort_values("infl").reset_index(drop=True)
ax = C.ax(30, 180, 44, 34)
for i, r in Cc.iterrows():
    ax.plot([r.conf, r.infl], [i, i], color=fs.LIGHT, lw=1.0, zorder=1)
    ax.scatter(r.conf, i, s=12, c=fs.VERM, lw=.3, ec="white", zorder=3)
    ax.scatter(r.infl, i, s=12, c=fs.BLUE, lw=.3, ec="white", zorder=4)
ax.set_yticks(range(len(Cc))); ax.set_yticklabels([c.replace("ICB_", "").replace("_", " ") for c in Cc.cohort], fontsize=5.0)
for t, tu in zip(ax.get_yticklabels(), Cc.tum): t.set_color(fs.TUMCOL[tu])
ax.tick_params(axis="y", length=0); ax.set_ylim(-.7, len(Cc) - .3)
ax.axvline(5, color=fs.INK, lw=.5, ls=(0, (3, 2)), zorder=0)
ax.set_xlim(-4, 104); ax.set_xticks([0, 50, 100]); ax.set_xlabel("signatures clearing the null (%)")
K.key(ax, [(fs.BLUE, "inflamed/IFN")], x=.02, y=1.07, s=12)
K.key(ax, [(fs.VERM, "confounder control")], x=.45, y=1.07, s=12)
K.frame(ax, left=False)
C.letter(1.5, 172, "e")
wins = int((Cc.infl > Cc.conf).sum())

# ============================================================ f  by tumor type
ax = C.ax(104, 180, 72, 34)
TT = [t for t in TUMS]
R["tum"] = R.cohort.map(fs.ICB_TUM)
for i, t in enumerate(TT):
    for j, f in enumerate(FAMORD):
        v = R[(R.tum == t) & (R.family == f)].groupby("signature").delta.median()
        xx = i + (j - 1) * .22
        ax.scatter(np.full(len(v), xx) + np.random.default_rng(i * 3 + j).uniform(-.06, .06, len(v)), v,
                   s=5, c=fs.FAM[f], lw=0, alpha=.45, zorder=2)
        ax.plot([xx - .09, xx + .09], [v.median()] * 2, color=fs.FAM[f], lw=2.0, zorder=4, solid_capstyle="butt")
ax.axhline(0, color=fs.INK, lw=.5)
ax.set_xticks(range(len(TT)))
ax.set_xticklabels([f"{t}\n({R[R.tum == t].cohort.nunique()})" for t in TT], fontsize=K.T_SMALL)
for tl, t in zip(ax.get_xticklabels(), TT): tl.set_color(fs.TUMCOL[t])
ax.tick_params(axis="x", length=0); ax.set_xlim(-.5, len(TT) - .5)
ax.set_ylim(-.16, .30); ax.set_yticks([-.1, 0, .1, .2, .3])
ax.set_ylabel("median ΔAUC per signature")
K.frame(ax)
K.key(ax, [(fs.FAM[f], fs.FAMLAB[f]) for f in FAMORD][:1], x=.02, y=1.07, s=12)
K.key(ax, [(fs.FAM["Other"], "other")], x=.30, y=1.07, s=12)
K.key(ax, [(fs.FAM["Confounder"], "confounder control")], x=.50, y=1.07, s=12)
C.letter(92, 172, "f")

C.save("Figure3")

k_ = {f: (int(R[R.family == f].clear.sum()), len(R[R.family == f])) for f in FAMORD}
print(f"Figure3 v6 | tests {len(R)} | {nsig}/52 at P<0.05 | perm null median {ph.perm_null_median:.0f} "
      f"p95 {ph.perm_null_p95:.0f} max {ph.perm_null_max:.0f} P {ph.perm_P:.4f}")
for f in FAMORD:
    k, n = k_[f]; lo, hi = proportion_confint(k, n, method="wilson")
    print(f"  {f}: {k}/{n} = {100*k/n:.0f}% ({100*lo:.0f}-{100*hi:.0f}%); family median gain "
          f"{med[[s for s in med.index if famd[s]==f]].median():+.3f}")
brk = sigp[sigp < .05].index.map(famd).value_counts().to_dict()
print(f"  e: inflamed > confounder in {wins}/12 cohorts; f medians", {t: {f: round(R[(R.tum == t) & (R.family == f)].groupby("signature").delta.median().median(), 3) for f in FAMORD} for t in TT})
print(f"  significant by family {brk}; AUC median {ma:.3f} vs random {mr:.3f}; best signature "
      f"{R.groupby('signature').AUC.median().max():.3f}")
