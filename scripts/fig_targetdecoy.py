"""Target-decoy false discovery rate for the signature library.

a  expected true and false discoveries for each null x statistic, FDR above each bar
b  decoys reproduce their target's co-expression cohort by cohort
c  how often decoys are called, at the paper's own unit (one Wilcoxon per module over 12 cohorts)
d  calibrated per test is not calibrated per library (folded statistic, published null)
e  why: co-expression widens the spread of a module's AUC
f  each signature against its own 32 decoys (directed, published null)
g  matched nulls drawn from a signature-enriched pool favour one declared direction
h  reading a decoy against its target's bank or against its own: the check of the shortcut

Intervals: cluster bootstrap over the 52 targets, each carried with its decoys.
"""
import glob, os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D + "/icb/"
RNG = np.random.default_rng(20261019)
C_FOLD, C_DIR, C_SIG, C_DEC = fs.GREY, fs.BLUE, fs.TEAL, fs.SLATE
TT = 6.9

U = pd.read_csv(D + "calib_v3_units.csv")
F = pd.read_csv(D + "calib_v3_fdr.csv")
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(D + "calib_v3/ICB_*_v3.csv"))])
DIRN = X.drop_duplicates("unit").set_index("unit").dirn
U["dirn"] = U.unit.map(DIRN); U["tdirn"] = U.target.map(DIRN)
OWN = pd.read_csv(D + "calib_v3_ownbank_compare.csv") if os.path.exists(D + "calib_v3_ownbank_compare.csv") else None
NULLS = [n for n in ("N0", "N1", "N2", "N1nat", "N2nat") if n in set(F["null"])]
COMBO = [(n, s) for n in NULLS for s in ("fold", "dir")]
NLAB = {"N0": "N0", "N1": "N1", "N2": "N2", "N1nat": "N1", "N2nat": "N2"}

def boot_rate(df, col="call", nb=4000):
    per = df.groupby("target")[col].mean().to_numpy(float)
    bs = [per[RNG.integers(0, len(per), len(per))].mean() for _ in range(nb)]
    return 100 * per.mean(), 100 * np.percentile(bs, 2.5), 100 * np.percentile(bs, 97.5)

C = K.Canvas(208)

# ================================================================= a  headline
ax = C.ax(15, 17, 99, 32)
fa = F[(F.pool == "A") & (F.decoys == "coex")].set_index(["null", "stat"])
xs, w = [], .72
gap = {"N0": 0, "N1": .5, "N2": .5, "N1nat": 1.0, "N2nat": 1.0}
for i, (nl, st) in enumerate(COMBO):
    x = i + gap[nl]; xs.append(x); r = fa.loc[(nl, st)]
    cc = C_FOLD if st == "fold" else C_DIR
    et, ef = r.called - r.exp_false, r.exp_false
    elo, ehi = r.decoy_lo * 52 / 100, r.decoy_hi * 52 / 100
    ax.bar(x, et, w, color=cc, lw=0, zorder=3)
    ax.bar(x, ef, w, bottom=et, color="white", edgecolor=cc, hatch="/////", lw=.6, zorder=3)
    ax.errorbar(x, et, yerr=[[max(0, ehi - ef)], [max(0, ef - elo)]], fmt="none", ecolor=fs.INK,
                lw=.6, capsize=1.4, zorder=4)
    ax.text(x, r.called + .8, f"{r.fdr:.0f}%", ha="center", va="bottom", fontsize=K.T_SMALL,
            color=cc, fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([("folded" if s == "fold" else "directed") for _, s in COMBO], fontsize=5.0)
ax.tick_params(axis="x", length=0, pad=1.5)
tx = ax.get_xaxis_transform()
for nl in NULLS:
    ii = [x for x, (n, _) in zip(xs, COMBO) if n == nl]
    ax.text(np.mean(ii), -.115, NLAB[nl], transform=tx, ha="center", va="top", fontsize=K.T_SMALL,
            color=fs.INK, fontweight="bold")
for grp, lab in [(("N0",), "published null"), (("N1", "N2"), "enriched pool"),
                 (("N1nat", "N2nat"), "natural pool")]:
    ii = [x for x, (n, _) in zip(xs, COMBO) if n in grp]
    if not ii: continue
    ax.plot([min(ii) - .36, max(ii) + .36], [-.215, -.215], transform=tx, color=fs.SLATE, lw=.5, clip_on=False)
    ax.text(np.mean(ii), -.235, lab, transform=tx, ha="center", va="top", fontsize=K.T_SMALL, color=fs.SLATE)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=fs.SLATE, lw=0, label="expected true"),
                   Patch(facecolor="white", edgecolor=fs.SLATE, hatch="/////", lw=.6, label="expected false"),
                   Patch(facecolor="none", lw=0, label="%  target-decoy FDR")],
          fontsize=K.T_SMALL, loc="upper left", ncol=3, handlelength=1.2, borderpad=.15, columnspacing=1.0,
          handletextpad=.4, frameon=False)
ax.set_ylabel("signatures called (of 52)")
ax.set_ylim(0, 40); ax.set_yticks([0, 10, 20, 30, 40]); ax.set_xlim(min(xs) - .7, max(xs) + .7)
K.frame(ax)
C.title(15, 11, "Directed scoring finds more signatures at a lower false discovery rate",
        size=TT, fontweight="bold")
C.letter(7.5, 10.5, "a")

# ================================================================= b  matching
ax = C.ax(131, 17, 34, 32)
x0 = X[X["null"] == "N0"]
tr = x0[x0.kind == "target"].set_index(["target", "cohort"]).r_int
for (p, s), cc, lab in [(("A", "rand"), C_DEC, "random-gene"), (("A", "coex"), C_SIG, "co-expressed")]:
    dd = x0[(x0.kind == "decoy") & (x0.pool == p) & (x0.structure == s)].groupby(["target", "cohort"]).r_int.median()
    j = pd.concat([tr, dd], axis=1, keys=["t", "d"]).dropna()
    ax.scatter(j.t, j.d, s=3, color=cc, alpha=.55, lw=0, zorder=3, label=lab)
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, .3, .6, .9]); ax.set_yticks([0, .3, .6, .9])
ax.set_xlabel("signature internal $r$"); ax.set_ylabel("its decoys' internal $r$")
K.key(ax, [(C_SIG, "co-expressed decoys"), (C_DEC, "random-gene decoys")], x=.04, y=.95, s=9, size=K.T_SMALL)
fs.note(ax, .04, .80, "one point per\nsignature × cohort", size=K.T_SMALL, color=fs.SLATE)
K.frame(ax)
C.title(131, 11, "Matched in every cohort", size=TT, fontweight="bold")
C.letter(123.5, 10.5, "b")

# ================================================================= c  heatmap of decoy call rates
ax = C.ax(42, 73, 72, 21)
ROWS = [("A", "coex", "co-expressed decoys"), ("B", "coex", "co-expressed, no TME gene"),
        ("A", "rand", "random-gene decoys")]
M = np.full((len(ROWS), len(COMBO)), np.nan)   # the 208-decoy own-bank check is panel h
for i, (p, s, _) in enumerate(ROWS):
    for j, (nl, st) in enumerate(COMBO):
        d = U[(U.kind == "decoy") & (U.pool == p) & (U.structure == s) & (U["null"] == nl) & (U.stat == st)]
        M[i, j] = 100 * d.call.mean()
rl = [r[2] for r in ROWS]
cm = LinearSegmentedColormap.from_list("td", ["#2F6DB5", "#F4F6F8", "#E9A07E", "#C1401A", "#7A1F0A"])
nrm = TwoSlopeNorm(vmin=0, vcenter=2.5, vmax=20)
ax.imshow(M, cmap=cm, norm=nrm, aspect="auto", interpolation="nearest")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isfinite(M[i, j]):
            ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=K.T_SMALL,
                    color="white" if M[i, j] > 11 or M[i, j] < .6 else fs.INK)
        else:
            ax.text(j, i, "–", ha="center", va="center", fontsize=K.T_SMALL, color=fs.GREY)
ax.set_yticks(range(len(rl))); ax.set_yticklabels(rl, fontsize=K.T_SMALL)
ax.set_xticks(range(len(COMBO)))
ax.set_xticklabels([f"{NLAB[n]}{'·nat' if 'nat' in n else ''}\n{'fold' if s=='fold' else 'dir'}" for n, s in COMBO], fontsize=K.T_SMALL)
ax.tick_params(length=0)
for x in np.arange(1.5, len(COMBO) - 1, 2):
    ax.axvline(x, color="white", lw=1.6)
for sp in ax.spines.values(): sp.set_visible(False)
cax = C.ax(42, 102.5, 26, 1.8)
cb = plt.colorbar(plt.cm.ScalarMappable(norm=nrm, cmap=cm), cax=cax, orientation="horizontal")
cb.set_ticks([0, 2.5, 10, 20]); cb.ax.tick_params(labelsize=K.T_SMALL, length=1.5); cb.outline.set_visible(False)
C.text(70, 103.6, "% of decoys called; nominal 2.5% (two-sided P < 0.05, positive)", size=K.T_SMALL, va="center")
C.title(15, 67, "Decoys called at the paper's own unit of analysis", size=TT, fontweight="bold")
C.letter(7.5, 66.5, "c")

# ================================================================= e  spread
ax = C.ax(131, 73, 34, 21)
x0 = X[X["null"] == "N0"].copy()
x0["signed"] = np.where(x0.dirn == 1, x0.dir_obs, 1 - x0.dir_obs)
sp = x0[x0.kind == "decoy"].groupby(["cohort", "structure"]).signed.std().unstack()
sp = sp[sp.index.isin(X.cohort.unique())]
ax.scatter(sp.rand, sp.coex, s=10, color=C_SIG, lw=.4, edgecolor="white", zorder=3)
ANN = {"ICB_Miao1": (3, -2, "left"), "ICB_Jung": (-3, 2, "right"), "ICB_Kim": (3, 0, "left"),
       "ICB_Van_Allen": (3, 0, "left"), "ICB_Mariathasan": (3, 0, "left")}
for c, r in sp.iterrows():
    if c in ANN:
        dx, dy, ha = ANN[c]
        ax.annotate(c.replace("ICB_", "").replace("_", " "), (r.rand, r.coex), xytext=(dx, dy), textcoords="offset points",
                    fontsize=5.0, color=fs.SLATE, va="center", ha=ha)
lo_, hi_ = .04, .14
ax.plot([lo_, hi_], [lo_, hi_], color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_xticks([.05, .1]); ax.set_yticks([.05, .1])
ax.set_xlabel("sd of AUC, random-gene decoys"); ax.set_ylabel("sd, co-expressed")
K.frame(ax)
C.title(131, 67, "Co-expression widens the spread", size=TT, fontweight="bold")
C.letter(123.5, 66.5, "e")

# ================================================================= d  per test vs per library
ax = C.ax(15, 121, 41, 26)
x0 = X[(X["null"] == "N0") & (X.kind == "decoy") & (X.pool == "A")].copy()
x0["clear"] = (x0.fold_obs > x0.fold_p95).astype(float)
vals = {}
for s in ("rand", "coex"):
    vals[("test", s)] = boot_rate(x0[x0.structure == s], "clear")
    vals[("lib", s)] = boot_rate(U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == s) & (U["null"] == "N0") & (U.stat == "fold")])
for gi, (g, nom, lab) in enumerate([("test", 5, "each test:\nabove the null 95th pct"), ("lib", 2.5, "each module:\nWilcoxon over 12 cohorts")]):
    for si, (s, cc) in enumerate([("rand", C_DEC), ("coex", C_SIG)]):
        v = vals[(g, s)]; x = gi * 1.0 + (si - .5) * .36
        ax.bar(x, v[0], .34, color=cc, lw=0, zorder=3)
        ax.errorbar(x, v[0], yerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt="none", ecolor=fs.INK, lw=.6, capsize=1.4, zorder=4)
        ax.text(x, v[2] + .25, f"{v[0]:.1f}", ha="center", va="bottom", fontsize=K.T_SMALL, color=cc, fontweight="bold")
    ax.plot([gi - .42, gi + .42], [nom, nom], color=fs.INK, lw=.6, ls=(0, (3, 2)), zorder=5)
    ax.text(gi + .44, nom, f"{nom:g}%", fontsize=K.T_SMALL, va="center", ha="left", color=fs.INK)
ax.set_xticks([0, 1]); ax.set_xticklabels(["each test:\nabove null 95th pct", "each module:\nWilcoxon, 12 cohorts"], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0)
ax.set_ylabel("decoys called (%)"); ax.set_ylim(0, 11); ax.set_xlim(-.6, 1.75)
K.key(ax, [(C_DEC, "random-gene decoys"), (C_SIG, "co-expressed decoys")], x=.04, y=.96, s=9, size=K.T_SMALL)
K.frame(ax)
C.title(15, 115, "Calibrated per test, not per library", size=TT, fontweight="bold")
C.letter(7.5, 114.5, "d")
DVAL = vals

# ================================================================= g  direction
ax = C.ax(73, 121, 41, 26)
gn = [n for n in NULLS]
for gi, nl in enumerate(gn):
    for si, (dn, cc, lab) in enumerate([(1, C_DIR, "sensitivity"), (-1, fs.VERM, "resistance")]):
        d = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == nl) & (U.stat == "dir") & (U.dirn == dn)]
        v = boot_rate(d); x = gi + (si - .5) * .36
        ax.bar(x, v[0], .34, color=cc, lw=0, zorder=3, label=lab if gi == 0 else None)
        ax.errorbar(x, v[0], yerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt="none", ecolor=fs.INK, lw=.6, capsize=1.2, zorder=4)
ax.axhline(2.5, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
ax.set_xticks(range(len(gn)))
ax.set_xticklabels([NLAB[n] for n in gn], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0, pad=1.5)
tx = ax.get_xaxis_transform()
for grp, lab in [(("N0",), "published"), (("N1", "N2"), "enriched pool"), (("N1nat", "N2nat"), "natural pool")]:
    ii = [i for i, n in enumerate(gn) if n in grp]
    if not ii: continue
    ax.plot([min(ii) - .36, max(ii) + .36], [-.13, -.13], transform=tx, color=fs.SLATE, lw=.5, clip_on=False)
    ax.text(np.mean(ii), -.15, lab, transform=tx, ha="center", va="top", fontsize=K.T_SMALL, color=fs.SLATE)
ax.set_ylabel("decoys called (%)"); ax.set_ylim(0, 26)
ax.legend(fontsize=K.T_SMALL, loc="upper right", handlelength=1.0, borderpad=.15, labelspacing=.3,
          handletextpad=.4, title="declared direction", title_fontsize=K.T_SMALL, frameon=False)
K.frame(ax)
C.title(73, 115, "An enriched pool favours one direction", size=TT, fontweight="bold")
C.letter(65.5, 114.5, "g")

# ================================================================= h  shared vs own bank
ax = C.ax(131, 121, 34, 26)
if OWN is not None:
    for nl, cc, mk, lab in [("N0", fs.GREY, "o", "N0"), ("N2", fs.VERM, "s", "N2 enriched"), ("N2nat", C_DIR, "^", "N2 natural")]:
        q = OWN[(OWN["null"] == nl) & (OWN.stat == "dir")].pivot_table(index="unit", columns="src", values="med").dropna()
        if not len(q): continue
        ax.scatter(q.shared, q.own, s=4, color=cc, marker=mk, lw=0, alpha=.7, zorder=3,
                   label=f"{lab}, r = {np.corrcoef(q.shared, q.own)[0, 1]:.2f}")
    lim = (-.12, .12)
    ax.plot(lim, lim, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xticks([-.1, 0, .1]); ax.set_yticks([-.1, 0, .1])
    ax.set_xticklabels([K.minus("-0.1"), "0", "0.1"]); ax.set_yticklabels([K.minus("-0.1"), "0", "0.1"])
    ax.legend(fontsize=K.T_SMALL, loc="upper left", handlelength=.8, borderpad=.15, labelspacing=.3, handletextpad=.3)
else:
    fs.note(ax, .5, .5, "pending", ha="center", va="center")
ax.set_xlabel("median ΔAUC, target's bank"); ax.set_ylabel("own bank")
K.frame(ax)
C.title(131, 115, "Target's bank vs its own", size=TT, fontweight="bold")
C.letter(123.5, 114.5, "h")

# ================================================================= f  per signature
ax = C.ax(15, 168, 150, 26)
t = U[(U.kind == "target") & (U["null"] == "N0") & (U.stat == "dir")].set_index("target")
d = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == "N0") & (U.stat == "dir")]
order = t.med_delta.sort_values(ascending=False).index.tolist()
for i, s in enumerate(order):
    v = d[d.target == s].med_delta.to_numpy()
    lo, hi = np.percentile(v, [5, 95])
    ax.plot([i, i], [lo, hi], color=fs.LIGHT, lw=2.2, solid_capstyle="butt", zorder=2)
    ax.plot([i - .3, i + .3], [np.median(v)] * 2, color=C_DEC, lw=.6, zorder=3)
    r = t.loc[s]
    fc = C_SIG if r.call else "white"; ec = fs.VERM if r.neg else C_SIG
    ax.scatter(i, r.med_delta, s=11, facecolor=fc, edgecolor=ec, lw=.7 if not r.neg else 1.0, zorder=4)
ax.axhline(0, color=fs.INK, lw=.5, zorder=1)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([K.lab(s) for s in order], rotation=90, fontsize=5.0)
ax.tick_params(axis="x", length=1.2, pad=1)
ax.set_xlim(-.8, len(order) - .2)
ax.set_ylabel("median ΔAUC,\ndirected, N0")
K.key(ax, [(C_SIG, "signature, called"), ("white", "signature, not called")], x=.80, y=.95, s=10, size=K.T_SMALL)
ax.scatter([.80], [.95 - .085], s=10, facecolor="white", edgecolor=C_SIG, lw=.7, transform=ax.transAxes, zorder=11, clip_on=False)
ax.scatter([.80], [.95 - .17], s=10, facecolor="white", edgecolor=fs.VERM, lw=1.0, transform=ax.transAxes, zorder=11, clip_on=False)
ax.text(.835, .95 - .17, "significantly negative", transform=ax.transAxes, va="center", fontsize=K.T_SMALL, color=fs.INK)
ax.plot([.785, .815], [.95 - .255] * 2, color=fs.LIGHT, lw=2.2, transform=ax.transAxes, clip_on=False)
ax.text(.835, .95 - .255, "its 32 decoys, 5–95%", transform=ax.transAxes, va="center", fontsize=K.T_SMALL, color=fs.INK)
K.frame(ax)
C.title(15, 162, "Each signature against its own decoys", size=TT, fontweight="bold")
C.letter(7.5, 161.5, "f")

C.save("Fig_targetdecoy")
print("Fig_targetdecoy |", ", ".join(f"{n}-{s} {fa.loc[(n,s)].called:.0f} FDR {fa.loc[(n,s)].fdr:.1f}" for n, s in COMBO))
print("  d:", {k: tuple(round(x, 1) for x in v) for k, v in DVAL.items()})
