"""Figure 4 (GB v10) - a target-decoy false discovery rate for 52 immune signatures and response to checkpoint blockade.

a  decoys reproduce their target's co-expression cohort by cohort (first decoy set)
b  signatures called, split into expected true and expected false discoveries, by null and statistic (first set)
c  calibrated per test, not per library (folded, N0)
d  the two statistics by the signatures' own co-expression (closest-matched decoys)
e  one decoy set against ten: the estimate of each set, pooled, and the first set
f  the estimate by decoy class, with intervals over signatures and over correlation clusters (ten sets)
g  where each signature sits among its own decoys: the decoy-referenced q behind pi0
h  each signature's median delta against its 320 co-expressed decoys (folded, ten sets)
i  signatures called and the estimate as the per-unit P threshold varies (ten sets; 0.05 is the one reported)
j  quantile-quantile plot of the per-unit signed-rank P: signatures, co-expressed and random-gene decoys
k  the decoys' call rate by their own internal correlation, under N0 and the natural-pool N1 and N2 (first set)
l  the estimate pooled over the first k of the ten decoy sets
Data: calib_v3_fdr.csv, calib_v3_units.csv, calib_v3/ICB_*_v3.csv, calib_v3_fdr_diff_close.csv (first set);
v12_pkg/response_fdr.csv, response_fdr_by_set.csv, response_calls.csv.gz (package, ten sets).
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import glob, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
D = K.D + "/icb/"; P = D + "v12_pkg/"
RNG = np.random.default_rng(20261019)

F = pd.read_csv(D + "calib_v3_fdr.csv")
U = pd.read_csv(D + "calib_v3_units.csv")
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(D + "calib_v3/ICB_*_v3.csv"))])
SC = pd.read_csv(D + "calib_v3_fdr_diff_close.csv")
RF = pd.read_csv(P + "response_fdr.csv"); RB = pd.read_csv(P + "response_fdr_by_set.csv")
RC = pd.read_csv(P + "response_calls.csv.gz")
FAM = pd.read_csv(D + "v5_signature_clusters.csv").set_index("signature").family


def boot_rate(df, col="call", nb=4000):
    per = df.groupby("target")[col].mean().to_numpy(float)
    bs = [per[RNG.integers(0, len(per), len(per))].mean() for _ in range(nb)]
    return 100 * per.mean(), 100 * np.percentile(bs, 2.5), 100 * np.percentile(bs, 97.5)


C = K.Canvas(210)

# ================================================================== a  matching
ax = C.panel("a", 14, 12, 32, 32, head="Decoys match co-expression")
x0 = X[X["null"] == "N0"]
tr = x0[x0.kind == "target"].set_index(["target", "cohort"]).r_int
for (p, s), cc in [(("A", "rand"), K.SLATE), (("A", "coex"), K.AQUA)]:
    dd = x0[(x0.kind == "decoy") & (x0.pool == p) & (x0.structure == s)].groupby(["target", "cohort"]).r_int.median()
    j = pd.concat([tr, dd], axis=1, keys=["t", "d"]).dropna()
    K.contour2d(ax, j.t, j.d, cc, (0, 1), (0, 1), levels=(0.5, 0.75, 0.9, 0.97), bw=0.18)
    ax.scatter(j.t, j.d, s=0.8, color=cc, alpha=0.5, lw=0, zorder=3)
K.identity(ax, 0, 1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
ax.set_xlabel("signature internal r"); ax.set_ylabel("decoys’ internal r (median)")
K.legend(ax, [(K.AQUA, "co-expressed"), (K.SLATE, "random-gene")], x=0.06, y=0.95, dy=0.09, s=8, size=K.T_SMALL)

# ================================================================== b  calls by null and statistic
ax = C.panel("b", 64, 12, 72, 32, head="Matching the null on co-expression halves the folded estimate")
fa = F[(F.pool == "A") & (F.decoys == "coex")].set_index(["null", "stat"])
COMBO = [(n, s) for n in ("N0", "N1", "N2", "N1nat", "N2nat") for s in ("fold", "dir")]
gap = {"N0": 0, "N1": .55, "N2": .55, "N1nat": 1.1, "N2nat": 1.1}
xs = []
for i, (nl, st) in enumerate(COMBO):
    x = i + gap[nl]; xs.append(x); r = fa.loc[(nl, st)]
    K.stacked_calls(ax, x, r.called, r.exp_false, K.STAT[st], width=0.66, text=f"{r.fdr:.0f}%")
ax.set_xticks(xs); ax.set_xticklabels(["F" if s == "fold" else "D" for _, s in COMBO]); ax.tick_params(axis="x", length=0)
tx = ax.get_xaxis_transform()
for nl in ("N0", "N1", "N2", "N1nat", "N2nat"):
    ii = [x for x, (n, _) in zip(xs, COMBO) if n == nl]
    ax.text(np.mean(ii), -0.13, nl[:2], transform=tx, ha="center", va="top", fontsize=K.T_SMALL, fontweight="bold")
for grp, lab in [(("N0",), "size, engine"), (("N1", "N2"), "enriched pool"), (("N1nat", "N2nat"), "natural pool")]:
    ii = [x for x, (n, _) in zip(xs, COMBO) if n in grp]
    ax.plot([min(ii) - .33, max(ii) + .33], [-0.235, -0.235], transform=tx, color=K.INK2, lw=0.5, clip_on=False)
    ax.text(np.mean(ii), -0.25, lab, transform=tx, ha="center", va="top", fontsize=K.T_SMALL, color=K.INK2)
ax.set_ylabel("signatures called (of 52)"); ax.set_ylim(0, 38); ax.set_xlim(min(xs) - .6, max(xs) + .6)
K.legend(ax, [(K.SLATE, "folded (F)", "box"), (K.BLUE, "directed (D)", "box"), (K.tint(K.SLATE, .62), "expected false", "box")],
         x=0.02, y=0.98, dy=0.085, s=8, size=K.T_SMALL)
ax.text(0.30, 0.98, "% estimated FDR", transform=ax.transAxes, fontsize=K.T_SMALL, va="center", color=K.INK2)

# ================================================================== c  per test vs per library
ax = C.panel("c", 152, 12, 26, 32, head="Per test vs per library")
xd = X[(X["null"] == "N0") & (X.kind == "decoy") & (X.pool == "A")].copy()
xd["clear"] = (xd.fold_obs > xd.fold_p95).astype(float)
for gi, (g, nom) in enumerate([("test", 5.0), ("lib", 2.5)]):
    ax.plot([gi - .42, gi + .42], [nom, nom], color=K.ORANGE, lw=0.7, ls=(0, (2, 1.4)), zorder=2)
    ax.text(gi + .44, nom, f"{nom:g}%", fontsize=K.T_SMALL - 0.3, ha="left", va="center", color=K.ORANGE)
    for si, (s, cc) in enumerate([("rand", K.SLATE), ("coex", K.AQUA)]):
        v = boot_rate(xd[xd.structure == s], "clear") if g == "test" else boot_rate(
            U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == s) & (U["null"] == "N0") & (U.stat == "fold")])
        x = gi + (si - .5) * .36
        ax.bar(x, v[0], .30, color=K.tint(cc, .55), edgecolor="none", zorder=3)
        K.vbar_interval(ax, x, v[1], v[2], cc, lw=1.0, z=4)
        K.dots(ax, [x], [v[0]], cc, s=11, z=5)
ax.set_xticks([0, 1]); ax.set_xticklabels(["each test", "each module"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("decoys called (%)"); ax.set_ylim(0, 10.5); ax.set_xlim(-.55, 1.75)
K.legend(ax, [(K.SLATE, "random-gene"), (K.AQUA, "co-expressed"), (K.ORANGE, "nominal", "line")], x=0.06, y=0.98, dy=0.09, s=7,
         size=K.T_SMALL - 0.3)

# ================================================================== d  by co-expression
ax = C.panel("d", 14, 64, 32, 32, head="Direction helps weak co-expression")
for gi, (sub, lab) in enumerate([("targets r <= 0.5", "r ≤ 0.5\n(30)"), ("targets r > 0.5", "r > 0.5\n(22)")]):
    q = SC[(SC.subset == sub) & (SC.alt == "N0-dir")]
    et = q[q.quantity == "exp_true"].iloc[0]; fd = q[q.quantity == "fdr"].iloc[0]
    for si, (called, true_, fdr_, st) in enumerate([(et.ref_called, et.ref_value, fd.ref_value, "fold"),
                                                    (et.alt_called, et.alt_value, fd.alt_value, "dir")]):
        x = gi + (si - .5) * .36
        K.stacked_calls(ax, x, called, called - true_, K.STAT[st], width=0.33, text=f"{fdr_:.0f}%")
ax.set_xticks([0, 1]); ax.set_xticklabels(["r ≤ 0.5 (30)", "r > 0.5 (22)"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("signatures called"); ax.set_ylim(0, 18); ax.set_xlim(-.55, 1.55)
ax.text(0.02, 0.98, "closest-matched decoys, N0", transform=ax.transAxes, fontsize=K.T_SMALL, color=K.INK2, va="top")

# ================================================================== e  ten decoy sets
ax = C.panel("e", 64, 64, 26, 32, head="Ten decoy sets")
F0 = F[(F["null"] == "N0") & (F.pool == "A") & (F.decoys == "coex")].set_index("stat")
for gi, st in enumerate(("fold", "dir")):
    b = RB[RB.stat == st].fdr.to_numpy()
    pooled = RF[(RF.stat == st) & (RF.decoys == "coex") & (RF.boot == "clusters")].iloc[0]
    K.dots(ax, gi + np.linspace(-.22, .22, len(b)), b, K.STAT[st], s=7, alpha=0.9)
    ax.plot([gi - .3, gi + .3], [pooled.fdr] * 2, color=K.INK, lw=1.0, zorder=5)
    K.dots(ax, [gi + .38], [F0.loc[st, "fdr"]], K.INK, s=12, filled=False)
ax.set_xticks([0, 1]); ax.set_xticklabels(["folded", "directed"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("estimated FDR (%)"); ax.set_ylim(0, 24); ax.set_xlim(-.55, 1.6)
K.legend(ax, [(K.SLATE, "each set"), (K.INK, "pooled", "line"), (K.INK, "first set", "o", False)], x=0.50, y=0.98, dy=0.085, s=8,
         size=K.T_SMALL)

# ================================================================== f  by decoy class
ax = C.panel("f", 106, 64, 30, 32, head="By decoy class")
rows = [("coex", "co-expressed"), ("coex closest", "closest-matched"), ("rand", "random-gene")]
for gi, st in enumerate(("fold", "dir")):
    for k, (dc, lab) in enumerate(rows):
        y = k + (gi - .5) * 0.32
        rc = RF[(RF.stat == st) & (RF.decoys == dc) & (RF.boot == "clusters")].iloc[0]
        rt = RF[(RF.stat == st) & (RF.decoys == dc) & (RF.boot == "targets")].iloc[0]
        K.hbar_interval(ax, y, rc.fdr_lo, min(rc.fdr_hi, 60), K.tint(K.STAT[st], .55), lw=0.8)
        K.hbar_interval(ax, y, rt.fdr_lo, rt.fdr_hi, K.STAT[st], lw=1.3)
        K.dots(ax, [rc.fdr], [y], K.STAT[st], s=10, z=6)
        K.dots(ax, [rc.fdr_pi0], [y], K.STAT[st], s=10, filled=False, z=6)
ax.set_yticks(range(3)); ax.set_yticklabels([r[1] for r in rows]); ax.set_ylim(2.6, -0.6); K.frame(ax, left=False)
ax.set_xlim(0, 60); ax.set_xlabel("estimated FDR (%)")


# ================================================================== g  q of each target among its decoys
ax = C.panel("g", 152, 64, 26, 32, head="π0 from the decoys")
qs = {}
for st in ("fold", "dir"):
    c = RC[RC.stat == st]
    t = c[c.kind == "target"].set_index("target")
    d = c[(c.kind == "decoy") & (c["class"] == "coex")]
    if st == "dir":
        d = d[d.dirn == d.target.map(t.dirn)]
    s_t = np.sign(t.med_delta) * -np.log10(t.p.clip(lower=1e-300))
    s_d = np.sign(d.med_delta) * -np.log10(d.p.clip(lower=1e-300))
    q = pd.Series({tg: (1 + (s_d[d.target == tg] >= s_t[tg] - 1e-12).sum()) / (1 + (d.target == tg).sum()) for tg in t.index})
    qs[st] = q
bins = np.linspace(0, 1, 11)
K.hist_steps(ax, qs["fold"], bins, K.SLATE, alpha=0.55, density_=False)
pi0 = RF[(RF.stat == "fold") & (RF.decoys == "coex") & (RF.boot == "clusters")].iloc[0].pi0
ax.axhline(pi0 * 52 / 10, color=K.ORANGE, lw=0.9, zorder=4)
ax.text(0.55, pi0 * 52 / 10 + 0.6, f"π0 = {pi0:.2f}", ha="center", va="bottom", fontsize=K.T_SMALL, color=K.INK)
ax.set_xlabel("q: rank among own decoys"); ax.set_ylabel("signatures (folded)")
ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1])
ax.text(0.98, 0.98, "null signatures\nspread evenly", transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL - 0.3,
        color=K.INK2, linespacing=1.05)

# ================================================================== h  each signature against its decoys
ax = C.panel("h", 14, 116, 164, 30, head="Each signature against its 320 co-expressed decoys (folded)")
c = RC[RC.stat == "fold"]
t = c[c.kind == "target"].set_index("target")
d = c[(c.kind == "decoy") & (c["class"] == "coex")]
qd = d.groupby("target").med_delta.quantile([.05, .25, .5, .75, .95]).unstack()
order = t.med_delta.sort_values(ascending=False).index
for i, s in enumerate(order):
    if i % 2 == 0:
        ax.axvspan(i - .5, i + .5, color=K.GRID, lw=0, zorder=0)
    K.band(ax, i, q=qd.loc[s].to_numpy(), color=K.AQUA, width=0.9, dot=False)
    ax.plot([i - .28, i + .28], [qd.loc[s, .5]] * 2, color=K.AQUA, lw=0.8, zorder=4)
    K.dots(ax, [i], [t.loc[s, "med_delta"]], K.FAM[FAM[s]], s=13 if t.loc[s, "call"] else 10, filled=bool(t.loc[s, "call"]), z=5)
K.ref(ax, 0)
ax.set_xticks(range(len(order))); ax.set_xticklabels([K.lab(s) for s in order], rotation=90, fontsize=4.2)
ax.tick_params(axis="x", length=0, pad=1); ax.set_xlim(-0.8, len(order) - 0.2)
ax.set_ylabel("median ΔAUC")
K.hlegend(C, 100, 113.0, [(K.BLUE, "inflamed/IFN"), (K.GRAY, "other"), (K.ORANGE, "confounder"), (K.INK, "called", "o", True),
                          (K.INK, "not called", "o", False), (K.tint(K.AQUA, .45), "decoys 5–95%", "s"), (K.tint(K.AQUA, .18), "25–75%", "s")], size=K.T_SMALL, s=7,
          gap=1.6)
# ================================================================== i  threshold curve
TH = np.exp(np.linspace(np.log(0.001), np.log(0.2), 40))
curves = {}
for st in ("fold", "dir"):
    c_ = RC[RC.stat == st]; t_ = c_[c_.kind == "target"]; d_ = c_[(c_.kind == "decoy") & (c_["class"] == "coex")]
    if st == "dir":
        d_ = d_[d_.dirn.to_numpy() == d_.target.map(t_.set_index("target").dirn).to_numpy()]
    Rk, Fk = [], []
    for th in TH:
        r_ = int(((t_.p < th) & (t_.med_delta > 0)).sum())
        e_ = d_.assign(k=(d_.p < th) & (d_.med_delta > 0)).groupby("target").k.mean().sum()
        Rk.append(r_); Fk.append(100 * e_ / r_ if r_ else np.nan)
    curves[st] = (np.array(Rk), np.array(Fk))
pd.DataFrame({"threshold": TH, "called_fold": curves["fold"][0], "fdr_fold": curves["fold"][1], "called_dir": curves["dir"][0],
              "fdr_dir": curves["dir"][1]}).to_csv(D + "v13_fig/threshold_curve.csv", index=False)
axr = C.ax(14, 175, 30, 11); C.letter(6.5, 167.8, "i"); C.head(11.5, 168.0, "The threshold")
axf = C.ax(14, 188, 30, 16)
for st in ("fold", "dir"):
    axr.plot(TH, curves[st][0], color=K.STAT[st], lw=1.0); axf.plot(TH, curves[st][1], color=K.STAT[st], lw=1.0)
for a_ in (axr, axf):
    a_.set_xscale("log"); a_.axvline(0.05, color=K.RULE, lw=0.6, zorder=0); a_.set_xlim(0.001, 0.2); a_.minorticks_off()
axr.set_xticks([]); axr.set_ylabel("called", fontsize=K.T_SMALL); axr.set_ylim(0, 40); axr.set_yticks([0, 20, 40])
axr.spines["bottom"].set_visible(False)
axf.set_xticks([0.001, 0.01, 0.05, 0.2]); axf.set_xticklabels(["0.001", "0.01", "0.05", "0.2"])
axf.set_ylim(0, 40); axf.set_yticks([0, 20, 40]); axf.set_ylabel("FDR (%)", fontsize=K.T_SMALL); axf.set_xlabel("P threshold")
K.legend(axr, [(K.SLATE, "folded", "line"), (K.BLUE, "directed", "line")], x=0.04, y=0.98, dy=0.26, size=K.T_SMALL)

# ================================================================== j  QQ
ax = C.panel("j", 60, 175, 27, 29, head="P values against uniform")
c_ = RC[RC.stat == "fold"]
for msk, cc, lab_, s_ in ((c_.kind == "decoy") & (c_["class"] == "rand"), K.GRAY, "random-gene", 1.5), \
                         ((c_.kind == "decoy") & (c_["class"] == "coex"), K.AQUA, "co-expressed", 1.5), \
                         ((c_.kind == "target"), K.INK, "signatures", 5):
    pv = np.sort(c_[msk].p.to_numpy()); n_ = len(pv); ex = (np.arange(1, n_ + 1) - 0.5) / n_
    ax.scatter(-np.log10(ex), -np.log10(pv), s=s_, color=cc, lw=0, alpha=0.8 if cc != K.INK else 1, zorder=3 if cc == K.INK else 2,
               rasterized=n_ > 2000)
K.identity(ax, 0, 4.5)
ax.set_xlim(0, 4.5); ax.set_ylim(0, 4.5); ax.set_xticks([0, 2, 4]); ax.set_yticks([0, 2, 4])
ax.set_xlabel("expected −log10 P"); ax.set_ylabel("observed −log10 P")
K.legend(ax, [(K.INK, "signatures"), (K.AQUA, "co-expressed decoys"), (K.GRAY, "random-gene decoys")], x=0.04, y=0.97, dy=0.09,
         s=7, size=K.T_SMALL)

# ================================================================== k  call rate by decoy co-expression
RCV = pd.read_csv(D + "calib_v3_rcurve.csv")
ax = C.panel("k", 102, 175, 27, 29, head="Decoy calls by co-expression")
bins_ = list(dict.fromkeys(RCV.rbin))
NC = {"N0": K.tint(K.ORANGE, .35), "N1nat": K.AQUA, "N2nat": K.INK}
for nl in ("N0", "N1nat", "N2nat"):
    g = RCV[(RCV["null"] == nl) & (RCV.stat == "fold")].set_index("rbin").loc[bins_]
    ax.plot(range(len(bins_)), g.call_rate, color=NC[nl], lw=0.9, zorder=2)
    K.dots(ax, range(len(bins_)), g.call_rate, NC[nl], s=7, z=3)
K.ref(ax, 2.5)
ax.text(5.3, 2.2, "nominal 2.5%", ha="right", va="top", fontsize=K.T_SMALL, color=K.INK2)
ax.set_xticks(range(len(bins_))); ax.set_xticklabels([b.replace("(", "").replace("]", "").replace(", ", "–") for b in bins_],
                                                     rotation=40, ha="right", fontsize=K.T_SMALL)
ax.set_xlabel("decoy internal r"); ax.set_ylabel("decoys called (%)"); ax.set_ylim(0, 13); ax.set_yticks([0, 4, 8, 12])
K.legend(ax, [(NC["N0"], "N0", "line"), (NC["N1nat"], "N1, natural pool", "line"), (NC["N2nat"], "N2, natural pool", "line")],
         x=0.04, y=0.97, dy=0.1, size=K.T_SMALL)
ndec = RCV[(RCV["null"] == "N0") & (RCV.stat == "fold")].decoys.sum()

# ================================================================== l  convergence over decoy sets
ax = C.panel("l", 150, 175, 28, 29, head="Pooling decoy sets")
conv = {}
for st in ("fold", "dir"):
    c_ = RC[RC.stat == st]; t_ = c_[c_.kind == "target"]; d_ = c_[(c_.kind == "decoy") & (c_["class"] == "coex")]
    if st == "dir":
        d_ = d_[d_.dirn.to_numpy() == d_.target.map(t_.set_index("target").dirn).to_numpy()]
    r_ = int(t_.call.sum())                     # the estimate sums the decoy call rate over all 52 targets
    conv[st] = [100 * d_[d_.set <= k].groupby("target").call.mean().sum() / r_ for k in range(1, 11)]
    ax.plot(range(1, 11), conv[st], color=K.STAT[st], lw=1.0, zorder=3); K.dots(ax, range(1, 11), conv[st], K.STAT[st], s=6, z=4)
pd.DataFrame({"sets_pooled": range(1, 11), "fdr_fold": conv["fold"], "fdr_dir": conv["dir"]}).to_csv(D + "v13_fig/convergence_sets.csv", index=False)
ax.set_xticks([1, 5, 10]); ax.set_xlim(0.5, 10.5); ax.set_ylim(0, 25); ax.set_yticks([0, 10, 20])
ax.set_xlabel("decoy sets pooled"); ax.set_ylabel("estimated FDR (%)")
K.legend(ax, [(K.SLATE, "folded", "line"), (K.BLUE, "directed", "line")], x=0.45, y=0.97, dy=0.1, size=K.T_SMALL)

C.save("Figure4_GB10")
print("Figure4_GB10 | b", {f"{n}-{s}": (int(fa.loc[(n, s)].called), round(fa.loc[(n, s)].fdr, 1)) for n, s in COMBO},
      "| e sets fold", RB[RB.stat == "fold"].fdr.round(1).tolist(), "| pi0 fold", round(pi0, 3),
      "| q>0.5 fold", int((qs["fold"] > .5).sum()), "| i at 0.05", {st: (int(np.interp(0.05, TH, curves[st][0])),
      round(float(np.interp(0.05, TH, curves[st][1])), 1)) for st in curves}, "| k decoys", int(ndec), "| l", {k: [round(v, 1) for v in vv] for k, vv in conv.items()})
