"""Figure 5 (GB v10) - what the estimate estimates, and how far it sits from it.

a  where every signature is null (1,000 within-cohort label permutations): realized false calls per permutation and the
   plain and pi0-corrected estimates, folded and directed
b  per signature under the same permutations: its own call rate against its decoys'
c  per unit: the call rate of signatures and of each decoy class against the nominal 2.5%
d  planted signal in the response cohorts, truth by gene content (>= 10% of a signature's genes causal)
e  the same replicates against the formal truth (called more often than its decoys over 200 further label draws)
f  how many signatures each truth counts as true
g  planted survival signal in 20 TCGA cancer types with new causal gene sets, formal truth
h  where the plain estimate's excess comes from
i  how often the intervals contain the setting's realized FDR
j  replicate by replicate: realized against estimated FDR (formal truth; folded, beta = 1, 12 cohorts; log axes, realized 0
   drawn at the bottom), with marginal densities
k  the estimate minus the realized FDR in every replicate, by scenario and statistic
l  the pi0 estimate in each of the 1,000 label permutations
Ratios of totals over replicates (200 per response setting, 100 per TCGA setting). validation/PLAN_round7_estimand_and_bias.md
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import glob, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
D = K.D + "/icb/"; T = D + "v12_truth/"
P0 = pd.read_csv(T + "null_pi0.csv"); SN = pd.read_csv(T + "summary_null.csv").set_index("stat")
NS = pd.read_csv(D + "v10_calib_null_summary.csv").set_index(["stat", "estimator"])
PT = pd.read_csv(D + "v10_calib_null_per_target.csv", index_col=0)
SG = pd.read_csv(D + "v10_calib_signal.csv"); SG = SG[SG.truth == "10%"]
SR = pd.read_csv(T + "summary_response.csv"); ST = pd.read_csv(T + "summary_tcga.csv")
SIG = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(T + "signal_*.csv"))])
SCOL = {"P1 inflamed": K.BLUE, "P2 non-inflamed": K.ORANGE, "P3 outside signatures": K.AQUA}
SLAB = {"P1 inflamed": "3 inflamed signatures", "P2 non-inflamed": "3 other signatures", "P3 outside signatures": "outside every signature"}
TCOL = {"T1 library": K.BLUE, "T2 outside": K.AQUA, "T3 mixed": K.ORANGE}
TLAB = {"T1 library": "2 signatures", "T2 outside": "1 Hallmark set (out of library)", "T3 mixed": "both"}

C = K.Canvas(205)

# ================================================================== a  global null
for gi, st in enumerate(("fold", "dir")):
    x = P0[P0.stat == st]
    ax = C.panel("a", 14, 12, 30, 28, head="All signatures null") if gi == 0 else C.ax(52, 12, 30, 28)
    mx = int(max(x.R.max(), 8))
    bins = np.arange(-0.5, mx + 1.5, 1)
    K.hist_steps(ax, x.R, bins, K.SLATE, alpha=0.35)
    for v, filled, lab in ((x.E.mean(), True, "plain"), (x.E0.mean(), False, "π0")):
        ax.axvline(v, color=K.AQUA, lw=1.0 if filled else 0.8, ls="-" if filled else (0, (1.2, 1.2)), zorder=4)
    ax.axvline(x.R.mean(), color=K.INK, lw=1.0, zorder=3)
    ax.set_xlim(-0.6, 8.6); ax.set_xticks(range(0, 9, 2))
    ax.set_xlabel("false calls per library")
    if gi == 0:
        ax.set_ylabel("permutations (share)")
    r = SN.loc[st]
    ax.text(0.97, 0.97, f"{K.STATLAB[st]}\nplain {r.ratio:.2f}\nπ0 {r.ratio_pi0:.2f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=K.T_SMALL, linespacing=1.15)
K.hlegend(C, 14, 50.5, [(K.INK, "realized (mean)", "line"), (K.AQUA, "estimate", "line")], size=K.T_SMALL, gap=2)
K.hlegend(C, 58, 50.5, [(K.AQUA, "π0-corrected", "line")], size=K.T_SMALL, gap=2)

# ================================================================== b  per signature under permutation
ax = C.panel("b", 98, 12, 28, 28, head="Per signature")
for st, mk in (("fold", "o"), ("dir", "^")):
    K.dots(ax, 100 * PT[f"{st}_decoy_callrate"], 100 * PT[f"{st}_target_callrate"], K.STAT[st], s=7, marker=mk)
K.identity(ax, 0, 9)
ax.set_xlim(0, 8); ax.set_ylim(0, 8); ax.set_xticks([0, 4, 8]); ax.set_yticks([0, 4, 8])
ax.set_xlabel("decoys called (%)"); ax.set_ylabel("signature called (%)")
rf, rd = NS.loc[("fold", "coex A"), "per_target_r"], NS.loc[("dir", "coex A"), "per_target_r"]
K.legend(ax, [(K.SLATE, f"folded, r = {rf:.2f}"), (K.BLUE, f"directed, r = {rd:.2f}", "^")], x=0.06, y=0.95, dy=0.1, s=7,
         size=K.T_SMALL)

# ================================================================== c  per unit
ax = C.panel("c", 142, 12, 36, 28, head="Per unit")
groups = [("sig.", None), ("co-expr.", "coex A"), ("no-TME", "coex B"), ("random", "random A")]
for gi, (lab, est) in enumerate(groups):
    for si, st in enumerate(("fold", "dir")):
        r = NS.loc[(st, "coex A" if est is None else est)]
        v = r.target_test_rate if est is None else r.decoy_test_rate
        ax.bar(gi + (si - .5) * .36, v, .34, color=K.STAT[st], edgecolor="white", linewidth=0.5, zorder=3)
K.ref(ax, 2.5, label="2.5%", pos=0.995)
ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("called (%)"); ax.set_ylim(0, 4.4); ax.set_xlim(-.6, 4.1)
K.legend(ax, [(K.SLATE, "folded", "box"), (K.BLUE, "directed", "box")], x=0.70, y=0.97, dy=0.1, s=8, size=K.T_SMALL)


# ================================================================== d, e, g  planted signal, log axes
def logax(ax, lo=0.05, hi=300):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo * 20, hi)
    tk = [0.1, 1, 10, 100]
    ax.set_xticks(tk); ax.set_xticklabels(["0.1", "1", "10", "100"])
    ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"])
    ax.minorticks_off()
    K.identity(ax, lo, hi)
    ax.set_xlabel("realized FDR (%)"); ax.set_ylabel("estimate (%)")


ax = C.panel("d", 14, 60, 34, 32, head="Truth by gene content")
for (sc, b, c, st), g in SG.groupby(["scenario", "beta", "cohorts", "stat"]):
    xr, ye = 100 * g.V.sum() / g.R.sum(), 100 * g["E coex A"].sum() / g["R coex A"].sum()
    K.dots(ax, [xr], [ye], SCOL[sc], s=13 if c == 12 else 7, marker="o" if st == "fold" else "^", alpha=1 if b == 1 else 0.55)
logax(ax)
ax.text(0.04, 0.96, "above the line:\nestimate high", transform=ax.transAxes, fontsize=K.T_SMALL - 0.3, color=K.INK2, va="top")

ax = C.panel("e", 62, 60, 34, 32, head="Formal truth")
for r in SR.itertuples():
    for v, filled in ((r.estimate, True), (r.estimate_pi0, False)):
        K.dots(ax, [max(r.realized, 0.06)], [v], SCOL[r.scenario], s=13 if r.cohorts == 12 else 7, filled=filled,
               marker="o" if r.stat == "fold" else "^", alpha=1 if r.beta == 1 else 0.6)
logax(ax)

ax = C.panel("g", 144, 60, 34, 32, head="TCGA survival, new causal sets")
for r in ST.itertuples():
    for v, filled in ((r.estimate, True), (r.estimate_pi0, False)):
        K.dots(ax, [r.realized], [v], TCOL[r.kind], s=14 if r.beta == 0.6 else 8, filled=filled, marker="s")
logax(ax)
ax.text(0.04, 0.96, "TCGA decoys, 0.34 vs 0.47\nco-expression", transform=ax.transAxes, fontsize=K.T_SMALL - 0.3,
        color=K.INK2, va="top", linespacing=1.05)
C.letter(110.5, 60.8, "g") if False else None

# legends for d, e, g
K.hlegend(C, 14, 101.5, [(SCOL[k], f"response: {SLAB[k]}") for k in SCOL], size=K.T_SMALL, s=8, gap=2.0)
K.hlegend(C, 14, 105.3, [(TCOL[k], f"TCGA: {TLAB[k]}", "s") for k in TCOL] + [(K.INK, "plain", "o", True),
                                                                               (K.INK, "π0-corrected", "o", False),
                                                                               (K.INK, "directed", "^", True)],
          size=K.T_SMALL, s=8, gap=2.0)

# ================================================================== f  how many are true
ax = C.panel("f", 112, 60, 17, 32, head="Signatures true")
nt = SIG.groupby("scenario")[["n_true_gc", "n_true_cp"]].mean()
nt.loc["TCGA, all"] = [np.nan, ST.n_true_cp.mean()]
yy = 0
labs = []
for sc in ["P1 inflamed", "P2 non-inflamed", "P3 outside signatures"]:
    ax.barh(yy - 0.18, nt.loc[sc, "n_true_gc"], 0.34, color=K.LGRAY, edgecolor="white", linewidth=0.4)
    ax.barh(yy + 0.18, nt.loc[sc, "n_true_cp"], 0.34, color=SCOL[sc], edgecolor="white", linewidth=0.4)
    labs.append(sc.split()[0]); yy += 1
ax.set_yticks(range(3)); ax.set_yticklabels(labs); ax.set_ylim(2.6, -0.6); K.frame(ax, left=False)
ax.set_xlim(0, 52); ax.set_xticks([0, 25, 50]); ax.set_xlabel("signatures (mean)")
K.legend(ax, [(K.LGRAY, "gene content", "box"), (K.INK2, "formal", "box")], x=0.50, y=0.16, dy=0.10, s=7, size=K.T_SMALL - 0.4)

# ================================================================== h  decomposition
ax = C.panel("h", 22, 118, 92, 34, head="Where the excess comes from")
rows = []
for r in SR[(SR.beta == 1.0) & (SR.cohorts == 12)].itertuples():
    rows.append((f"{r.scenario.split()[0]} {K.STATLAB[r.stat]}", r, SCOL[r.scenario]))
for r in ST[ST.beta == .6].itertuples():
    rows.append((f"TCGA {r.kind.split()[0]}", r, TCOL[r.kind]))
for i, (lab, r, cc) in enumerate(rows):
    parts = [(r.exp_V, K.INK), (r.excess_null, cc), (r.part_true, K.tint(cc, 0.6))]
    left = 0
    for v, col in parts:
        ax.barh(i, v, 0.62, left=left, color=col, edgecolor="white", linewidth=0.5, zorder=3)
        left += v
    K.dots(ax, [r.estimate_pi0], [i], K.INK, s=11, filled=False, marker="D", z=6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([x[0] for x in rows]); ax.set_ylim(len(rows) - 0.4, -0.6)
K.frame(ax, left=False); ax.set_xlabel("% of calls (expected)"); ax.set_xlim(0, 240)
K.legend(ax, [(K.INK, "false calls of null signatures", "box"), (K.GRAY, "decoys' excess over null signatures", "box"),
              (K.tint(K.GRAY, 0.6), "decoy rates of true signatures", "box"), (K.INK, "π0-corrected estimate", "D", False)],
         x=0.62, y=0.32, dy=0.085, s=8, size=K.T_SMALL)

# ================================================================== i  coverage
ax = C.panel("i", 140, 118, 34, 34, head="Interval coverage (%)")
groups = [("P1", SR[SR.scenario.str.startswith("P1")]), ("P2", SR[SR.scenario.str.startswith("P2")]),
          ("P3", SR[SR.scenario.str.startswith("P3")]), ("TCGA T1", ST[ST.kind.str.startswith("T1")]),
          ("TCGA T2", ST[ST.kind.str.startswith("T2")]), ("TCGA T3", ST[ST.kind.str.startswith("T3")])]
cov = np.array([[g.cov_t_fdr.mean(), g.cov_c_fdr.mean(), g.cov0_t_fdr.mean(), g.cov0_c_fdr.mean()] for _, g in groups]) * 100
ax.imshow(cov, aspect="auto", cmap=K.SEQ_BLUE, vmin=0, vmax=100)
for i in range(cov.shape[0]):
    for j in range(cov.shape[1]):
        ax.text(j, i, f"{cov[i, j]:.0f}", ha="center", va="center", fontsize=K.T_SMALL, color="white" if cov[i, j] > 55 else K.INK)
ax.set_yticks(range(len(groups))); ax.set_yticklabels([g[0] for g in groups])
ax.set_xticks(range(4)); ax.set_xticklabels(["sig.", "clust.", "sig.", "clust."]); ax.tick_params(length=0)
[sp.set_visible(False) for sp in ax.spines.values()]
ax.text(0.5, -0.62, "plain", ha="center", va="bottom", fontsize=K.T_SMALL)
ax.text(2.5, -0.62, "π0-corrected", ha="center", va="bottom", fontsize=K.T_SMALL)
ax.axvline(1.5, color="white", lw=1.5)
# ================================================================== j  replicate level
ax = C.ax(14, 174, 34, 25); C.letter(6.5, 166.0, "j"); C.head(11.5, 166.2, "Replicate by replicate")
sub = SIG[(SIG.stat == "fold") & (SIG.beta == 1.0) & (SIG.cohorts == 12) & (SIG.R > 0)].copy()
sub["est"] = 100 * sub.E / sub.R; sub["real"] = 100 * sub.V_cp / sub.R
FL, TOP = 0.5, 500.0                                               # log axes; a realized FDR of 0 is drawn at 0.5%
lg = lambda v: np.log10(np.clip(v, FL, TOP))
for sc in SCOL:
    g = sub[sub.scenario == sc]
    ax.scatter(lg(g.est), lg(g.real) + np.where(g.real <= 0, np.random.default_rng(2).uniform(-0.06, 0.06, len(g)), 0), s=2.2,
               color=SCOL[sc], alpha=0.45, lw=0, zorder=3)
K.identity(ax, lg(FL), lg(TOP))
tk = np.array([1, 10, 100])
ax.set_xlim(lg(FL) - 0.1, lg(TOP)); ax.set_ylim(lg(FL) - 0.12, lg(TOP))
ax.set_xticks(np.log10(tk)); ax.set_xticklabels(["1", "10", "100"])
ax.set_yticks(np.r_[lg(FL), np.log10(tk)]); ax.set_yticklabels(["0", "1", "10", "100"])
ax.axhline(lg(FL) + 0.15, color=K.GRID, lw=0.4, zorder=1)
ax.set_xlabel("estimated FDR (%)"); ax.set_ylabel("realized FDR (%)")
K.marginal(C, 14, 174, 34, 25, [lg(sub[sub.scenario == sc].est) for sc in SCOL], [lg(sub[sub.scenario == sc].real) for sc in SCOL],
           [SCOL[sc] for sc in SCOL], size=3.5, xlim=(lg(FL) - 0.1, lg(TOP)), ylim=(lg(FL) - 0.12, lg(TOP)))
above = 100 * (sub.est >= sub.real).mean()
ax.text(0.04, 0.97, f"estimate ≥ realized in\n{above:.0f}% of {len(sub)} replicates", transform=ax.transAxes, ha="left",
        va="top", fontsize=K.T_SMALL, color=K.INK2, linespacing=1.1)

# ================================================================== k  estimate minus realized
ax = C.panel("k", 72, 174, 56, 25, head="Estimate minus realized, every replicate")
allr = SIG[SIG.R > 0].copy(); allr["gap"] = 100 * (allr.E - allr.V_cp) / allr.R
CAP = 250
pos_ = 0; ticks, tlabs = [], []
for sc in SCOL:
    for st in ("fold", "dir"):
        v = allr[(allr.scenario == sc) & (allr.stat == st)].gap
        K.violin(ax, np.clip(v, None, CAP), pos_, SCOL[sc], width=0.4, alpha=0.85 if st == "fold" else 0.4)
        ax.text(pos_, CAP + 8, f"{100 * (v >= 0).mean():.0f}%", ha="center", va="bottom", fontsize=K.T_SMALL, color=K.INK2)
        ticks.append(pos_); tlabs.append(f"{sc.split()[0]} {'F' if st == 'fold' else 'D'}"); pos_ += 1
    pos_ += 0.4
K.ref(ax, 0)
ax.set_xticks(ticks); ax.set_xticklabels(tlabs, fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylabel("percentage points"); ax.set_ylim(-60, CAP + 32); ax.set_yticks([0, 100, 200])
ax.set_xlim(-0.6, pos_ - 0.8)

# ================================================================== l  pi0 under the global null
ax = C.panel("l", 148, 174, 30, 25, head="π0 when all are null")
bins = np.linspace(0, 1, 27)
for st in ("fold", "dir"):
    v = P0[P0.stat == st].pi0
    if st == "fold":
        K.hist_steps(ax, v, bins, K.STAT[st], alpha=0.5)
    else:
        ax.hist(v, bins=bins, histtype="step", color=K.STAT[st], lw=0.9, density=True, zorder=3)
    ax.axvline(v.mean(), color=K.STAT[st], lw=0.9, ls=(0, (2, 1.2)), zorder=4)
ax.set_xlim(0, 1.02); ax.set_xticks([0, 0.5, 1]); ax.set_xlabel("π0 estimate"); ax.set_ylabel("density")
K.legend(ax, [(K.SLATE, f"folded, mean {P0[P0.stat == 'fold'].pi0.mean():.2f}", "box"),
              (K.BLUE, f"directed, mean {P0[P0.stat == 'dir'].pi0.mean():.2f}", "line")], x=0.03, y=0.97, dy=0.1, s=7, size=K.T_SMALL)

C.save("Figure5_GB10")
print("Figure5_GB10 | null", SN[["ratio", "ratio_pi0"]].round(3).to_dict(), "| n true", nt.round(1).to_dict(), "| j above",
      round(above, 1), len(sub), "| k", allr.groupby(["scenario", "stat"]).gap.apply(lambda v: round(100 * (v >= 0).mean(), 1)).to_dict(),
      "| l pi0 mean", P0.groupby("stat").pi0.mean().round(3).to_dict())
