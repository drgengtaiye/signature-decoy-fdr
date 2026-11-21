"""Supplementary Figure 4: calibration of the target-decoy estimate in detail.   python run_restyled.py figS4.py

a  per signature: how often it is called under permutation against the mean call rate of its decoys
b  planted signal: the realized false discovery proportion under three gene-content definitions of a true signature,
   against the decoy estimate (beta = 1, 12 cohorts)
c  the same with decoys matched less or more closely (random-gene, co-expressed, closest-matched)
d  intervals for the observed FDR: the five resampling schemes of the plan, and subsampling 9 of 12 cohorts
e  coverage of the realized FDP by the target and cluster bootstrap intervals in the planted-signal replicates, truth by
   gene content
f  coverage of the setting's realized FDR under the call-probability truth, plain and pi0-corrected estimates, response
   and TCGA planted signal (validation/PLAN_round7_estimand_and_bias.md, analysis 6)
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v4kit as K
import figstyle as fs
K.apply()
D = K.D + "/icb/"
TT = 6.9
C_FOLD, C_DIR, C_SIG, C_DEC = fs.GREY, fs.BLUE, fs.TEAL, fs.SLATE
NS = pd.read_csv(D + "v10_calib_null_summary.csv").set_index(["stat", "estimator"])
PT = pd.read_csv(D + "v10_calib_null_per_target.csv", index_col=0)
SG = pd.read_csv(D + "v10_calib_signal.csv")
IV = pd.read_csv(D + "v10_fdr_intervals.csv")
DG = pd.read_csv(D + "v10_cohort_boot_diag.csv")
C = K.Canvas(160)

# ================================================================= b
ax = C.ax(15, 15, 55, 30)
for st, cc, mk in (("fold", C_FOLD, "o"), ("dir", C_DIR, "^")):
    ax.scatter(100 * PT[f"{st}_decoy_callrate"], 100 * PT[f"{st}_target_callrate"], s=7, marker=mk, color=cc, lw=0,
               alpha=.85, zorder=3, label=f"{'folded' if st == 'fold' else 'directed'}, r = "
                                         f"{np.corrcoef(PT[f'{st}_decoy_callrate'], PT[f'{st}_target_callrate'])[0, 1]:.2f}")
ax.plot([0, 9], [0, 9], color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_xlim(0, 9); ax.set_ylim(0, 9)
ax.set_xlabel("its decoys' call rate (%)"); ax.set_ylabel("signature called (%)")
ax.legend(fontsize=K.T_SMALL, loc="upper left", handletextpad=.3, borderpad=.2, frameon=False)
K.frame(ax)
C.title(15, 9, "Per signature, under permutation", size=TT, fontweight="bold")
C.letter(7.5, 8.5, "a")

# ================================================================= c
ax = C.ax(95, 15, 70, 30)
SC = [("P1 inflamed", "signal in 3\ninflamed sets"), ("P2 non-inflamed", "signal in 3\nother sets"),
      ("P3 outside signatures", "signal outside\nevery set")]
TRUTH = [("any gene", "o", "white"), ("10%", "o", fs.INK), ("25%", "s", fs.INK)]
for gi, (sc, lab) in enumerate(SC):
    for si, st in enumerate(("fold", "dir")):
        x = gi + (si - .5) * .38
        g0 = SG[(SG.scenario == sc) & (SG.beta == 1.0) & (SG.cohorts == 12) & (SG.stat == st)]
        for tdef, mk, fc in TRUTH:
            g = g0[g0.truth == tdef]
            ax.scatter(x, 100 * g.V.sum() / g.R.sum(), s=12, marker=mk, facecolor=fc, edgecolor=fs.INK, lw=.6, zorder=4)
        g = g0[g0.truth == "10%"]
        ax.scatter(x, min(100 * g["E coex A"].sum() / g["R coex A"].sum(), 245), s=20, marker="D",
                   color=C_FOLD if st == "fold" else C_DIR, zorder=5, clip_on=False)
ax.set_xticks(range(3)); ax.set_xticklabels([s[1] for s in SC], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylim(0, 250); ax.set_yticks([0, 50, 100, 150, 200, 250]); ax.set_yticklabels(["0", "50", "100", "150", "200", "≥250"])
ax.set_ylabel("false discoveries (%)")
for k, (mk, fc, ec, l) in enumerate([("D", C_FOLD, C_FOLD, "decoy estimate (grey folded, blue directed)"),
                                     ("o", "white", fs.INK, "realized; true if it shares any causal gene"),
                                     ("o", fs.INK, fs.INK, "realized; true if ≥ 10% of its genes are causal"),
                                     ("s", fs.INK, fs.INK, "realized; true if ≥ 25%")]):
    ax.scatter([.03], [.95 - .085 * k], s=10, marker=mk, facecolor=fc, edgecolor=ec, lw=.6, transform=ax.transAxes, clip_on=False, zorder=10)
    ax.text(.065, .95 - .085 * k, l, transform=ax.transAxes, va="center", fontsize=K.T_SMALL, color=fs.INK)
K.frame(ax)
C.title(95, 9, "What counts as false sets the bias", size=TT, fontweight="bold")
C.letter(87.5, 8.5, "b")

# ================================================================= d
ax = C.ax(15, 66, 55, 30)
EST = [("E random A", "random-gene", C_DEC), ("E coex A", "co-expressed", C_SIG), ("E coex A, closest", "closest-matched", fs.AMBER)]
for gi, (sc, lab) in enumerate(SC):
    g = SG[(SG.scenario == sc) & (SG.beta == 1.0) & (SG.cohorts == 12) & (SG.stat == "dir") & (SG.truth == "10%")]
    ax.plot([gi - .4, gi + .4], [100 * g.V.sum() / g.R.sum()] * 2, color=fs.INK, lw=1.0, zorder=3)
    for k, (col, l2, cc) in enumerate(EST):
        rk = col.replace("E ", "R ")
        ax.scatter(gi + (k - 1) * .22, min(100 * g[col].sum() / g[rk].sum(), 245), s=14, color=cc, zorder=4, clip_on=False)
ax.set_xticks(range(3)); ax.set_xticklabels([s[1] for s in SC], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_ylim(0, 250); ax.set_yticks([0, 100, 200]); ax.set_ylabel("false discoveries (%), directed")
K.key(ax, [(cc, l2) for _, l2, cc in EST], x=.03, y=.97, s=9, size=K.T_SMALL)
fs.note(ax, .03, .66, "— realized (≥ 10% of genes)", size=K.T_SMALL, color=fs.INK)
K.frame(ax)
C.title(15, 60, "Matching the decoys more closely", size=TT, fontweight="bold")
C.letter(7.5, 59.5, "c")

# ================================================================= e
ax = C.ax(95, 66, 70, 34)
ROWS = ["(a) targets", "(b) clusters, cut 0.4", "(b) clusters, cut 0.5", "(b) clusters, cut 0.6", "(c) decoys only",
        "(d) cohorts", "(e) clusters x cohorts"]
LAB = {"(a) targets": "targets (published)", "(b) clusters, cut 0.4": "clusters, 10", "(b) clusters, cut 0.5": "clusters, 13",
       "(b) clusters, cut 0.6": "clusters, 20", "(c) decoys only": "decoys only", "(d) cohorts": "cohorts",
       "(e) clusters x cohorts": "clusters × cohorts"}
for si, (an, cc, off) in enumerate([("N0 folded", C_FOLD, .17), ("N0 directed", C_DIR, -.17)]):
    q = IV[IV.analysis == an].set_index("interval")
    pt = q.fdr.iloc[0]
    for i, r in enumerate(ROWS):
        y = len(ROWS) - i + off
        lo, hi = q.loc[r, "lo"], q.loc[r, "hi"]
        ax.plot([lo, min(hi, 120)], [y, y], color=cc, lw=1.1, solid_capstyle="butt", zorder=3)
        if hi > 120: ax.text(121, y, f"{hi:.0f}", fontsize=4.6, va="center", color=cc)
    ax.scatter([pt] * len(ROWS), [len(ROWS) - i + off for i in range(len(ROWS))], s=6, color=cc, zorder=4)
    dg = DG[(DG.stat == ("fold" if "folded" in an else "dir")) & (DG.resampling == "9 of 12 without replacement")].iloc[0]
    ax.scatter(dg.fdr_ratio_of_means, 0 + off, s=10, marker="D", color=cc, zorder=4)
ax.set_yticks(list(range(len(ROWS), 0, -1)) + [0])
ax.set_yticklabels([LAB[r] for r in ROWS] + ["9 of 12 cohorts, no repeats"], fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0); ax.set_xlim(0, 130); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("target-decoy FDR, 95% interval (%)")
K.key(ax, [(C_FOLD, "folded, 17%"), (C_DIR, "directed, 7%")], x=.45, y=.985, s=9, size=K.T_SMALL)
K.frame(ax, left=False)
C.title(95, 60, "What the interval resamples changes it", size=TT, fontweight="bold")
C.letter(87.5, 59.5, "d")

# ================================================================= f
ax = C.ax(15, 117, 55, 34)
S1 = SG[(SG.truth == "10%") & (SG.R > 0)]
cv = S1.groupby(["scenario", "beta", "cohorts", "stat"])[["cover_target", "cover_cluster"]].mean().reset_index()
cv = cv.sort_values(["scenario", "cohorts", "beta", "stat"])
for k, (_, r) in enumerate(cv.iterrows()):
    ax.scatter(100 * r.cover_target, k, s=9, color=C_DEC, zorder=4)
    ax.scatter(100 * r.cover_cluster, k, s=9, color=C_SIG, zorder=4)
    ax.plot([100 * r.cover_target, 100 * r.cover_cluster], [k, k], color=fs.LIGHT, lw=.8, zorder=2)
ax.axvline(95, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_yticks(range(len(cv)))
ax.set_yticklabels([f"{r.scenario.split()[0]} β{r.beta:g} {r.cohorts}c {('fold' if r.stat == 'fold' else 'dir')}" for _, r in cv.iterrows()],
                   fontsize=4.4)
ax.tick_params(axis="y", length=0); ax.set_xlim(0, 100)
ax.set_xlabel("interval covers the realized FDP (%)\ndark: target bootstrap; teal: cluster bootstrap", fontsize=K.T_SMALL)
K.frame(ax, left=False)
C.title(15, 111, "Coverage in the planted replicates", size=TT, fontweight="bold")
C.letter(7.5, 110.5, "e")

# ================================================================= f  coverage, call-probability truth
SR = pd.read_csv(D + "v12_truth/summary_response.csv"); STC = pd.read_csv(D + "v12_truth/summary_tcga.csv")
ax = C.ax(105, 117, 60, 34)
rows = [(f"{r.scenario.split()[0]} β{r.beta:g} {r.cohorts}c {('fold' if r.stat == 'fold' else 'dir')}", r)
        for r in SR.sort_values(["scenario", "cohorts", "beta", "stat"]).itertuples()]
rows += [(f"TCGA {r.kind.split()[0]} β{r.beta:g} fold", r) for r in STC.itertuples()]
for k, (lab, r) in enumerate(rows):
    for v, cc, mk in ((r.cov_t_fdr, C_DEC, "o"), (r.cov_c_fdr, C_SIG, "o"), (r.cov0_t_fdr, C_DEC, "D"), (r.cov0_c_fdr, C_SIG, "D")):
        ax.scatter(100 * v, k, s=8, marker=mk, facecolor=cc if mk == "o" else "white", edgecolor=cc, lw=.6, zorder=4)
ax.axvline(95, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_yticks(range(len(rows))); ax.set_yticklabels([x[0] for x in rows], fontsize=4.0)
ax.tick_params(axis="y", length=0); ax.set_xlim(-2, 102)
ax.set_xlabel("interval covers the realized FDR (%)\nfilled: plain; open: π0-corrected; dark target, teal cluster", fontsize=K.T_SMALL)
K.frame(ax, left=False)
C.title(105, 111, "Coverage, truth by call probability", size=TT, fontweight="bold")
C.letter(97.5, 110.5, "f")

C.save("SupplementaryFigure4")
print("FigS_calib |", NS.round(3).to_string())
