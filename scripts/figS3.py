"""Additional file 1, Fig. S2 (GB v11): the panels of the v10 Figure 5 moved out of the main figure, unchanged in content
(validation/PLAN_round9_standard_fdr_comparison.md, B), in the order the text cites them.

a  per unit: the call rate of signatures and of each decoy class under the global null against the nominal 2.5% (v10 5c)
b  the pi0 estimate in each of the 1,000 label permutations (v10 5l)
c  how many signatures each truth counts as true (v10 5f)
d  how often the intervals contain the setting's realized FDR (v10 5i)
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import glob, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
D = K.D + "/icb/"; T = D + "v12_truth/"
P0 = pd.read_csv(T + "null_pi0.csv")
NS = pd.read_csv(D + "v10_calib_null_summary.csv").set_index(["stat", "estimator"])
SR = pd.read_csv(T + "summary_response.csv"); ST = pd.read_csv(T + "summary_tcga.csv")
SIG = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(T + "signal_*.csv"))])
SCOL = {"P1 inflamed": K.BLUE, "P2 non-inflamed": K.ORANGE, "P3 outside signatures": K.AQUA}

C = K.Canvas(52)

# ================================================================== a  per unit
ax = C.panel("a", 14, 12, 34, 28, head="Per unit, all null")
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

# ================================================================== b  how many are true
ax = C.panel("c", 108, 12, 22, 28, head="Signatures true")
nt = SIG.groupby("scenario")[["n_true_gc", "n_true_cp"]].mean()
for yy, sc in enumerate(["P1 inflamed", "P2 non-inflamed", "P3 outside signatures"]):
    ax.barh(yy - 0.18, nt.loc[sc, "n_true_gc"], 0.34, color=K.LGRAY, edgecolor="white", linewidth=0.4)
    ax.barh(yy + 0.18, nt.loc[sc, "n_true_cp"], 0.34, color=SCOL[sc], edgecolor="white", linewidth=0.4)
ax.set_yticks(range(3)); ax.set_yticklabels(["P1", "P2", "P3"]); ax.set_ylim(2.6, -0.6); K.frame(ax, left=False)
ax.set_xlim(0, 52); ax.set_xticks([0, 25, 50]); ax.set_xlabel("signatures (mean)")
K.legend(ax, [(K.LGRAY, "gene content", "box"), (K.SLATE, "formal (scenario hue)", "box")], x=0.38, y=0.16, dy=0.10, s=7,
         size=K.T_SMALL - 0.4)

# ================================================================== c  coverage
ax = C.panel("d", 150, 12, 28, 28, head="Interval coverage (%)")
grp = [("P1", SR[SR.scenario.str.startswith("P1")]), ("P2", SR[SR.scenario.str.startswith("P2")]),
       ("P3", SR[SR.scenario.str.startswith("P3")]), ("TCGA T1", ST[ST.kind.str.startswith("T1")]),
       ("TCGA T2", ST[ST.kind.str.startswith("T2")]), ("TCGA T3", ST[ST.kind.str.startswith("T3")])]
cov = np.array([[g.cov_t_fdr.mean(), g.cov_c_fdr.mean(), g.cov0_t_fdr.mean(), g.cov0_c_fdr.mean()] for _, g in grp]) * 100
ax.imshow(cov, aspect="auto", cmap=K.SEQ_BLUE, vmin=0, vmax=100)
for i in range(cov.shape[0]):
    for j in range(cov.shape[1]):
        ax.text(j, i, f"{cov[i, j]:.0f}", ha="center", va="center", fontsize=K.T_SMALL, color="white" if cov[i, j] > 55 else K.INK)
ax.set_yticks(range(len(grp))); ax.set_yticklabels([g[0] for g in grp])
ax.set_xticks(range(4)); ax.set_xticklabels(["sig.", "clust.", "sig.", "clust."]); ax.tick_params(length=0)
[sp.set_visible(False) for sp in ax.spines.values()]
ax.text(0.5, -0.62, "plain", ha="center", va="bottom", fontsize=K.T_SMALL)
ax.text(2.5, -0.62, "π0-corrected", ha="center", va="bottom", fontsize=K.T_SMALL)
ax.axvline(1.5, color="white", lw=1.5)

# ================================================================== d  pi0 under the global null
ax = C.panel("b", 64, 12, 28, 28, head="π0 when all are null")
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
              (K.BLUE, f"directed, mean {P0[P0.stat == 'dir'].pi0.mean():.2f}", "line")], x=0.03, y=0.97, dy=0.1, s=7,
         size=K.T_SMALL)

C.save("FigS2_calib_extra_GB11")
print("FigS2_calib_extra_GB11 | n true", nt.round(1).to_dict(), "| pi0 mean", P0.groupby("stat").pi0.mean().round(3).to_dict())
