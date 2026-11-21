"""Supplementary Figure 5: the package rerun and the variability of one draw of decoys (validation V3 of
validation/PLAN_round6_applications.md; analysis 4 of PLAN_round7_estimand_and_bias.md): ten decoy sets, seed 20261022.
python run_restyled.py figS5.py

a  each signature's median delta over the 12 cohorts, published pipeline (calib_v3, N0) against the package rerun,
   folded and directed
b  the estimate from each of the ten decoy sets (points), pooled over the ten (diamond), and from the first decoy set of
   the response analysis (ring), per application
c  decoys' internal correlation against their target's, rerun decoys (median over cohorts)
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
U = pd.read_csv(D + "calib_v3_units.csv"); U = U[(U["null"] == "N0") & (U.kind == "target")]
R = pd.read_csv(D + "v12_pkg/response_calls.csv.gz"); R = R[R.kind == "target"]
F0 = pd.read_csv(D + "calib_v3_fdr.csv"); F0 = F0[F0["null"] == "N0"]
P = D + "v12_pkg/"
I1 = pd.read_csv(P + "response_decoy_info.csv.gz"); I1 = I1[I1.set == 1]
C = K.Canvas(62)

ax = C.ax(15, 13, 40, 40)
out = {}
for st, cc, mk in (("fold", C_FOLD, "o"), ("dir", C_DIR, "^")):
    a = U[U.stat == st].set_index("unit").med_delta; b = R[R.stat == st].set_index("unit").med_delta.reindex(a.index)
    ax.scatter(a, b, s=8, marker=mk, color=cc, lw=0, zorder=3, label=f"{'folded' if st == 'fold' else 'directed'}, r = {np.corrcoef(a, b)[0, 1]:.3f}")
    ca = U[U.stat == st].set_index("unit").call.astype(str).eq("True"); cb = R[R.stat == st].set_index("unit").call.reindex(ca.index).astype(bool)
    out[st] = (int(ca.sum()), int(cb.sum()), int((ca & cb).sum()), float(np.corrcoef(a, b)[0, 1]))
lo, hi = -.17, .14
ax.plot([lo, hi], [lo, hi], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.axhline(0, color=fs.LIGHT, lw=.6); ax.axvline(0, color=fs.LIGHT, lw=.6)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("median Δ, published pipeline"); ax.set_ylabel("median Δ, package rerun")
ax.legend(fontsize=K.T_SMALL, loc="upper left", frameon=False, handletextpad=.3)
K.frame(ax)
C.title(15, 7, "Each signature, published and rerun", size=TT, fontweight="bold"); C.letter(7.5, 6.5, "a")

ax = C.ax(75, 13, 50, 40)
apps = [("response\nfolded", "response", "fold"), ("response\ndirected", "response", "dir"), ("TCGA OS\nfolded", "tcga_OS", "fold"),
        ("breast PC1\nfolded", "breast_pc1", "fold"), ("breast\nmean z", "breast_mean", "fold")]
rows = []
for k, (lab, app, st) in enumerate(apps):
    bs = pd.read_csv(P + f"{app}_fdr_by_set.csv"); bs = bs[bs.stat == st]
    pool = pd.read_csv(P + f"{app}_fdr.csv"); pool = pool[(pool.stat == st) & (pool.decoys == "coex") & (pool.boot == "targets")].iloc[0]
    cc = C_DIR if st == "dir" else C_FOLD
    ax.scatter(k + np.linspace(-.18, .18, len(bs)), bs.fdr, s=6, color=cc, lw=0, zorder=3)
    ax.scatter(k, pool.fdr, s=22, marker="D", color=cc, edgecolor=fs.INK, lw=.5, zorder=4)
    if app == "response":
        f0 = F0[(F0.stat == st) & (F0.pool == "A") & (F0.decoys == "coex")].iloc[0]
        ax.scatter(k, f0.fdr, s=26, facecolor="none", edgecolor=fs.INK, lw=.8, zorder=5)
    rows.append((lab.replace("\n", " "), round(pool.fdr, 1), round(bs.fdr.min(), 1), round(bs.fdr.max(), 1)))
ax.set_xticks(range(len(apps))); ax.set_xticklabels([a[0] for a in apps], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_ylabel("estimated FDR (%)"); ax.set_xlim(-.6, len(apps) - .4); ax.set_ylim(0, 80)
fs.note(ax, .02, .97, "small points, each decoy set; diamonds, pooled over the ten;\nring, the first decoy set of the response analysis", size=K.T_SMALL, color=fs.SLATE)
K.frame(ax)
C.title(75, 7, "One draw of decoys against ten", size=TT, fontweight="bold"); C.letter(67.5, 6.5, "b")

ax = C.ax(140, 13, 32, 40)
for cls, cc in (("rand", C_DEC), ("coex", C_SIG)):
    g = I1[I1["class"] == cls]
    ax.scatter(g.r_target_med, g.r_achieved_med, s=3, color=cc, alpha=.5, lw=0, zorder=3)
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(0, 1); ax.set_ylim(-.05, 1); ax.set_xlabel("target internal $r$"); ax.set_ylabel("decoy internal $r$")
K.key(ax, [(C_SIG, "co-expressed"), (C_DEC, "random-gene")], x=.04, y=.97, s=9, size=K.T_SMALL)
K.frame(ax)
C.title(140, 7, "Rerun decoys (set 1)", size=TT, fontweight="bold"); C.letter(132.5, 6.5, "c")
C.save("SupplementaryFigure5")
print("FigS3 |", out, "|", rows)
