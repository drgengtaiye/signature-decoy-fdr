"""Figure 2 (GB v10) - the 52 signatures are one correlated library, not 52 independent tests.

a  every signature in one view, as concentric rings: family, scoring engine, gene count, response and prognostic
   effects (null standard deviations), declared direction
b  shared genes: signatures joined where the Jaccard index of their gene lists is at least 0.08
c  family to correlation cluster to response call (folded statistic, ten decoy sets pooled)
d  correlation between signature scores, averaged over the 12 checkpoint cohorts (Fisher z, weights n - 3), clustered
e  how alike the correlation structure is between cohorts
f  eigenvalues of the 52-signature correlation matrix in each cohort
g  effective number of independent tests in each checkpoint cohort (Li and Ji)
h  pairwise score correlations within and between clusters
i  shared genes against score correlation, every pair of signatures
j  the genes shared by the most signatures, and which signatures contain them
Data as in fig2_cns.py, plus v12_pkg/response_calls.csv.gz (calls).
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, glob, os, warnings
from collections import Counter
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle, FancyArrowPatch
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import networkx as nx
import cns10 as K
import figstyle as fs
K.apply()
I = f"{K.D}/icb"
EFF = pd.read_csv(f"{I}/v4_standardised_effects.csv")
X = pd.read_csv(f"{I}/v3_dissociation.csv")[["signature", "family", "n_genes", "method"]].merge(
    EFF[["signature", "z_resp", "z_prog"]], on="signature")
DIRN = pd.concat([pd.read_csv(f, usecols=["unit", "kind", "dirn"]) for f in sorted(glob.glob(f"{I}/calib_v3/ICB_*_v3.csv"))])
DIRN = DIRN[DIRN.kind == "target"].drop_duplicates("unit").rename(columns={"unit": "signature"})[["signature", "dirn"]]
X = X.merge(DIRN, on="signature", how="left"); assert X.dirn.notna().all() and len(X) == 52
RC = pd.read_csv(f"{I}/v12_pkg/response_calls.csv.gz")
CALL = RC[(RC.kind == "target") & (RC.stat == "fold")].set_index("target").call
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
fam_of = X.set_index("signature").family; ng = X.set_index("signature").n_genes
TUMC = {"melanoma": K.BLUE, "kidney": K.ORANGE, "lung": K.AQUA, "bladder": K.SLATE, "gastric": "#b8862b"}

# ---------------------------------------------------------------- score correlation and clusters (as fig2_cns.py)
cohorts = sorted(os.path.basename(p)[:-11] for p in glob.glob(f"{I}/scores/*_sig.csv.gz"))
cohorts = [c for c in cohorts if c != "ICB_VanDenEnde"]
sigs = sorted(EFF.signature); S_IDX = {s: i for i, s in enumerate(sigs)}; nS = len(sigs)
zsum = np.zeros((nS, nS)); wsum = np.zeros((nS, nS))
for co in cohorts:
    S = pd.read_csv(f"{I}/scores/{co}_sig.csv.gz"); cols = [c for c in S.columns if c in S_IDX]
    r = np.clip(np.nan_to_num(np.corrcoef(S[cols].to_numpy(float), rowvar=False), nan=0.0), -.999999, .999999)
    idx = [S_IDX[c] for c in cols]; w = max(len(S) - 3, 1)
    zsum[np.ix_(idx, idx)] += np.arctanh(r) * w; wsum[np.ix_(idx, idx)] += w
Rm = np.tanh(np.divide(zsum, wsum, out=np.zeros_like(zsum), where=wsum > 0)); np.fill_diagonal(Rm, 1.0)
dist = 1 - np.abs(Rm); np.fill_diagonal(dist, 0); dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average"); cl = fcluster(Z, t=.5, criterion="distance")
CLS = pd.Series(cl, index=sigs)
leaves = dendrogram(Z, no_plot=True)["leaves"]

C = K.Canvas(230)

# ================================================================== a  circular overview
order = []
for f in FAMORD:
    order += X[X.family == f].sort_values("z_resp", ascending=False).signature.tolist()
Xo = X.set_index("signature").loc[order].reset_index()
n = len(Xo); GAP = np.deg2rad(5)
fc = [int((Xo.family == f).sum()) for f in FAMORD]
span = 2 * np.pi - np.deg2rad(26) - 2 * GAP
step = span / n; t0 = np.deg2rad(13)
th = np.array([t0 + (i + .5) * step + GAP * (int(i >= fc[0]) + int(i >= fc[0] + fc[1])) for i in range(n)])
wd = step * .84
ax = C.fig.add_axes([1 / C.w, 1 - (3 + 100) / C.h, 100 / C.w, 100 / C.h], projection="polar")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_ylim(0, 3.75); ax.axis("off")
ENG = {"GSVA": K.INK, "ssGSEA": K.SLATE, "Weighted Mean": K.LGRAY}


def ring(r0, h, cols):
    ax.bar(th, h, width=wd, bottom=r0, color=cols, lw=0)


def guide(r, lw=0.35, c=K.RULE):
    for a0, a1 in ((t0 - step * .3, t0 + fc[0] * step + step * .3),
                   (t0 + fc[0] * step + GAP - step * .3, t0 + (fc[0] + fc[1]) * step + GAP + step * .3),
                   (t0 + (fc[0] + fc[1]) * step + 2 * GAP - step * .3, t0 + n * step + 2 * GAP + step * .3)):
        tt = np.linspace(a0, a1, 200); ax.plot(tt, np.full_like(tt, r), color=c, lw=lw, zorder=1)


ring(0.92, 0.10, [K.FAM[f] for f in Xo.family])
ring(1.06, 0.10, [ENG[m] for m in Xo.method])
lg = np.log10(Xo.n_genes.clip(lower=3))
guide(1.20)
ring(1.20, .04 + .30 * (lg - np.log10(3)) / (np.log10(430) - np.log10(3)), K.SLATE)
for base, v, cpos in ((1.86, Xo.z_resp.values, K.BLUE), (2.30, Xo.z_prog.values, K.ORANGE)):
    hgt = np.clip(v, -1.2, 2.4) * .12
    ax.bar(th, hgt, width=wd, bottom=base, color=[cpos if x >= 0 else K.GRAY for x in v], lw=0, zorder=3)
    guide(base, lw=0.4, c=K.INK)
    guide(base + 2 * .12, lw=0.3); guide(base - 1 * .12, lw=0.3)
ring(2.66, 0.10, [K.INK if d == 1 else K.LGRAY for d in Xo.dirn])
for t, s, f, c_ in zip(th, Xo.signature, Xo.family, Xo.signature.map(CALL)):
    deg = np.rad2deg(t); rot = 90 - deg; ha = "left"
    if 180 < deg < 360:
        rot += 180; ha = "right"
    ax.text(t, 2.86, K.lab(s), rotation=rot, rotation_mode="anchor", ha=ha, va="center", fontsize=4.5,
            color=K.FAM[f] if f != "Other" else K.INK2, fontweight="bold" if c_ else "normal")
for rr, lab_ in [(0.97, "family"), (1.11, "engine"), (1.32, "genes"), (1.90, "response"), (2.34, "prognosis"),
                 (2.71, "direction")]:
    ax.text(0, rr, lab_, ha="center", va="center", fontsize=4.4, color=K.INK)
ax.text(0, 0, "52\nsignatures", ha="center", va="center", fontsize=6.5, color=K.INK, fontweight="bold", linespacing=1.0)
C.letter(1.5, 3.5, "a")
# key
kx = C.ax(4, 106.5, 102, 9); kx.axis("off"); kx.set_xlim(0, 102); kx.set_ylim(9, 0)
nd = X.dirn.value_counts()
rows_ = [[("family", [(K.FAM["Inflamed/IFN"], "inflamed/IFN"), (K.FAM["Other"], "other"), (K.FAM["Confounder"], "confounder")]),
          ("engine", [(K.INK, "GSVA"), (K.SLATE, "ssGSEA"), (K.LGRAY, "weighted mean")])],
         [("effect, null SD", [(K.BLUE, "response > 0"), (K.ORANGE, "prognosis > 0"), (K.GRAY, "< 0")]),
          ("direction", [(K.INK, f"higher in responders ({nd[1]})"), (K.LGRAY, f"lower ({nd[-1]})")])]]
rend = C.fig.canvas.get_renderer()
mm_per_px = 102 / kx.get_window_extent(rend).width
for ri, row in enumerate(rows_):
    yy = 2.0 + ri * 4.2; xx = 0
    for head, its in row:
        t_ = kx.text(xx, yy, head, fontsize=K.T_SMALL, va="center", fontweight="bold")
        xx += t_.get_window_extent(rend).width * mm_per_px + 1.6
        for c, l_ in its:
            kx.add_patch(Rectangle((xx, yy - 1.0), 2.0, 2.0, fc=c, ec="none"))
            t_ = kx.text(xx + 2.8, yy, l_, fontsize=K.T_SMALL, va="center")
            xx += 2.8 + t_.get_window_extent(rend).width * mm_per_px + 2.2
        xx += 2.5
C.text(4, 115.8, "bold: called for response (folded, ten decoy sets)", size=K.T_SMALL, color=K.INK2)

# ================================================================== b  shared genes
genes = {}
for s in X.signature:
    f = f"{I}/SignatureSets-main/data-raw/IO/{s}.csv"
    if os.path.exists(f):
        genes[s] = set(pd.read_csv(f, encoding="latin-1").get("gene_name", pd.Series(dtype=str)).dropna())
Gx = nx.Graph(); Gx.add_nodes_from(sorted(genes))
sl = sorted(genes)
for i in range(len(sl)):
    for j in range(i + 1, len(sl)):
        a_, b_ = genes[sl[i]], genes[sl[j]]
        if a_ and b_:
            jac = len(a_ & b_) / len(a_ | b_)
            if jac >= .08:
                Gx.add_edge(sl[i], sl[j], w=jac)
iso = sorted([s for s in Gx.nodes if Gx.degree(s) == 0], key=lambda s: (FAMORD.index(fam_of[s]), s))
H = Gx.subgraph([s for s in Gx.nodes if Gx.degree(s) > 0]).copy()
comps = sorted(nx.connected_components(H), key=lambda c: (-len(c), sorted(c)))
big = nx.Graph(); big.add_nodes_from(sorted(comps[0])); big.add_edges_from(sorted(tuple(sorted(e)) for e in H.subgraph(comps[0]).edges))
dsp = dict(nx.shortest_path_length(big))
pos = nx.kamada_kawai_layout(big, dist={u: {v: np.sqrt(dsp[u][v]) for v in big} for u in big})


def fit(p_, x0, x1, y0, y1):
    a_ = np.array(list(p_.values())); lo, hi = a_.min(0), a_.max(0)
    sc = min((x1 - x0) / max(hi[0] - lo[0], 1e-9), (y1 - y0) / max(hi[1] - lo[1], 1e-9)); mid = (lo + hi) / 2
    return {k_: np.array([(x0 + x1) / 2, (y0 + y1) / 2]) + (v_ - mid) * sc for k_, v_ in p_.items()}


pos = fit(pos, 0.04, 0.70, 0.30, 0.98)
small = comps[1:]
for k, c in enumerate(small):
    sub = nx.Graph(); sub.add_nodes_from(sorted(c)); sub.add_edges_from(sorted(tuple(sorted(e)) for e in H.subgraph(c).edges))
    yy = 0.88 - k * 0.20
    if len(c) > 2:
        pos.update(fit(nx.circular_layout(sub), 0.80, 0.98, yy - 0.08, yy + 0.08))
    else:
        for i_, u in enumerate(sorted(c)):
            pos[u] = np.array([0.83 + 0.12 * i_, yy])
for k, s in enumerate(iso):
    pos[s] = np.array([0.04 + (k % 10) * 0.105, 0.15 - (k // 10) * 0.075])
ax = C.panel("b", 116, 9, 62, 50, head="Shared genes")
for u, v, dd in sorted(H.edges(data=True), key=lambda e: e[2]["w"]):
    p0, p1 = pos[u], pos[v]
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle="arc3,rad=0.12", arrowstyle="-", lw=0.25 + 2.4 * dd["w"],
                                 color=K.tint(K.SLATE, 0.75 - 0.6 * min(dd["w"] / 0.5, 1)), zorder=1))
for s in Gx.nodes:
    K.dots(ax, [pos[s][0]], [pos[s][1]], K.FAM[fam_of[s]], s=5 + 4.2 * np.log2(ng[s]), z=3)
deg = dict(H.degree()); hubs = sorted(deg, key=lambda k: (-deg[k], k))[:2]
cen = np.mean([pos[s] for s in big.nodes], axis=0)
for s in hubs:
    d_ = pos[s] - cen; d_ = d_ / max(np.linalg.norm(d_), 1e-9)
    ax.annotate(K.lab(s), pos[s], xytext=pos[s] + np.array([0.20 * np.sign(d_[0]), 0.10]), ha="right" if d_[0] < 0 else "left",
                va="center", fontsize=K.T_SMALL, arrowprops=dict(arrowstyle="-", lw=0.3, color=K.INK, shrinkA=1, shrinkB=2.5),
                annotation_clip=False, zorder=6)
ax.plot([0.02, 1.0], [0.215, 0.215], color=K.RULE, lw=0.35)
ax.text(0.02, 0.205, f"no Jaccard ≥ 0.08 with any other ({len(iso)})", fontsize=K.T_SMALL - 0.3, color=K.INK2, va="top")
ax.set_xlim(0, 1.02); ax.set_ylim(0.04, 1.02); ax.axis("off")
# size key
for k, gsz in enumerate((5, 20, 100)):
    K.dots(ax, [0.80 + k * 0.075], [0.30], K.GRAY, s=5 + 4.2 * np.log2(gsz), z=3)
    ax.text(0.80 + k * 0.075, 0.25, str(gsz), ha="center", va="top", fontsize=K.T_SMALL - 0.5, color=K.INK2)
ax.text(0.86, 0.36, "genes", ha="center", va="bottom", fontsize=K.T_SMALL - 0.3, color=K.INK2)

# ================================================================== c  alluvial: family -> cluster -> call
csize = CLS.value_counts()
big_cl = [c for c in csize.index if csize[c] >= 3]
def cl_name(c):
    return f"C{c}" if c in big_cl else "small"
rows = [(fam_of[s], cl_name(CLS[s]), "called" if CALL[s] else "not called") for s in order]
cl_order = [f"C{c}" for c in sorted(big_cl, key=lambda c: -csize[c])] + ["small"]
ax = C.panel("c", 124, 71, 46, 38, head="Family, cluster and call")
CLC = {f"C{c}": K.SLATE for c in big_cl}; CLC["small"] = K.LGRAY
colors_of = lambda k, c: (K.FAM[c] if k == 0 else CLC[c] if k == 1 else (K.INK if c == "called" else K.LGRAY))
pos_c = K.alluvial(ax, rows, [FAMORD, cl_order, ["called", "not called"]], colors_of, xs=[0.05, 0.5, 0.95], width=0.07,
                   gap=0.03, flow_color=lambda r: K.FAM[r[0]])
cnt = Counter(r[1] for r in rows); cntf = Counter(r[0] for r in rows); cntc = Counter(r[2] for r in rows)
for (k, c), (y0, y1) in pos_c.items():
    if k == 0:
        ax.text(0.05 - 0.05, (y0 + y1) / 2, f"{K.FAMLAB[c]}\n{cntf[c]}", ha="right", va="center", fontsize=K.T_SMALL - 0.2,
                linespacing=1.0)
    elif k == 1:
        ax.text(0.5 - 0.055, (y0 + y1) / 2, "≤ 2" if c == "small" else c, ha="right", va="center", fontsize=K.T_SMALL - 0.8,
                color=K.INK, zorder=5, bbox=dict(fc="white", ec="none", pad=0.15, alpha=0.8))
    else:
        ax.text(0.95 + 0.05, (y0 + y1) / 2, f"{c}\n{cntc[c]}", ha="left", va="center", fontsize=K.T_SMALL - 0.2, linespacing=1.0)
ax.set_xlim(-0.02, 1.02)
for xx_, t_ in ((0.05, "family"), (0.5, "cluster"), (0.95, "response")):
    ax.text(xx_, -0.03, t_, ha="center", va="top", fontsize=K.T_SMALL - 0.3, color=K.INK2)

# ================================================================== d  correlation heatmap
x0d, y0d, wd_, hd_ = 16, 135, 46, 46
C.letter(4.5, 123, "d"); C.head(9.5, 123.2, "Score correlation, 12 cohorts")
axd = C.ax(x0d, 127, wd_, 7.5); dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color=K.SLATE, no_labels=True)
for coll in axd.collections:
    coll.set_linewidth(0.35)
axd.axis("off"); axd.axhline(0.5, color=K.ORANGE, lw=0.5)
ax = C.ax(x0d, y0d, wd_, hd_)
M = Rm[np.ix_(leaves, leaves)]
ax.imshow(M, cmap=K.DIV, norm=Normalize(-1, 1), aspect="auto", interpolation="nearest")
lab_cl = cl[leaves]; start = 0
for k in range(1, len(leaves) + 1):
    if k == len(leaves) or lab_cl[k] != lab_cl[start]:
        if k - start > 1:
            ax.add_patch(Rectangle((start - .5, start - .5), k - start, k - start, fill=False, ec=K.INK, lw=0.55))
            if k - start >= 3:
                ax.text(k - .5, start - .7, f"C{lab_cl[start]}", fontsize=K.T_SMALL - 0.6, va="bottom", ha="right", color=K.INK)
        start = k
K.clean(ax)
for (xx, yy, ww, hh, vert) in ((x0d - 3.0, y0d, 2.2, hd_, True), (x0d, y0d + hd_ + 0.8, wd_, 2.2, False)):
    fsx = C.ax(xx, yy, ww, hh); fsx.axis("off")
    if vert:
        fsx.set_xlim(0, 1); fsx.set_ylim(len(leaves) - .5, -.5)
    else:
        fsx.set_ylim(0, 1); fsx.set_xlim(-.5, len(leaves) - .5)
    for i, li in enumerate(leaves):
        fsx.add_patch(Rectangle((0, i - .5) if vert else (i - .5, 0), 1, 1, fc=K.FAM[fam_of[sigs[li]]], ec="none"))
K.cbar(C, x0d + wd_ + 2.5, y0d + 12, 1.5, 22, K.DIV, -1, 1, [-1, 0, 1], label="r", horizontal=False)

# ================================================================== e  similarity of the correlation structure
iu = np.triu_indices(nS, 1)
EIG = {}; RV = {}
for co in cohorts:
    S = pd.read_csv(f"{I}/scores/{co}_sig.csv.gz"); cols = [c for c in sigs if c in S.columns]
    A = S[cols].to_numpy(float); ok = np.nanstd(A, 0) > 0; A = A[:, ok]
    r = np.corrcoef(A, rowvar=False)
    EIG[co] = np.sort(np.linalg.eigvalsh(r))[::-1]
    full = pd.DataFrame(np.nan, index=sigs, columns=sigs); keep = [c for c, o in zip(cols, ok) if o]
    full.loc[keep, keep] = r; RV[co] = full.to_numpy()[iu]
corder = sorted(cohorts, key=lambda c: (fs.ICB_TUM[c], K.COHORT[c]))
SIM = np.array([[pd.Series(RV[a]).corr(pd.Series(RV[b])) for b in corder] for a in corder])
ax = C.panel("e", 98, 135, 38, 38, head="Same structure in every cohort")
Ml = np.where(np.tril(np.ones_like(SIM), 0) > 0, SIM, np.nan)
ax.imshow(Ml, cmap=K.SEQ_BLUE, vmin=0.4, vmax=1, interpolation="nearest")
for i in range(len(corder)):
    for j in range(i):
        ax.text(j, i, f"{SIM[i, j]:.2f}"[1:], ha="center", va="center", fontsize=3.4,
                color="white" if SIM[i, j] > 0.78 else K.INK)
ax.set_xticks(range(len(corder))); ax.set_yticks(range(len(corder)))
ax.set_yticklabels([K.COHORT[c] for c in corder], fontsize=K.T_SMALL - 0.3)
ax.set_xticklabels([K.COHORT[c] for c in corder], fontsize=K.T_SMALL - 0.3, rotation=90)
ax.tick_params(length=0, pad=1); [sp.set_visible(False) for sp in ax.spines.values()]
for i, c in enumerate(corder):
    ax.add_patch(Rectangle((len(corder) - 0.35, i - .5), 0.6, 1, fc=TUMC[fs.ICB_TUM[c]], ec="white", lw=0.3, clip_on=False))
K.cbar(C, 122, 139, 12, 1.4, K.SEQ_BLUE, 0.4, 1.0, [0.4, 0.7, 1.0], label="r")
for k, t in enumerate(["melanoma", "kidney", "lung", "bladder", "gastric"]):
    C.fig.patches.append(Rectangle(((122 + 0) / C.w, 1 - (147.5 + k * 2.6) / C.h), 1.6 / C.w, 1.6 / C.h, fc=TUMC[t], ec="none",
                                   transform=C.fig.transFigure, figure=C.fig))
    C.text(124.4, 146.9 + k * 2.6, t, size=K.T_SMALL - 0.4)
offd = SIM[np.triu_indices(len(corder), 1)]

# ================================================================== f  eigenvalue spectra
ax = C.panel("f", 152, 135, 26, 38, head="Eigenvalues")
mx = max(len(v) for v in EIG.values())
E = np.array([np.r_[v, [np.nan] * (mx - len(v))] for v in EIG.values()])
xs_ = np.arange(1, mx + 1)
ax.fill_between(xs_, np.clip(np.nanmin(E, 0), 1e-3, None), np.nanmax(E, 0), color=K.tint(K.SLATE, 0.75), lw=0, zorder=1)
medev = np.nanmedian(E, 0)
ax.plot(xs_, np.clip(medev, 1e-3, None), color=K.INK, lw=1.0, zorder=3)
ax.axhline(1, color=K.ORANGE, lw=0.6, zorder=2)
ax.text(52, 1.4, "independent", ha="right", fontsize=K.T_SMALL - 0.3, color=K.ORANGE)
ax.set_yscale("log"); ax.set_ylim(1e-3, 50); ax.set_xlim(0.5, 52.5); ax.set_xticks([1, 26, 52])
ax.set_yticks([1e-3, 1e-1, 10]); ax.set_yticklabels(["0.001", "0.1", "10"]); ax.minorticks_off()
ax.set_xlabel("component"); ax.set_ylabel("eigenvalue")
pc1 = np.median([v[0] / v.sum() for v in EIG.values()])
ax.scatter([1], [medev[0]], s=10, color=K.INK, zorder=4)
ax.text(3, medev[0], f"PC1 {100 * pc1:.0f}%", va="center", fontsize=K.T_SMALL - 0.3)

# ================================================================== g  effective tests
P = pd.read_csv(f"{I}/v4_permutation_inference.csv")
pc = P[P.analysis.str.startswith("per-cohort") & (P.cohort != "ICB_VanDenEnde")].sort_values("meff_li_ji")
ax = C.panel("g", 21, 199, 26, 25, head="Effective tests")
y = np.arange(len(pc)); med = pc.meff_li_ji.median()
for yi, (c, v) in enumerate(zip(pc.cohort, pc.meff_li_ji)):
    ax.plot([0, v], [yi, yi], color=K.tint(TUMC[fs.ICB_TUM[c]], 0.35), lw=0.9, zorder=2, solid_capstyle="butt")
    K.dots(ax, [v], [yi], TUMC[fs.ICB_TUM[c]], s=12, z=3)
ax.axvline(med, color=K.INK, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
ax.text(med, len(pc) - 0.3, f"median {med:.0f}", ha="center", va="bottom", fontsize=K.T_SMALL - 0.3, color=K.INK2)
ax.set_yticks(y); ax.set_yticklabels([f"{K.COHORT[c]} ({int(n_)})" for c, n_ in zip(pc.cohort, pc.n)], fontsize=K.T_SMALL - 0.3)
K.frame(ax, left=False); ax.set_ylim(-0.6, len(pc) + 0.3)
ax.set_xlim(0, 52); ax.set_xticks([0, 26, 52]); ax.set_xlabel("effective tests (of 52)")

# ================================================================== h  within and between clusters
same = (cl[:, None] == cl[None, :])[iu]; rv = Rm[iu]
ax = C.panel("h", 61, 199, 28, 25, head="Pairs of signatures")
K.density(ax, rv[~same], K.GRAY, lo=-0.7, hi=1.0, alpha=0.6)
K.density(ax, rv[same], K.INK, lo=-0.2, hi=1.0, alpha=0.55)
for v_, c_ in ((np.median(rv[~same]), K.GRAY), (np.median(rv[same]), K.INK)):
    ax.axvline(v_, color=c_, lw=0.6, ls=(0, (2, 1.5)))
ax.set_xlabel("score correlation r"); ax.set_ylabel("density"); ax.set_ylim(0, None)
ax.set_xlim(-0.7, 1.0); ax.set_xticks([-0.5, 0, 0.5, 1])
ax.text(0.66, 0.97, f"same cluster\n({same.sum()})", transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL - 0.3,
        color=K.INK, linespacing=1.0)
ax.text(0.20, 0.75, f"different\n({(~same).sum()})", transform=ax.transAxes, ha="center", va="top", fontsize=K.T_SMALL - 0.3,
        color=K.INK2, linespacing=1.0)

# ================================================================== i  gene overlap vs score correlation
JAC = np.zeros(len(iu[0]))
for k, (a_, b_) in enumerate(zip(*iu)):
    A_, B_ = genes.get(sigs[a_], set()), genes.get(sigs[b_], set())
    JAC[k] = len(A_ & B_) / len(A_ | B_) if (A_ and B_) else np.nan
ax = C.panel("i", 101, 199, 32, 25, head="Correlated without shared genes")
z0 = JAC == 0
rng = np.random.default_rng(4)
ax.axvspan(-0.075, -0.005, color=K.GRID, lw=0, zorder=0)
for msk, c, a, s_ in ((~same, K.GRAY, 0.55, 1.6), (same, K.INK, 0.85, 2.2)):
    m0 = z0 & msk
    ax.scatter(rng.uniform(-0.068, -0.012, m0.sum()), rv[m0], s=s_, color=c, alpha=a, lw=0, zorder=3 if c == K.INK else 2)
    m1 = (JAC > 0) & msk; ax.scatter(JAC[m1], rv[m1], s=s_ + 1, color=c, alpha=a, lw=0, zorder=3 if c == K.INK else 2)
K.ref(ax, 0); ax.axhline(0.5, color=K.ORANGE, lw=0.5, ls=(0, (2, 1.5)))
ax.set_xlim(-0.08, 0.45); ax.set_ylim(-0.6, 1.05); ax.set_xticks([0, 0.2, 0.4])
ax.set_xlabel("shared genes (Jaccard)"); ax.set_ylabel("score correlation r")
ax.text(-0.04, -0.57, "none", ha="center", va="bottom", fontsize=K.T_SMALL - 0.5, color=K.INK2)
hi_r0 = int(((rv > 0.5) & z0).sum()); hi_r = int((rv > 0.5).sum())
ax.text(0.97, 0.05, f"{hi_r0} of {hi_r} pairs with\nr > 0.5 share no gene", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=K.T_SMALL - 0.3, color=K.INK, linespacing=1.0)

# ================================================================== j  most-shared genes: membership matrix
cnt_g = Counter(g for s_ in genes.values() for g in s_)
top = sorted(cnt_g.items(), key=lambda kv: (-kv[1], kv[0]))[:12]       # ties broken by name, reproducibly
holders = [s for s in order if s in genes and any(g in genes[s] for g, _ in top)]
ax = C.panel("j", 144, 199, 27, 25, head="Most-shared genes")
for i, (g, n_) in enumerate(top):
    ax.plot([-0.5, len(holders) - 0.5], [i, i], color=K.GRID, lw=0.4, zorder=0)
    for j, s in enumerate(holders):
        if g in genes[s]:
            ax.scatter([j], [i], s=3.2, color=K.FAM[fam_of[s]], lw=0, zorder=2)
ax.set_yticks(range(len(top))); ax.set_yticklabels([g for g, _ in top], fontsize=K.T_SMALL - 0.4, style="italic")
ax.set_ylim(len(top) - 0.5, -0.7); ax.set_xlim(-0.8, len(holders) - 0.2)
ax.set_xticks([]); K.frame(ax, left=False, bottom=False); ax.tick_params(axis="y", length=0, pad=1)
ax.set_xlabel(f"signatures containing them ({len(holders)})", labelpad=2)
axb = C.ax(172.5, 199, 6, 25); axb.barh(range(len(top)), [n_ for _, n_ in top], height=0.6, color=K.SLATE, lw=0)
axb.set_ylim(len(top) - 0.5, -0.7); axb.set_xlim(0, 15); axb.set_yticks([]); axb.set_xticks([])
K.frame(axb, left=False, bottom=False)
for i, (_, n_) in enumerate(top):
    axb.text(n_ + 0.4, i, str(n_), va="center", fontsize=K.T_SMALL - 0.9, color=K.INK2)

pd.DataFrame({"signature": order, "family": [fam_of[s] for s in order], "cluster": [int(CLS[s]) for s in order],
              "cluster_size": [int(csize[CLS[s]]) for s in order], "called_response_folded": [bool(CALL[s]) for s in order]}
             ).to_csv(f"{I}/v13_fig/family_cluster_call.csv", index=False)
pd.DataFrame([{"gene": g, "n_signatures": n_, "signatures": ";".join(s for s in order if s in genes and g in genes[s])}
              for g, n_ in top]).to_csv(f"{I}/v13_fig/most_shared_genes.csv", index=False)
C.save("Figure2_GB10")
print(f"Figure2_GB10 | called {int(CALL.sum())} | nodes {H.number_of_nodes()} edges {H.number_of_edges()} isolates {len(iso)} hubs {hubs} "
      f"| meff median {med} range {pc.meff_li_ji.min()}-{pc.meff_li_ji.max()} | clusters {len(set(cl))} sizes {sorted(csize.tolist(), reverse=True)} "
      f"| alluvial {Counter(rows)} | r same {np.median(rv[same]):.2f} diff {np.median(rv[~same]):.2f} | pc1 {pc1:.3f} "
      f"| sim median {np.median(offd):.2f} range {offd.min():.2f}-{offd.max():.2f} | r>0.5 no gene {hi_r0}/{hi_r} | top {top} "
      f"| holders {len(holders)}")
