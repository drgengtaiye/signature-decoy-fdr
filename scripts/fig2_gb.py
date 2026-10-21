"""Figure 2 (Genome Biology) - the signature library in one view.

The v6 figure without the attribution and triage content, which went to a separate paper:
the atlas-cell and triage rings of panel a and the family-to-outcome alluvial (old e) are gone,
a ring for each signature's declared direction takes their place, and panel c moves to the
bottom row, drawn as bars.

a  every signature across every analysis, as concentric tracks (hero)
b  the genes the signatures share: network, edges weighted by Jaccard overlap
c  effective number of independent tests in each checkpoint cohort (Li and Ji)
d  the correlation between signature scores, averaged over the checkpoint cohorts, clustered
"""
import numpy as np, pandas as pd, glob, os, warnings; warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import networkx as nx
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
I = f"{K.D}/icb"
EFF = pd.read_csv(f"{I}/v4_standardised_effects.csv")
X = pd.read_csv(f"{I}/v3_dissociation.csv")[["signature", "family", "n_genes", "method", "d_resp"]].merge(EFF[["signature", "z_resp", "z_prog", "z_shift"]], on="signature")
# declared direction, as carried into the target-decoy analysis (+1 higher in responders, -1 lower)
DIRN = pd.concat([pd.read_csv(f, usecols=["unit", "kind", "dirn"]) for f in sorted(glob.glob(f"{I}/calib_v3/ICB_*_v3.csv"))])
DIRN = DIRN[DIRN.kind == "target"].drop_duplicates("unit").rename(columns={"unit": "signature"})[["signature", "dirn"]]
X = X.merge(DIRN, on="signature", how="left")
assert X.dirn.notna().all() and len(X) == 52
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
ENGC = {"GSVA": fs.SLATE, "ssGSEA": fs.TEAL, "Weighted Mean": fs.AMBER}
DIRC = {1: fs.INK, -1: "#B9C2CB"}
C = K.Canvas(212)

# ============================================================ a  the library in one view
order = []
for f in FAMORD: order += X[X.family == f].sort_values("z_resp", ascending=False).signature.tolist()
Xo = X.set_index("signature").loc[order].reset_index()
n = len(Xo); GAP = np.deg2rad(4.5)
fam_counts = [int((Xo.family == f).sum()) for f in FAMORD]
span = 2 * np.pi - np.deg2rad(22) - 2 * GAP             # leave a wedge at the top for the track key
step = span / n; th = []; t0 = np.deg2rad(11)
for i in range(n):
    extra = GAP * (int(i >= fam_counts[0]) + int(i >= fam_counts[0] + fam_counts[1]))
    th.append(t0 + (i + .5) * step + extra)
th = np.array(th); wdt = step * .86
ax = C.fig.add_axes([2 / C.w, 1 - (4 + 118) / C.h, 118 / C.w, 118 / C.h], projection="polar")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_ylim(0, 4.05); ax.axis("off")
def ring(r0, h, cols): ax.bar(th, h, width=wdt, bottom=r0, color=cols, lw=0)
ring(.95, .08, [fs.FAM[f] for f in Xo.family])
ring(1.06, .08, [ENGC[m] for m in Xo.method])
lg = np.log10(Xo.n_genes.clip(lower=3)); ring(1.17, .06 + .22 * (lg - np.log10(3)) / (np.log10(430) - np.log10(3)), fs.SLATE)
def diverging(base, v, scale, c_pos, c_neg):
    h = np.clip(v, -1.2, 2.4) * scale
    ax.bar(th, h, width=wdt, bottom=base, color=[c_pos if x >= 0 else c_neg for x in v], lw=0)
    tt = np.linspace(0, 2 * np.pi, 400); ax.plot(tt, np.full_like(tt, base), color=fs.INK, lw=.35)
diverging(1.62, Xo.z_resp.values, .10, fs.BLUE, fs.GREY)
diverging(2.02, Xo.z_prog.values, .10, fs.VERM, fs.GREY)
SHIFTC = "#6A51A3"
diverging(2.42, Xo.z_shift.values, .10, SHIFTC, fs.GREY)
ring(2.74, .10, [DIRC[d] for d in Xo.dirn])
for t, s, f in zip(th, Xo.signature, Xo.family):
    deg = np.rad2deg(t); rot = 90 - deg
    ha = "left"
    if 180 < deg < 360: rot += 180; ha = "right"
    ax.text(t, 2.95, K.lab(s), rotation=rot, rotation_mode="anchor", ha=ha, va="center", fontsize=5.0,
            color=fs.FAM[f] if f != "Other" else fs.SLATE)
for rr, lab in [(.99, "family"), (1.10, "engine"), (1.26, "genes"), (1.62, "response"), (2.02, "prognosis"), (2.42, "shift"),
                (2.79, "direction")]:
    ax.text(0, rr, lab, ha="center", va="center", fontsize=4.4, color=fs.INK, rotation=0)
ax.text(0, 0, "52\nsignatures", ha="center", va="center", fontsize=K.T_NOTE, color=fs.INK)
C.letter(1.5, 3, "a")
kx = C.ax(4, 123, 114, 8); kx.axis("off"); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
nd = X.dirn.value_counts()
rows = [[("family", [(fs.FAM["Inflamed/IFN"], "inflamed/IFN"), (fs.FAM["Other"], "other"), (fs.FAM["Confounder"], "confounder")]),
         ("engine", [(fs.SLATE, "GSVA"), (fs.TEAL, "ssGSEA"), (fs.AMBER, "weighted mean")])],
        [("z > 0", [(fs.BLUE, "response"), (fs.VERM, "prognosis"), (SHIFTC, "shift"), (fs.GREY, "any z < 0")]),
         ("declared direction", [(DIRC[1], f"higher in responders ({nd[1]})"), (DIRC[-1], f"lower ({nd[-1]})")])]]
rend = C.fig.canvas.get_renderer(); kw_ = kx.get_window_extent(rend).width
def tw(t_): return t_.get_window_extent(rend).width / kw_       # measured text width in key-axis units
for ri, row in enumerate(rows):
    yy = .80 - ri * .55; xx = 0
    for head, its in row:
        if head: xx += tw(kx.text(xx, yy, head, fontsize=5.0, va="center", fontweight="bold")) + .012
        for c, lab_ in its:
            kx.scatter([xx + .005], [yy], s=9, c=c, marker="s", lw=0)
            xx += .016 + tw(kx.text(xx + .016, yy, lab_, fontsize=5.0, va="center")) + .016
        xx += .025

# ============================================================ b  shared genes
genes = {}
for s in X.signature:
    f = f"{I}/SignatureSets-main/data-raw/IO/{s}.csv"
    if os.path.exists(f):
        g = pd.read_csv(f, encoding="latin-1").get("gene_name", pd.Series(dtype=str)).dropna()
        genes[s] = set(g)
Gx = nx.Graph()
for s in genes: Gx.add_node(s)
sl = list(genes); J = []
for i in range(len(sl)):
    for j in range(i + 1, len(sl)):
        a, b = genes[sl[i]], genes[sl[j]]
        if not a or not b: continue
        jac = len(a & b) / len(a | b)
        if jac >= .08: Gx.add_edge(sl[i], sl[j], w=jac); J.append(jac)
comps = sorted([c for c in nx.connected_components(Gx) if len(c) > 1], key=len, reverse=True)
famidx = lambda s: FAMORD.index(X.set_index("signature").family[s])
def fit(p_, x0, x1, y0, y1):
    a_ = np.array(list(p_.values())); lo, hi = a_.min(0), a_.max(0); sc = min((x1 - x0) / max(hi[0] - lo[0], 1e-9), (y1 - y0) / max(hi[1] - lo[1], 1e-9))
    mid = (lo + hi) / 2; return {k_: np.array([(x0 + x1) / 2, (y0 + y1) / 2]) + (v_ - mid) * sc for k_, v_ in p_.items()}
# largest component: stress layout with graph distance; smaller components: a row beneath; isolates: two rows at the bottom
Hc = nx.Graph(); Hc.add_nodes_from(sorted(comps[0])); Hc.add_edges_from(sorted(tuple(sorted(e)) for e in Gx.subgraph(comps[0]).edges))
dsp = dict(nx.shortest_path_length(Hc))  # stress layout on sqrt graph distance: spreads the dense core, shortens pendant chains
pos = fit(nx.kamada_kawai_layout(Hc, dist={u: {v: np.sqrt(dsp[u][v]) for v in Hc} for u in Hc}), .06, .94, .50, .98)  # sorted input keeps the layout reproducible
small = sorted(comps[1:], key=lambda c: (-len(c), sorted(c))); xw = 1 / max(len(small), 1)
for k, c in enumerate(small):
    sub = nx.Graph(); sub.add_nodes_from(sorted(c)); sub.add_edges_from(sorted(tuple(sorted(e)) for e in Gx.subgraph(c).edges)); p_ = nx.circular_layout(sub) if len(c) > 2 else {u: np.array([i_, 0.]) for i_, u in enumerate(sorted(c))}
    pos.update(fit(p_, k * xw + xw * .2, (k + 1) * xw - xw * .2, .31, .42 if len(c) > 2 else .365))
iso = sorted([s for s in Gx.nodes if Gx.degree(s) == 0], key=lambda s: (famidx(s), s))
for k, s in enumerate(iso):
    pos[s] = np.array([.06 + .88 * (k % 5) / 4, .22 - .055 * (k // 5)])
xs_ = np.array([p_[0] for p_ in pos.values()]); ys_ = np.array([p_[1] for p_ in pos.values()])
ax = C.ax(124, 6, 54, 112)
for u, v, dd in Gx.edges(data=True):
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color=fs.SLATE, lw=.25 + 2.2 * dd["w"], alpha=.25 + .6 * dd["w"], zorder=1)
fam_s = X.set_index("signature").family; ng = X.set_index("signature").n_genes
for s in Gx.nodes:
    ax.scatter(pos[s][0], pos[s][1], s=8 + 5 * np.log2(ng[s]), c=fs.FAM[fam_s[s]], lw=.4, ec="white", zorder=3)
deg = dict(Gx.degree()); hubs = sorted(deg, key=lambda k: -deg[k])[:2]   # the two most connected (17 and 13 edges; the next three tie at 12)
cen = np.mean([pos[s] for s in Hc.nodes], axis=0); Rc = np.percentile([np.linalg.norm(pos[s] - cen) for s in Hc.nodes], 80)
for k, s in enumerate(hubs):
    d_ = pos[s] - cen
    ax.annotate(K.lab(s), pos[s], xytext=cen + (Rc + .04) * d_ / max(np.linalg.norm(d_), 1e-9), ha="right" if d_[0] < 0 else "left",
                va="center", fontsize=5.0, zorder=5, arrowprops=dict(arrowstyle="-", lw=.3, color=fs.INK, shrinkA=1, shrinkB=2),
                bbox=dict(fc="white", ec="none", pad=.3, alpha=.85), annotation_clip=False)
ax.text(.5, .015, "isolated (no Jaccard ≥ 0.08 with any other set)", ha="center", va="top", fontsize=5.0, color=fs.SLATE)
ax.set_xlim(-.04, 1.04); ax.set_ylim(-.01, 1.0)
fs.clean(ax)
C.letter(118, 4, "b")

# ============================================================ c  effective tests
P = pd.read_csv(f"{I}/v4_permutation_inference.csv"); pc = P[P.analysis.str.startswith("per-cohort") & (P.cohort != "ICB_VanDenEnde")].sort_values("meff_li_ji")
ax = C.ax(25, 142, 58, 56); y = np.arange(len(pc))
tum = [fs.ICB_TUM[c] for c in pc.cohort]
ax.barh(y, pc.meff_li_ji, height=.66, color=[fs.TUMCOL[t] for t in tum], lw=0)
for yi, v in zip(y, pc.meff_li_ji): ax.text(v - .8, yi, f"{v:.0f}", va="center", ha="right", fontsize=5.0, color="white", fontweight="bold")
ax.axvline(pc.meff_li_ji.median(), color=fs.INK, lw=.6, ls=(0, (3, 2)))
ax.set_yticks(y); ax.set_yticklabels([f"{c.replace('ICB_', '').replace('_', ' ')} ({int(n_)})" for c, n_ in zip(pc.cohort, pc.n)], fontsize=5.0)
ax.tick_params(axis="y", length=0); ax.set_ylim(-.6, len(pc) - .4)
ax.set_xlim(0, 52); ax.set_xticks([0, 17, 52]); ax.set_xlabel("effective tests (of 52)", fontsize=K.T_SMALL); K.frame(ax, left=False)
K.key(ax, [(fs.TUMCOL[t], t) for t in ["melanoma", "kidney", "lung", "bladder", "gastric"]], .70, .40, dy=.075, size=5.0, marker="s", s=10, xpad=.03)
C.letter(1.5, 134, "c")

# ============================================================ d  score correlation, clustered
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
dn = dendrogram(Z, no_plot=True); leaves = dn["leaves"]
axd = C.ax(104, 140, 66, 8)
dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color=fs.SLATE, no_labels=True); fs.clean(axd)
for coll in axd.collections: coll.set_linewidth(.4)
axd.axhline(.5, color=fs.VERM, lw=.5, ls=(0, (2, 1.5)))
ax = C.ax(104, 149, 66, 56)
M = Rm[np.ix_(leaves, leaves)]
ax.imshow(M, cmap=fs.DIV, norm=Normalize(-1, 1), aspect="auto", interpolation="nearest")
lab_cl = cl[leaves]; start = 0
for k in range(1, len(leaves) + 1):
    if k == len(leaves) or lab_cl[k] != lab_cl[start]:
        if k - start > 1: ax.add_patch(Rectangle((start - .5, start - .5), k - start, k - start, fill=False, ec=fs.INK, lw=.6))
        start = k
fs.clean(ax)
fsx = C.ax(100.5, 149, 2.2, 56); fsx.set_xlim(0, 1); fsx.set_ylim(len(leaves) - .5, -.5); fsx.axis("off")
for i, li in enumerate(leaves): fsx.add_patch(Rectangle((0, i - .5), 1, 1, fc=fs.FAM[X.set_index("signature").family[sigs[li]]], ec="none"))
cax = C.ax(172, 166, 1.3, 22)
cb = C.fig.colorbar(plt.cm.ScalarMappable(Normalize(-1, 1), fs.DIV), cax=cax, ticks=[-1, 0, 1]); cb.outline.set_linewidth(.3)
cax.tick_params(labelsize=5.0, length=1.2, pad=1); cax.set_title("r", fontsize=5.0, pad=2)
C.letter(94, 134, "d")

C.save("Figure2_GB")
print(f"Figure2_GB | n {n} | direction {nd.to_dict()} | network nodes {Gx.number_of_nodes()} edges {Gx.number_of_edges()} Jaccard range {min(J):.2f}-{max(J):.2f} hubs {hubs}"
      f" | meff median {pc.meff_li_ji.median()} range {pc.meff_li_ji.min()}-{pc.meff_li_ji.max()} | clusters {len(set(cl))}")
