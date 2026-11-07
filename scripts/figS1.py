"""Figure 2 redrawn at the density and in the look of the AML multi-omic figures (pixl_nc.py).
The 52 signatures as one correlated library. Old sheet: fig2_v10.py -> Figure2_GB10.
Tables: v3_dissociation.csv, v4_standardised_effects.csv (standardised_effects.py), v4_permutation_inference.csv (perm_inference.py),
calib_v3/ICB_*_v3.csv (declared direction), v12_pkg/response_calls.csv.gz (calls), scores/*_sig.csv.gz (signature scores),
SignatureSets-main/data-raw/IO/*.csv (gene lists). Correlation and clusters computed as in v5_clusters.py.
New panels are descriptive reads of those tables (validation/NOTE_figures_v17_descriptive_only.md)."""
import glob
import os
import sys
import warnings
from collections import Counter
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
import networkx as nx

warnings.filterwarnings("ignore")
I = D + "/icb"
COHORT = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Kim": "Kim", "ICB_Limagne1": "Limagne",
          "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao", "ICB_Puch": "Puch", "ICB_Riaz": "Riaz", "ICB_Van_Allen": "Van Allen"}
TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma", "ICB_Van_Allen": "melanoma",
       "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney", "ICB_Jung": "lung", "ICB_Limagne1": "lung",
       "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)",
         "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)", "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
FAMLAB = {"Inflamed/IFN": "inflamed/IFN", "Confounder": "confounder", "Other": "other"}
ENGORD = ["GSVA", "ssGSEA", "Weighted Mean"]
ENGLAB = {"GSVA": "GSVA", "ssGSEA": "ssGSEA", "Weighted Mean": "wt. mean"}
ENGC = {"GSVA": "#E9ECD2", "ssGSEA": "#B5C26C", "Weighted Mean": "#566F2F"}
assert check_palette("engine", ENGC)


def lab(s):
    if s in SHORT:
        return SHORT[s]
    a, b = s.rsplit("_", 1)
    return f"{a.replace('_', ' ')} ({b})"


def sgn(v, nd=2):
    return minus(f"{v:+.{nd}f}")


def vbar(pg, x, y, w, h, cmap, lo, hi, ticks, label):
    a = pg.ax(x, y, w, h)
    a.imshow(np.linspace(1, 0, 200)[:, None], aspect="auto", cmap=cmap, extent=(0, 1, lo, hi), rasterized=True)
    a.set_xticks([]); a.yaxis.tick_right()
    a.set_yticks(ticks, [minus(t) for t in ticks], fontsize=5)
    a.tick_params(length=1.5, pad=1, width=0.4)
    for s in a.spines.values():
        s.set_visible(False)
    pg.txt(x + w / 2, y - 1.5, label, ha="center")


# ---------------------------------------------------------------- tables, as in fig2_v10.py
EFF = pd.read_csv(f"{I}/v4_standardised_effects.csv")
T = pd.read_csv(f"{I}/v3_dissociation.csv")[["signature", "family", "n_genes", "method"]].merge(
    EFF[["signature", "z_resp", "z_prog"]], on="signature")
DIRN = pd.concat([pd.read_csv(f, usecols=["unit", "kind", "dirn"]) for f in sorted(glob.glob(f"{I}/calib_v3/ICB_*_v3.csv"))])
DIRN = DIRN[DIRN.kind == "target"].drop_duplicates("unit").rename(columns={"unit": "signature"})[["signature", "dirn"]]
T = T.merge(DIRN, on="signature", how="left")
assert T.dirn.notna().all() and len(T) == 52
RC = pd.read_csv(f"{I}/v12_pkg/response_calls.csv.gz")
CALL = RC[(RC.kind == "target") & (RC.stat == "fold")].set_index("target").call
assert int(CALL.sum()) == 21
T["called"] = T.signature.map(CALL).astype(bool)
fam_of = T.set_index("signature").family
ng = T.set_index("signature").n_genes
order = []
for f in FAMORD:
    order += T[T.family == f].sort_values("z_resp", ascending=False).signature.tolist()
TO = T.set_index("signature").loc[order]

# ---------------------------------------------------------------- score correlation and clusters (v5_clusters.py, fig2_v10.py)
cohorts = sorted(os.path.basename(p)[:-11] for p in glob.glob(f"{I}/scores/*_sig.csv.gz"))
cohorts = [c for c in cohorts if c != "ICB_VanDenEnde"]
assert len(cohorts) == 12
sigs = sorted(EFF.signature)
S_IDX = {s: i for i, s in enumerate(sigs)}
nS = len(sigs)
iu = np.triu_indices(nS, 1)
zsum, wsum = np.zeros((nS, nS)), np.zeros((nS, nS))
EIG, RV, RFULL = {}, {}, {}
for co in cohorts:
    S = pd.read_csv(f"{I}/scores/{co}_sig.csv.gz")
    cols = [c for c in S.columns if c in S_IDX]
    A = S[cols].to_numpy(float)
    r0 = np.corrcoef(A, rowvar=False)
    r = np.clip(np.nan_to_num(r0, nan=0.0), -.999999, .999999)
    idx = [S_IDX[c] for c in cols]
    w = max(len(S) - 3, 1)
    zsum[np.ix_(idx, idx)] += np.arctanh(r) * w
    wsum[np.ix_(idx, idx)] += w
    cols2 = [c for c in sigs if c in S.columns]
    A2 = S[cols2].to_numpy(float)
    ok = np.nanstd(A2, 0) > 0
    A2 = A2[:, ok]
    r2 = np.corrcoef(A2, rowvar=False)
    EIG[co] = np.sort(np.linalg.eigvalsh(r2))[::-1]
    full = pd.DataFrame(np.nan, index=sigs, columns=sigs)
    keep = [c for c, o in zip(cols2, ok) if o]
    full.loc[keep, keep] = r2
    RFULL[co] = full.to_numpy()
    RV[co] = RFULL[co][iu]
Rm = np.tanh(np.divide(zsum, wsum, out=np.zeros_like(zsum), where=wsum > 0))
np.fill_diagonal(Rm, 1.0)
dist = 1 - np.abs(Rm)
np.fill_diagonal(dist, 0)
dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average")
cl = fcluster(Z, t=.5, criterion="distance")
CLS = pd.Series(cl, index=sigs)
assert CLS.nunique() == 13
assert (CLS == pd.read_csv(f"{I}/v5_signature_clusters.csv").set_index("signature").cluster_fixed.reindex(CLS.index)).all()
leaves = dendrogram(Z, no_plot=True)["leaves"]
csize = CLS.value_counts()
big_cl = sorted([c for c in csize.index if csize[c] >= 3], key=lambda c: -csize[c])
rv = Rm[iu]
same = (cl[:, None] == cl[None, :])[iu]

# ---------------------------------------------------------------- gene lists and Jaccard (fig2_v10.py)
genes = {}
for s in T.signature:
    genes[s] = set(pd.read_csv(f"{I}/SignatureSets-main/data-raw/IO/{s}.csv", encoding="latin-1").get("gene_name", pd.Series(dtype=str)).dropna())
assert len(genes) == 52
sl = sorted(genes)
Gx = nx.Graph()
Gx.add_nodes_from(sl)
JM = pd.DataFrame(0.0, index=sl, columns=sl)
for i in range(len(sl)):
    for j in range(i + 1, len(sl)):
        a_, b_ = genes[sl[i]], genes[sl[j]]
        jac = len(a_ & b_) / len(a_ | b_) if (a_ and b_) else 0.0
        JM.loc[sl[i], sl[j]] = JM.loc[sl[j], sl[i]] = jac
        if jac >= .08:
            Gx.add_edge(sl[i], sl[j], w=jac)
iso = sorted([s for s in Gx.nodes if Gx.degree(s) == 0], key=lambda s: (FAMORD.index(fam_of[s]), s))
H = Gx.subgraph([s for s in Gx.nodes if Gx.degree(s) > 0]).copy()
comps = sorted(nx.connected_components(H), key=lambda c: (-len(c), sorted(c)))
deg = dict(H.degree())
hubs = sorted(deg, key=lambda k: (-deg[k], k))[:2]
JAC = np.array([JM.loc[sigs[a], sigs[b]] for a, b in zip(*iu)])
z0 = JAC == 0
hi_r0, hi_r = int(((rv > 0.5) & z0).sum()), int((rv > 0.5).sum())
cnt_g = Counter(g for s_ in genes.values() for g in s_)
top = sorted(cnt_g.items(), key=lambda kv: (-kv[1], kv[0]))[:12]
holders = [s for s in order if any(g in genes[s] for g, _ in top)]

P = pd.read_csv(f"{I}/v4_permutation_inference.csv")
PC = P[P.analysis.str.startswith("per-cohort") & (P.cohort != "ICB_VanDenEnde")].copy()
PC["lab"] = PC.cohort.map(COHORT)
PC = PC.sort_values(["meff_li_ji", "lab"], ascending=[False, True]).reset_index(drop=True)
med_meff = PC.meff_li_ji.median()
pc1 = np.median([v[0] / v.sum() for v in EIG.values()])
corder = sorted(cohorts, key=lambda c: (TUM[c], COHORT[c]))
SIM = np.array([[pd.Series(RV[a]).corr(pd.Series(RV[b])) for b in corder] for a in corder])

# ---- agreement with the deposited tables and with the old sheet
_s = pd.read_csv(f"{I}/v13_fig/cohort_similarity.csv", index_col=0).loc[corder, corder].to_numpy()
assert np.allclose(_s, SIM, atol=1e-6)
_e = pd.read_csv(f"{I}/v13_fig/eigen_by_cohort.csv")
for c in cohorts:
    assert np.allclose(_e[c].dropna().to_numpy()[:len(EIG[c])], EIG[c], atol=1e-6), c
_f = pd.read_csv(f"{I}/v13_fig/family_cluster_call.csv").set_index("signature")
assert (_f.cluster.reindex(sigs) == CLS).all() and (_f.called_response_folded.reindex(sigs) == T.set_index("signature").called.reindex(sigs)).all()
_m = pd.read_csv(f"{I}/v13_fig/most_shared_genes.csv")
assert list(_m.gene) == [g for g, _ in top] and list(_m.n_signatures) == [n for _, n in top]
for _, r_ in PC.iterrows():    # the table's mean |r| is the mean over the same per-cohort matrix drawn here
    assert abs(np.nanmean(np.abs(RV[r_.cohort])) - r_.mean_abs_r) < 1e-6, r_.cohort
assert set(iso) == {s for s in sl if JM.loc[s].drop(s).max() < .08}
assert (len(iso), round(med_meff, 6), round(pc1, 2), int(same.sum()), int((~same).sum()), hi_r0, hi_r, len(holders)) == (20, 17.0, 0.47, 353, 973, 183, 427, 23)

pg = Page(240)

# ================================================================ a  the 52 signatures, one row each
pg.panel(1, 3, "a", "The 52 signatures")
Y0, RH, GAP = 13.5, 2.0, 0.9
ROW, yy = {}, Y0
for f in FAMORD:
    y_start = yy
    for s in TO[TO.family == f].index:
        ROW[s] = yy
        yy += RH
    pg.txt(1.7, (y_start + yy) / 2, FAMLAB[f], ha="center", rotation=90)
    pg.line(3.0, y_start + 0.1, 3.0, yy - 0.25, color=INK, lw=0.4)
    yy += GAP
YEND = yy - GAP
XC = {"family": (26.6, 2.3), "engine": (28.9, 7.6), "genes": (36.5, 5.6), "direction": (42.1, 2.5), "response": (44.6, 6.6),
      "prognosis": (51.2, 6.6), "called": (57.8, 2.8), "cluster": (60.6, 4.4)}
for k, (x, w) in XC.items():
    pg.txt(x + w / 2 - 0.08, Y0 - 0.7, k, ha="center", va="bottom", rotation=90, weight="bold")
LOG0, LOG1 = np.log10(3), np.log10(434)
for s, r in TO.iterrows():
    y = ROW[s]
    pg.txt(26.2, y + RH / 2 - 0.07, lab(s), ha="right", weight="bold" if r.called else "normal")
    tile(pg, *[XC["family"][0], y], XC["family"][1], RH, FAM[r.family])
    tile(pg, XC["engine"][0], y, XC["engine"][1], RH, ENGC[r.method], ENGLAB[r.method])
    tile(pg, XC["genes"][0], y, XC["genes"][1], RH, BLUE((np.log10(r.n_genes) - LOG0) / (LOG1 - LOG0)), f"{r.n_genes:g}")
    tile(pg, XC["direction"][0], y, XC["direction"][1], RH, PUOR(0.80 if r.dirn == 1 else 0.20), "+" if r.dirn == 1 else "−")
    for k, v in (("response", r.z_resp), ("prognosis", r.z_prog)):
        tile(pg, XC[k][0], y, XC[k][1], RH, PUOR(float(np.clip(0.5 + v / 5, 0, 1))), sgn(v))
    point(pg.c, XC["called"][0] + XC["called"][1] / 2 - 0.08, y + RH / 2 - 0.07, SIG, ms=2.5, filled=bool(r.called))
    pg.txt(XC["cluster"][0] + XC["cluster"][1] / 2, y + RH / 2 - 0.07, f"C{CLS[s]}", ha="center")
key(pg, 3.5, YEND + 2.2, [(SIG, f"called ({int(TO.called.sum())})", "o"), (SIG, f"not called ({int((~TO.called).sum())})", "open-o")])
pg.cbar(37.5, YEND + 1.5, 16, 1.3, PUOR, -2.5, 2.5, [-2.5, 0, 2.5], "null SD")
pg.src("a", TO.reset_index().assign(cluster=lambda d: ["C%d" % CLS[s] for s in d.signature],
                                    label=lambda d: [lab(s) for s in d.signature])[
    ["signature", "label", "family", "method", "n_genes", "dirn", "z_resp", "z_prog", "called", "cluster"]])

# ================================================================ b  each signature's r with the other 51
RDF = pd.DataFrame(Rm, index=sigs, columns=sigs)
pg.panel(66.4, 3, "b", "Score r, other 51")
ax = pg.ax(68.5, Y0, 19, YEND - Y0)
out = []
for s in order:
    v = RDF.loc[s].drop(s).to_numpy()
    y = ROW[s] + RH / 2 - 0.07
    q = np.quantile(v, [0, 0.25, 0.5, 0.75, 1])
    ax.plot([q[0], q[4]], [y, y], color=OTHER, lw=0.35, solid_capstyle="butt", zorder=2)
    ax.plot([q[1], q[3]], [y, y], color=MUTED, lw=1.0, solid_capstyle="butt", zorder=3)
    point(ax, q[2], y, FAM[fam_of[s]], ms=2.3, zorder=5)
    out.append((s, *q))
zero(ax, 0)
ax.set_xlim(-0.8, 1.0); ax.set_ylim(YEND, Y0); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([-0.5, 0, 0.5, 1], ["−0.5", "0", "0.5", "1"]); ax.set_xlabel("score r")
pg.src("b", pd.DataFrame(out, columns=["signature", "r_min", "r_q25", "r_median", "r_q75", "r_max"]))

# ================================================================ c  each signature's largest Jaccard index
pg.panel(89.8, 3, "c", "Largest Jaccard")
ax = pg.ax(92, Y0, 15.5, YEND - Y0)
out = []
for s in order:
    v = JM.loc[s].drop(s)
    y = ROW[s] + RH / 2 - 0.07
    point(ax, v.max(), y, FAM[fam_of[s]], ms=2.3)
    out.append((s, v.max(), v.idxmax()))
zero(ax, 0.08)
ax.set_xlim(0, 0.62); ax.set_ylim(YEND, Y0); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([0, 0.2, 0.4, 0.6], ["0", "0.2", "0.4", "0.6"]); ax.set_xlabel("Jaccard")
pg.src("c", pd.DataFrame(out, columns=["signature", "largest_jaccard", "partner"]))


# ================================================================ d  network of shared genes
def fit(p_, x0, x1, y0, y1):
    a_ = np.array(list(p_.values()))
    lo, hi = a_.min(0), a_.max(0)
    sc = np.array([(x1 - x0) / max(hi[0] - lo[0], 1e-9), (y1 - y0) / max(hi[1] - lo[1], 1e-9)])
    mid = (lo + hi) / 2
    return {k_: np.array([(x0 + x1) / 2, (y0 + y1) / 2]) + (v_ - mid) * sc for k_, v_ in p_.items()}


big = nx.Graph()
big.add_nodes_from(sorted(comps[0]))
big.add_edges_from(sorted(tuple(sorted(e)) for e in H.subgraph(comps[0]).edges))
dsp = dict(nx.shortest_path_length(big))
pos = nx.kamada_kawai_layout(big, dist={u: {v: np.sqrt(dsp[u][v]) for v in big} for u in big})
pos = fit(pos, 0.03, 0.70, 0.04, 0.95)
for k, c in enumerate(comps[1:]):
    sub = nx.Graph()
    sub.add_nodes_from(sorted(c))
    sub.add_edges_from(sorted(tuple(sorted(e)) for e in H.subgraph(c).edges))
    yy_ = 0.86 - k * 0.24
    if len(c) > 2:
        pos.update(fit(nx.circular_layout(sub), 0.80, 0.97, yy_ - 0.10, yy_ + 0.10))
    else:
        for i_, u in enumerate(sorted(c)):
            pos[u] = np.array([0.82 + 0.13 * i_, yy_])
pg.panel(112, 3, "d", "Shared genes")
ax = pg.ax(113, 6, 66, 42)
segs = sorted(H.edges(data=True), key=lambda e: e[2]["w"])
ax.add_collection(LineCollection([[pos[u], pos[v]] for u, v, _ in segs], linewidths=[0.25 + 2.4 * d["w"] for _, _, d in segs],
                                 colors=[tint(MUTED, 0.3 + 0.7 * min(d["w"] / 0.5, 1)) for _, _, d in segs], zorder=1))
AREA = lambda n: 5 + 4.2 * np.log2(n)  # noqa: E731
for s in H.nodes:
    ax.scatter([pos[s][0]], [pos[s][1]], s=AREA(ng[s]), color=FAM[fam_of[s]], edgecolors="white", linewidths=0.3, zorder=3)
cen = np.mean([pos[s] for s in big.nodes], axis=0)
for s in hubs:
    d_ = pos[s] - cen
    sx = np.sign(d_[0]) if abs(d_[0]) > 1e-6 else 1
    ax.annotate(lab(s), pos[s], xytext=pos[s] + np.array([0.20 * sx, 0.20]), ha="right" if sx < 0 else "left", va="center", fontsize=5,
                arrowprops=dict(arrowstyle="-", lw=0.3, color=INK, shrinkA=1, shrinkB=2.5), annotation_clip=False, zorder=6)
for k, gsz in enumerate((5, 20, 100)):
    ax.scatter([0.80 + k * 0.075], [0.12], s=AREA(gsz), color=OTH, edgecolors="white", linewidths=0.3, zorder=3)
    ax.text(0.80 + k * 0.075, 0.045, str(gsz), ha="center", va="center", fontsize=5)
ax.text(0.875, 0.20, "genes", ha="center", va="center", fontsize=5)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
for yk, (c_, l_) in zip((39.6, 42.1, 44.6), ((INFL, "inflamed/IFN"), (OTH, "other"), (CONF, "confounder"))):
    key(pg, 146.5, yk, [(c_, l_, "o")])
pg.src("d", pd.DataFrame([{"signature_a": u, "signature_b": v, "jaccard": d["w"]} for u, v, d in segs]))

# ================================================================ e  signatures sharing none
pg.panel(112, 50.5, "e", "No Jaccard ≥ 0.08 with any other (20)")
for k, s in enumerate(iso):
    cx, cy = (113.6, 133.6, 156.0)[k // 7], 55.0 + (k % 7) * 2.6
    pg.c.scatter([cx], [cy], s=AREA(ng[s]), color=FAM[fam_of[s]], edgecolors="white", linewidths=0.3, zorder=3)
    pg.txt(cx + 1.9, cy, lab(s))
pg.src("e", pd.DataFrame({"signature": iso, "family": [fam_of[s] for s in iso], "n_genes": [ng[s] for s in iso],
                          "largest_jaccard": [JM.loc[s].drop(s).max() for s in iso]}))

# ================================================================ f  family, cluster and call: tile cross-table
pg.panel(112, 73.5, "f", "Family, cluster and call")
cname = {c: f"C{c}" for c in big_cl}
TO["cl_row"] = [cname.get(CLS[s], "≤ 2") for s in TO.index]
crow = [cname[c] for c in big_cl] + ["≤ 2"]
x0, y0, cw, ch = 119.6, 81.4, 6.5, 2.5
FX = {}
for i, f in enumerate(FAMORD):
    xf = x0 + cw * 2 * i
    pg.txt(xf + cw - 0.08, 77.2, f"{FAMLAB[f]}, {int((TO.family == f).sum())}", ha="center")
    pg.line(xf + 0.2, 78.2, xf + cw * 2 - 0.35, 78.2, lw=0.4)
    FX[f] = xf
    pg.txt(xf + cw / 2 - 0.08, 79.9, "called", ha="center", weight="bold")
    pg.txt(xf + cw * 1.5 - 0.08, 79.9, "not", ha="center", weight="bold")
xm = x0 + cw * 6 + 1.0
pg.txt(xm + cw - 0.08, 77.2, "cluster", ha="center")
pg.line(xm + 0.2, 78.2, xm + cw * 3 - 0.35, 78.2, lw=0.4)
for j, h_ in enumerate(("called", "not", "n")):
    pg.txt(xm + cw * j + cw / 2 - 0.08, 79.9, h_, ha="center", weight="bold")
out = []
for i, c in enumerate(crow):
    y = y0 + ch * i
    pg.txt(x0 - 0.8, y + ch / 2 - 0.07, c, ha="right")
    g = TO[TO.cl_row == c]
    rec = {"cluster": c}
    for f in FAMORD:
        for k, v in ((0, True), (1, False)):
            n = int(((g.family == f) & (g.called == v)).sum())
            tile(pg, FX[f] + cw * k, y, cw, ch, BLUE(float(np.sqrt(n / 16))), str(n))
            rec[f"{f} {'called' if v else 'not called'}"] = n
    for j, n in enumerate((int(g.called.sum()), int((~g.called).sum()), len(g))):
        tile(pg, xm + cw * j, y, cw, ch, PURPLE(float(np.sqrt(n / 26))), str(n))
    rec.update({"cluster called": int(g.called.sum()), "cluster not called": int((~g.called).sum()), "cluster n": len(g)})
    out.append(rec)
y = y0 + ch * len(crow)
pg.txt(x0 - 0.8, y + ch / 2 - 0.07, "all", ha="right", weight="bold")
for f in FAMORD:
    for k, v in ((0, True), (1, False)):
        pg.txt(FX[f] + cw * k + cw / 2 - 0.08, y + ch / 2 - 0.07, str(int(((TO.family == f) & (TO.called == v)).sum())), ha="center", weight="bold")
for j, n in enumerate((int(TO.called.sum()), int((~TO.called).sum()), len(TO))):
    pg.txt(xm + cw * j + cw / 2 - 0.08, y + ch / 2 - 0.07, str(n), ha="center", weight="bold")
pg.src("f", pd.DataFrame(out))

# ================================================================ g  family, engine and direction
pg.panel(112, 106.3, "g", "Family, engine and direction")
x0, y0, cw, ch = 130.6, 113.6, 8.4, 2.45
hdr = [ENGLAB[e] for e in ENGORD] + ["higher", "lower"]
pg.txt(x0 + cw * 1.5 - 0.08, 109.9, "engine", ha="center")
pg.line(x0 + 0.2, 110.8, x0 + cw * 3 - 0.35, 110.8, lw=0.4)
xdir = x0 + cw * 3 + 1.0
pg.txt(xdir + cw - 0.08, 109.9, "direction", ha="center")
pg.line(xdir + 0.2, 110.8, xdir + cw * 2 - 0.35, 110.8, lw=0.4)
xs = [x0 + cw * k for k in range(3)] + [xdir, xdir + cw]
for xx, h_ in zip(xs, hdr):
    pg.txt(xx + cw / 2 - 0.08, 112.3, h_, ha="center", weight="bold")
xn = xdir + cw * 2 + 0.6
pg.txt(xn + 2.5, 112.3, "n", ha="center", weight="bold")
out = []
for i, f in enumerate(FAMORD):
    y = y0 + ch * i
    g = TO[TO.family == f]
    pg.txt(x0 - 0.8, y + ch / 2 - 0.07, FAMLAB[f], ha="right")
    vals = [int((g.method == e).sum()) for e in ENGORD] + [int((g.dirn == 1).sum()), int((g.dirn == -1).sum())]
    for xx, n in zip(xs, vals):
        tile(pg, xx, y, cw, ch, BLUE(float(np.sqrt(n / 18))), str(n))
    pg.txt(xn + 2.5, y + ch / 2 - 0.07, str(len(g)), ha="center")
    out.append({"family": f, **dict(zip(hdr, vals)), "n": len(g)})
y = y0 + ch * 3
pg.txt(x0 - 0.8, y + ch / 2 - 0.07, "all", ha="right", weight="bold")
vals = [int((TO.method == e).sum()) for e in ENGORD] + [int((TO.dirn == 1).sum()), int((TO.dirn == -1).sum())]
for xx, n in zip(xs, vals):
    pg.txt(xx + cw / 2 - 0.08, y + ch / 2 - 0.07, str(n), ha="center", weight="bold")
pg.txt(xn + 2.5, y + ch / 2 - 0.07, str(len(TO)), ha="center", weight="bold")
out.append({"family": "all", **dict(zip(hdr, vals)), "n": len(TO)})
pg.src("g", pd.DataFrame(out))

# ================================================================ h  correlation matrix with 13 clusters
XH, YH, SH = 5.0, 137.0, 52.0
pg.panel(1, 127.5, "h", "Score correlation, 12 cohorts")
axd = pg.ax(XH, 130.6, SH, YH - 130.6 - 0.3)
dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color=MUTED, no_labels=True)
for coll in axd.collections:
    coll.set_linewidth(0.35)
axd.axis("off")
axd.axhline(0.5, color=INK, lw=0.5, ls=(0, (1.0, 1.5)))
ax = pg.ax(XH, YH, SH, SH)
M = Rm[np.ix_(leaves, leaves)]
ax.imshow(M, cmap=PUOR, norm=Normalize(-1, 1), aspect="auto", interpolation="nearest")
lab_cl = cl[leaves]
start = 0
blocks = []
for k in range(1, len(leaves) + 1):
    if k == len(leaves) or lab_cl[k] != lab_cl[start]:
        blocks.append((int(lab_cl[start]), start, k))
        if k - start > 1:
            ax.add_patch(Rectangle((start - .5, start - .5), k - start, k - start, fill=False, ec=INK, lw=0.55))
            ax.text(k - .5 + 0.2, start - .7, f"C{lab_cl[start]}", fontsize=5, va="bottom", ha="right", color=INK)
        start = k
# clusters of one have no box: named beside their cells, in leaf order
ax.text(7.6, 1.0, ", ".join(f"C{c}" for c, a_, b_ in blocks if b_ <= 3), fontsize=5, va="center", ha="left", color=INK)
ax.text(11.4, 8.5, ", ".join(f"C{c}" for c, a_, b_ in blocks if 7 <= a_ and b_ <= 11), fontsize=5, va="center", ha="left", color=INK)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ss = [sigs[li] for li in leaves]
for i, s in enumerate(ss):
    pg.box(XH - 2.6, YH + i, 2.0, 1.0, FAM[fam_of[s]], z=2)
    pg.box(XH + i, YH + SH + 0.6, 1.0, 2.0, FAM[fam_of[s]], z=2)
vbar(pg, XH + SH + 1.0, YH + 21, 1.4, 22, PUOR, -1, 1, [-1, 0, 1], "r")
pg.src("h", pd.DataFrame(M, index=ss, columns=ss).reset_index().rename(columns={"index": "signature"}))

# ================================================================ i  similarity between cohorts
pg.panel(63.5, 127.5, "i", "Similarity between cohorts")
x0, y0, cs = 77.4, 133.0, 3.72
out = []
for i, a in enumerate(corder):
    pg.txt(x0 - 1.8, y0 + cs * i + cs / 2 - 0.07, COHORT[a], ha="right")
    pg.box(x0 - 1.3, y0 + cs * i, 1.0, cs - 0.15, TYPE[TUM[a]], z=2)
    pg.txt(x0 + cs * i + cs / 2 - 0.07, y0 + cs * 12 + 0.5, COHORT[a], ha="center", va="top", rotation=90)
    for j, b in enumerate(corder[: i + 1]):
        if i == j:
            tile(pg, x0 + cs * j, y0 + cs * i, cs, cs, BLUE(1.0))
        else:
            tile(pg, x0 + cs * j, y0 + cs * i, cs, cs, BLUE(float(np.clip((SIM[i, j] - 0.4) / 0.6, 0, 1))), f"{SIM[i, j]:.2f}"[1:])
            out.append({"cohort_a": a, "cohort_b": b, "r": SIM[i, j]})
pg.cbar(x0 + cs * 5.8, y0 + 0.6, 18, 1.4, BLUE, 0.4, 1.0, [0.4, 0.7, 1.0], "r")
for k, t_ in enumerate(["melanoma", "kidney", "lung", "bladder", "gastric"]):
    pg.box(x0 + cs * 8.0, y0 + 6.6 + k * 2.3 - 0.8, 1.6, 1.6, TYPE[t_], z=2)
    pg.txt(x0 + cs * 8.0 + 2.2, y0 + 6.6 + k * 2.3, t_)
pg.src("i", pd.DataFrame(out))

# ================================================================ j  eigenvalue spectra
pg.panel(124.3, 127.5, "j", "Eigenvalues")
ax = pg.ax(131.5, 132, 18.5, 22.4)
mx = max(len(v) for v in EIG.values())
E = np.array([np.r_[v, [np.nan] * (mx - len(v))] for v in EIG.values()])
xs_ = np.arange(1, mx + 1)
medev = np.nanmedian(E, 0)
for v in E:
    ax.plot(xs_, np.clip(v, 1e-3, None), color=OTHER, lw=0.35, zorder=1)
ax.plot(xs_, np.clip(medev, 1e-3, None), color=INK, lw=0.9, zorder=3)
zero(ax, 1, axis="y")
ax.set_yscale("log"); ax.set_ylim(1e-3, 50); ax.set_xlim(0.5, 52.5); ax.set_xticks([1, 26, 52])
ax.set_yticks([1e-3, 1e-1, 10], ["0.001", "0.1", "10"]); ax.minorticks_off()
ax.set_xlabel("component"); ax.set_ylabel("eigenvalue")
point(ax, 1, medev[0], INK, ms=2.6, zorder=5)
ax.text(4.5, medev[0], f"PC1 {100 * pc1:.0f}%", va="center", fontsize=5)
ax.text(52.5, 1.35, "independent", ha="right", va="bottom", fontsize=5)
pg.src("j", pd.DataFrame(E.T, columns=list(EIG)).assign(component=xs_, median=medev))

# ================================================================ k  pairs in the same and in different clusters
pg.panel(153.3, 127.5, "k", "Pairs of signatures")
ax = pg.ax(160.5, 132, 18.5, 22.4)
bins = np.arange(-0.8, 1.0001, 0.05)
hist(ax, [rv[~same], rv[same]], bins, color=[NULLC, MUTED], stacked=True)
top_y = ax.get_ylim()[1]
for v_, c_ in ((np.median(rv[~same]), MUTED), (np.median(rv[same]), INK)):
    ax.plot([v_, v_], [0, top_y * 1.02], color=c_, lw=0.5, ls=(0, (1.0, 1.5)), zorder=5)
ax.set_ylim(0, top_y * 1.32)
ax.text(-0.74, top_y * 1.30, f"different\n({int((~same).sum())})", ha="left", va="top", fontsize=5, linespacing=1.05)
ax.text(1.0, top_y * 1.30, f"same cluster\n({int(same.sum())})", ha="right", va="top", fontsize=5, linespacing=1.05, color=INK)
ax.set_xlabel("score r"); ax.set_ylabel("pairs"); ax.set_xlim(-0.8, 1.0); ax.set_xticks([-0.5, 0, 0.5, 1], ["−0.5", "0", "0.5", "1"])
pg.src("k", pd.DataFrame({"a": [sigs[a] for a in iu[0]], "b": [sigs[b] for b in iu[1]], "r": rv, "same_cluster": same}))

# ================================================================ l  shared genes against score correlation
pg.panel(124.3, 162.3, "l", "Score r, shared genes")
ax = pg.ax(131.5, 167, 21, 20.6)
rng = np.random.default_rng(4)
jit = np.where(z0, rng.uniform(-0.068, -0.012, z0.size), JAC)
for msk, c_, s_, zz in ((~same, OTHER, 0.9, 2), (same, INK, 1.3, 3)):
    ax.scatter(jit[msk], rv[msk], s=s_, color=c_, lw=0, zorder=zz, rasterized=True)
zero(ax, 0, axis="y")
zero(ax, 0.5, axis="y")
ax.set_xlim(-0.08, 0.45); ax.set_ylim(-0.98, 1.05); ax.set_xticks([0, 0.2, 0.4], ["0", "0.2", "0.4"]); ax.set_yticks([-0.5, 0, 0.5, 1], ["−0.5", "0", "0.5", "1"])
ax.set_xlabel("Jaccard"); ax.set_ylabel("score r")
ax.text(-0.074, -0.96, "none", ha="left", va="bottom", fontsize=5)
ax.text(0.45, -0.96, f"{hi_r0} of {hi_r} pairs,\nr > 0.5, none", ha="right", va="bottom", fontsize=5, linespacing=1.05)
pg.src("l", pd.DataFrame({"a": [sigs[a] for a in iu[0]], "b": [sigs[b] for b in iu[1]], "jaccard": JAC, "r": rv, "same_cluster": same}))

# ================================================================ m  how many signatures each gene sits in
pg.panel(155, 162.3, "m", "Gene sharing")
ax = pg.ax(161.5, 167, 17.5, 20.6)
cv = np.array(list(cnt_g.values()))
hist(ax, cv, np.arange(0.5, cv.max() + 1.5), log=True)
ax.set_xlim(0.5, cv.max() + 0.5); ax.set_xticks([1, 5, 9, 13]); ax.set_ylim(0.7, len(cv) * 3)
ax.set_yticks([1, 10, 100, 1000], ["1", "10", "100", "1,000"]); ax.minorticks_off()
ax.set_xlabel("signatures"); ax.set_ylabel("genes")
ax.text(cv.max() + 0.5, len(cv) * 2.6, f"n = {len(cv):,}", ha="right", va="top", fontsize=5)
pg.src("m", pd.DataFrame({"gene": list(cnt_g), "n_signatures": cv}).sort_values(["n_signatures", "gene"], ascending=[False, True]))

# ================================================================ cohort band: rows shared by n, o and p
YC, RC_ = 200.8, 2.68
HC = RC_ * 12
pg.panel(1, 195.5, "n", "Effective tests")
ax = pg.ax(13.5, YC, 24, HC)
for i, r in PC.iterrows():
    y = YC + RC_ * i + RC_ / 2 - 0.07
    pg.txt(12.9, y, r.lab, ha="right")
    point(ax, r.meff_li_ji, y, TYPE[TUM[r.cohort]], ms=3.0)
    pg.txt(42.0, y, f"{r.n:.0f}", ha="right")
    pg.txt(47.0, y, f"{r.meff_li_ji:.0f}", ha="right")
zero(ax, med_meff)
pg.txt(13.5 + 24 * med_meff / 52, YC - 1.2, f"median {med_meff:.0f}", ha="center")
pg.txt(42.0, YC - 1.2, "n", ha="right", weight="bold")
pg.txt(47.0, YC - 1.2, "tests", ha="right", weight="bold")
ax.set_xlim(0, 52); ax.set_ylim(YC + HC, YC); ax.set_xticks([0, 26, 52]); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("effective tests (of 52)")
pg.src("n", PC[["cohort", "n", "meff_li_ji"]].assign(tumor=PC.cohort.map(TUM)))

pg.panel(48.8, 195.5, "o", "Abs. r, PC share")
x0, cw = 50.0, 5.65
hdr = ["abs. r", "PC1", "PC2", "PC3"]
for j, h_ in enumerate(hdr):
    pg.txt(x0 + cw * j + cw / 2 - 0.08, YC - 1.2, h_, ha="center", weight="bold")
out = []
for i, r in PC.iterrows():
    y = YC + RC_ * i
    v = EIG[r.cohort]
    sh = 100 * v[:3] / v.sum()
    tile(pg, x0, y, cw, RC_, BLUE(float(np.clip((r.mean_abs_r - 0.2) / 0.4, 0, 1))), f"{r.mean_abs_r:.2f}")
    for j, s_ in enumerate(sh):
        tile(pg, x0 + cw * (j + 1), y, cw, RC_, PURPLE(float(np.clip(s_ / 60, 0, 1))), f"{s_:.0f}")
    out.append({"cohort": r.cohort, "mean_abs_r": r.mean_abs_r, "PC1_pct": sh[0], "PC2_pct": sh[1], "PC3_pct": sh[2]})
pg.src("o", pd.DataFrame(out))

pg.panel(74.3, 195.5, "p", "Within-cluster r, by cohort")
x0, cw = 75.6, 5.65
cr = [cname[c] for c in big_cl]
for j, h_ in enumerate(cr + ["same", "diff."]):
    pg.txt(x0 + cw * j + (0.5 if j >= len(cr) else 0) + cw / 2 - 0.08, YC - 1.2, h_, ha="center", weight="bold")
out = []
for i, r in PC.iterrows():
    y = YC + RC_ * i
    rf = RFULL[r.cohort]
    rec = {"cohort": r.cohort}
    for j, c in enumerate(big_cl):
        m_ = np.flatnonzero(cl == c)
        ij = np.triu_indices(len(m_), 1)
        v = np.nanmedian(rf[np.ix_(m_, m_)][ij])
        tile(pg, x0 + cw * j, y, cw, RC_, PUOR(float((v + 1) / 2)), minus(f"{v:.2f}"))
        rec[f"C{c}"] = v
    for j, (name_, msk) in enumerate((("same", same), ("different", ~same))):
        v = np.nanmedian(rf[iu][msk])
        tile(pg, x0 + cw * (len(cr) + j) + 0.5, y, cw, RC_, PUOR(float((v + 1) / 2)), minus(f"{v:.2f}"))
        rec[name_] = v
    out.append(rec)
pg.src("p", pd.DataFrame(out))
pg.cbar(50, YC + HC + 1.3, 22.4, 1.2, PURPLE, 0, 60, [0, 30, 60], "% of variance")
pg.cbar(90, YC + HC + 1.3, 22, 1.2, PUOR, -1, 1, [-1, 0, 1], "median r")

# ================================================================ q  most-shared genes
pg.panel(123.8, 195.5, "q", "Most-shared genes")
ax = pg.ax(134.5, YC, 39.5, HC)
for i, (g, n_) in enumerate(top):
    y = YC + RC_ * i + RC_ / 2 - 0.07
    pg.txt(133.9, y, g, ha="right", style="italic")
    pg.txt(178.6, y, str(n_), ha="right")
    ax.plot([-0.5, len(holders) - 0.5], [y, y], color="#d9d9d9", lw=0.3, zorder=0)
    for j, s in enumerate(holders):
        if g in genes[s]:
            point(ax, j, y, FAM[fam_of[s]], ms=2.2, zorder=3)
pg.txt(178.6, YC - 1.2, "n", ha="right", weight="bold")
ax.set_ylim(YC + HC, YC); ax.set_xlim(-0.7, len(holders) - 0.3); ax.set_yticks([]); ax.set_xticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel(f"signatures containing them ({len(holders)})", labelpad=2)
pg.src("q", pd.DataFrame([{"gene": g, "n_signatures": n_, "signatures": ";".join(s for s in holders if g in genes[s])} for g, n_ in top]))

pg.save("Fig2")
print(f"called {int(CALL.sum())} | nodes {H.number_of_nodes()} edges {H.number_of_edges()} isolates {len(iso)} hubs {hubs} | meff median {med_meff} "
      f"| clusters {CLS.nunique()} sizes {sorted(csize.tolist(), reverse=True)} | r same {np.median(rv[same]):.3f} diff {np.median(rv[~same]):.3f} | pc1 {pc1:.3f} "
      f"| r>0.5 no gene {hi_r0}/{hi_r} | holders {len(holders)} | genes {len(cnt_g)} | blocks {blocks}")
