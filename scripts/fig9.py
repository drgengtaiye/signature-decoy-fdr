"""Figure 10, in the look of the AML multi-omic figures (pixl_nc.py). Round 15 parts B and C: the false-validation rate of
validation cohorts and what accounts for it (validation/PLAN_round15_literature_cohorts_controls.md, amendment 3; icb/v20_round15
from applications/v20_cohort_scores.R, v20_law.R and decoys/v20_cohorts_audit.py), with the law tested out of sample in round 18
(validation/PLAN_round18_law_out_of_sample.md; icb/v23_law_split from applications/v23_law_split.R and decoys/v23_law_split_audit.py).

Under the reading rule of round 18, H6c at K = 10 is not met, so the sheet keeps "predicted" for the out-of-sample panels alone.
Every P value printed comes from cohorts_h4.csv, cohorts_h5.csv, cohorts_h6.csv or law_split_h6.csv. The other coefficients are
descriptive (printed with n, no P). Nothing is resampled here. The script stops if a table disagrees with the quoted values."""
import glob
import re
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt, fm, TextPath
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

A = D + "/icb/v20_round15/"
V = D + "/signature-decoy-fdr/validation/"
F = pd.read_csv(A + "cohorts_fvr.csv")
H4 = pd.read_csv(A + "cohorts_h4.csv")
H5 = pd.read_csv(A + "cohorts_h5.csv").set_index("quantity")
H6 = pd.read_csv(A + "cohorts_h6.csv").set_index("K")
PAIRS = pd.read_csv(A + "cohorts_h5_pairs.csv")
LAW = {k: pd.read_csv(A + f"cohorts_law_K{k}.csv") for k in (10, 5, 20)}
CAND = pd.read_csv(V + "round15_cohorts/candidates_r15.csv")
PM = pd.read_csv(V + "published_models_r15.csv").set_index("model_id")
S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "cohort_scores_*.csv.gz"))], ignore_index=True)
S["kind"] = S["kind"].fillna("null")          # pandas reads the label "null" as missing, as in v20_cohorts_audit.py
S = S[S.scorable.astype(bool)].copy()
S["validated"] = S.validated.astype(float)
L = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "law_predictions*.csv.gz"))]).drop_duplicates(["cohort", "id", "K"])
OBS = S[S.kind == "null"][["cohort", "id", "validated"]]

# ---- the values the text quotes: stop rather than choose
Qs = H4[H4.type == "all (sum)"].iloc[0]
SP = H4[H4.type.str.startswith("Spearman")].iloc[0]
for got, want, nd in ((Qs.Q, 227.6, 1), (Qs.df, 24, 0), (SP.Q, 0.24, 2), (SP.p, 0.08, 2), (100 * H5.estimate["replication rate, rho_A >= 0.95"], 30.5, 1),
                      (100 * H5.estimate["replication rate, rho_A < 0.95"], 21.9, 1), (100 * H5.estimate["difference"], 8.6, 1),
                      (100 * H5.lo["difference"], 3.6, 1), (100 * H5.hi["difference"], 13.6, 1), (H5.estimate["H5 one-sided P (difference <= 0)"], 0.0007, 4),
                      (H5.estimate["pairs"], 3846, 0), (H5.estimate["models"], 322, 0), (100 * H5.estimate["replication rate in B of null models validated in A"], 37, 0),
                      (100 * H5.estimate["replication rate in B of co-expressed decoys validated in A"], 18.5, 1),
                      (H6.spearman[10], 0.89, 2), (H6.spearman[5], 0.75, 2), (H6.spearman[20], 0.85, 2), (H6.mae_law[10], 0.073, 3), (H6.mae_events[10], 0.169, 3),
                      (H6.mae_diff[10], -0.096, 3), (H6.mae_diff_lo[10], -0.143, 3), (H6.mae_diff_hi[10], -0.052, 3), (H6.auc_null_models[10], 0.85, 2)):
    assert round(float(got), nd) == want, (got, want)
# round 18, rounded from the table
A18 = D + "/icb/v23_law_split/"
HS8 = pd.read_csv(A18 + "law_split_h6.csv").set_index("K")
OOS = {k: pd.read_csv(A18 + f"law_split_K{k}.csv") for k in (10, 5, 20)}
SPLIT = pd.read_csv(A18 + "law_split.csv.gz")
for got, want, nd in ((HS8.spearman[10], 0.79, 2), (HS8.mae_law[10], 0.091, 3), (HS8.mae_events[10], 0.110, 3), (HS8.mae_diff[10], -0.019, 3), (HS8.mae_lo[10], -0.063, 3), (HS8.mae_hi[10], 0.025, 3),
                      (HS8.p_law_not_better[10], 0.20, 2), (HS8.spearman[5], 0.74, 2), (HS8.mae_law[5], 0.072, 3), (HS8.mae_events[5], 0.110, 3), (HS8.mae_diff[5], -0.038, 3), (HS8.mae_lo[5], -0.065, 3),
                      (HS8.mae_hi[5], -0.010, 3), (HS8.p_law_not_better[5], 0.0045, 4), (HS8.spearman[20], 0.59, 2), (HS8.cohorts[20], 28, 0), (HS8.mae_law[20], 0.182, 3), (HS8.spearman_events[10], 0.04, 2),
                      (HS8.auc_pooled[10], 0.68, 2), (HS8.auc_within_median[10], 0.53, 2), (100 * HS8.pred_mean[10], 21.5, 1), (100 * HS8.obs_mean[10], 15.9, 1)):
    assert round(float(got), nd) == want, (got, want)
assert HS8.p_perm[10] == 0 and HS8.p_perm[5] == 0 and HS8.mae_hi[10] > 0          # H6c at K = 10 is not met
for k_ in (10, 5, 20):
    assert abs(stats.spearmanr(OOS[k_].pred, OOS[k_].obs).statistic - HS8.spearman[k_]) < 1e-9 and abs(np.abs(OOS[k_].pred - OOS[k_].obs).mean() - HS8.mae_law[k_]) < 1e-9
# within-cohort AUC in sample, as decoys/v23_law_split_audit.py does for K = 10
AUCW = {}
for k_ in (10, 5, 20):
    _y = L[L.K == k_].merge(OBS, on=["cohort", "id"])
    _w = [roc_auc_score(g.validated, g.p_pred) for _, g in _y.groupby("cohort") if 0 < g.validated.sum() < len(g)]
    AUCW[k_] = (float(np.median(_w)), len(_w))
assert round(AUCW[10][0], 2) == 0.80
assert f"{Qs.p:.0e}" == "4e-35" and len(F) == 43 and F.independent_set.sum() == 35 and (H6.cohorts == 35).all() and H6.p_perm_one_sided[10] == 0

# ---- rows of the lookup table
F["acc"] = F.cohort.str.split("_", n=1).str[1]
SIB = {"GSE84433": "GSE84437", "GSE17536": "GSE17538", "GSE14333": "GSE17538", "GSE42127": "GSE41271", "GSE11969": "GSE13213",
       "GSE21050": "GSE71118", "GSE40967": "GSE39582", "GSE103091": "GSE58812"}     # amendment 3, decision 6
assert set(F.acc[~F.independent_set]) == set(SIB)
F["same_patients_as"] = F.acc.map(SIB)
F["twin"] = F.note.astype(str).str.startswith("rate of")                            # same patients and data: no rate of its own
F["rank"] = F.fvr.rank(method="min", ascending=False).astype(int)
_c = CAND.set_index(["type", "cohort"])
def _named(r):
    """papers naming the cohort. OV-AU is also the 'ICGC' ovarian cohort (amendment 3); the two sets of records do not overlap."""
    if r.cohort == "OV_OV-AU":
        return _c.records_naming[("OV", "OV-AU")] + _c.records_naming[("OV", "ICGC")]
    if r.cohort == "CESC_HTMCP-CC":
        return _c.records_naming[("CESC", "CGCI-HTMCP-CC")]
    return _c.records_naming[(r.type, r.acc)]
F["records_naming"] = F.apply(_named, axis=1).astype(int)
F["records_training_tcga_model"] = F.type.map(CAND.groupby("type").records_training_tcga_model.first())

# ---- random-gene decoys per cohort: the column the audit computes for co-expressed decoys, on the rand class
GROUPS = [{"GSE84437", "GSE84433", "GSE84426", "GSE147163"}, {"GSE39582", "GSE40967"},
          {"GSE17536", "GSE17537", "GSE17538", "GSE29621", "GSE14333"}, {"GSE58812", "GSE103091", "GSE83937"},
          {"GSE20711", "GSE20713"}, {"GSE96058", "GSE81540"}, {"GSE157011", "GSE157009", "GSE157010"},
          {"GSE41271", "GSE42127"}, {"GSE13213", "GSE11969"}, {"GSE26939", "GSE36471"}, {"GSE21050", "GSE71118"},
          {"OV-AU", "ICGC"}, {"HTMCP-CC", "CGCI-HTMCP-CC", "HTMCP", "CGCI"}, {"METABRIC"}]    # as v20_cohorts_audit.py


def paper_names(model, acc):
    names = next((g for g in GROUPS if acc in g), {acc})
    txt = str(PM.at[model, "validation_cohorts"]) if model in PM.index else ""
    return any(re.search(r"(?<![A-Za-z0-9])" + re.escape(n) + r"(?![A-Za-z0-9])", txt) for n in names)


drate = {}
for cid, acc in zip(F.cohort[~F.twin], F.acc[~F.twin]):
    x = S[S.cohort == cid]
    tg = x[x.kind == "target"]
    un = [m for m in tg.id if not paper_names(m, acc)]
    d = x[(x.kind == "decoy") & x.model.isin(un)]
    drate[cid] = (d[d["class"] == "coex"].validated.mean(), d[d["class"] == "rand"].validated.mean(), (d["class"] == "rand").sum())
F["rand_rate_not_naming"] = [drate[n.split("rate of ")[1] if t else c][1] for c, t, n in zip(F.cohort, F.twin, F.note.astype(str))]
_chk = F[~F.twin].assign(c2=lambda d: [drate[c][0] for c in d.cohort])
assert (np.abs(_chk.c2 - _chk.decoy_rate_not_naming) < 1e-9).all(), "co-expressed decoy rate does not reproduce the table"
DIRECT = F[~F.twin]                 # the 41 cohorts with a rate of their own
IND = F[F.independent_set]          # the 35 cohorts of the independent set


# ---- helpers
def tw(s, size=5.0):
    return TextPath((0, 0), s, prop=fm.FontProperties(family="Helvetica", size=size)).get_extents().width * 25.4 / 72


def rho(a, b):
    return stats.spearmanr(a, b).statistic


def pow10(pg, xr, y, p, ha="right", color=INK):
    """a P value as in the other sheets, without a superscript: two figures down to 0.001, below that p_small
    (0.0004, or < 0.0001). xr is the right end (or left, ha='left')."""
    if p >= 0.001:
        s = f"{p:.2g}" if p < 0.1 else f"{p:.2f}"
    else:
        op, v = p_small(p)
        s = f"< {v}" if op == "<" else v
    pg.txt(xr, y, s, ha=ha, color=color)


SQRT = ("function", {"functions": (np.sqrt, np.square)})


def place(pg, box, lim, xy, labs, scale=("linear", "linear"), rings=(1.1, 1.5, 2.0, 2.8, 3.6, 4.6, 5.8, 7.2, 9, 11, 13.5, 16), tries=60, seed=1, avoid=(), others=(), pad=0.3):
    """direct labels in mm: each takes the nearest free spot around its point, with a hairline when it had to move away.
    Layout only: seeded so that the sheet rebuilds identically; no data are touched."""
    x, y, w, h = box
    tf = {"linear": lambda v: v, "log": np.log10, "sqrt": np.sqrt}
    fr = lambda v, l, sc: (tf[sc](v) - tf[sc](l[0])) / (tf[sc](l[1]) - tf[sc](l[0]))  # noqa: E731
    mm = lambda q: np.array([[x + w * fr(a, lim[0], scale[0]), y + h * (1 - fr(b, lim[1], scale[1]))] for a, b in q]).reshape(-1, 2)  # noqa: E731
    P, OB = mm(xy), mm(others)
    OB = np.vstack([P, OB])
    wid = [tw(s_) for s_ in labs]
    TH = 1.75
    angs = np.deg2rad([0, 180, 35, 145, -35, -145, 90, -90, 18, 162, -18, -162, 60, 120, -60, -120])
    rng = np.random.default_rng(seed)
    best, bestc = None, 1e9
    dens = [-(np.hypot(*(P - p).T) < 8).sum() for p in P]
    for t in range(tries):
        order = np.argsort(dens) if t == 0 else rng.permutation(len(P))
        boxes, res, cost, segs = list(avoid), {}, 0.0, []
        for i in order:
            done = False
            for r in rings:
                for a_ in angs:
                    cx = P[i, 0] + np.cos(a_) * r + np.sign(np.round(np.cos(a_), 6)) * wid[i] / 2
                    cy = P[i, 1] - np.sin(a_) * r
                    b_ = (cx - wid[i] / 2 - pad, cy - TH / 2, cx + wid[i] / 2 + pad, cy + TH / 2)
                    if b_[0] < x + 0.3 or b_[2] > x + w + 1.5 or b_[1] < y - 1 or b_[3] > y + h - 0.2:
                        continue
                    if any(b_[0] < q[2] and q[0] < b_[2] and b_[1] < q[3] and q[1] < b_[3] for q in boxes):
                        continue
                    if ((OB[:, 0] > b_[0] - 0.55) & (OB[:, 0] < b_[2] + 0.55) & (OB[:, 1] > b_[1] - 0.55) & (OB[:, 1] < b_[3] + 0.55)).any():
                        continue
                    if any(b_[0] < u < b_[2] and b_[1] < v < b_[3] for u, v in segs):
                        continue
                    ex = min(max(P[i, 0], b_[0]), b_[2]); ey = min(max(P[i, 1], b_[1]), b_[3])
                    line = np.linspace(0.05, 0.95, 12)
                    sample = [(P[i, 0] + (ex - P[i, 0]) * f, P[i, 1] + (ey - P[i, 1]) * f) for f in line] if r > 1.6 else []
                    if any(q[0] < u < q[2] and q[1] < v < q[3] for u, v in sample for q in boxes):
                        continue
                    if r > 1.6 and any(np.hypot(u - o[0], v - o[1]) < 0.4 for u, v in sample[2:] for o in OB):
                        continue
                    boxes.append(b_); segs += sample; res[i] = (cx, cy, r, b_); cost += r
                    done = True
                    break
                if done:
                    break
            if not done:
                cost += 100
        if cost < bestc:
            best, bestc = res, cost
    for i, (cx, cy, r, b_) in best.items():
        pg.txt(cx, cy, labs[i], ha="center", bbox={"fc": "white", "ec": "none", "pad": 0.15}, zorder=4)
        if r > 1.6:
            ex = min(max(P[i, 0], b_[0]), b_[2]); ey = min(max(P[i, 1], b_[1]), b_[3])
            pg.line(P[i, 0], P[i, 1], ex, ey, color=OTHER, lw=0.3)
    missing = [labs[i] for i in range(len(P)) if i not in best]
    if missing:
        print("UNPLACED", missing)
    return [best[i][3] for i in best]


pg = Page(240)
HAIR = "#d4d4d4"

# ================================================================== a  the lookup table
ORDER = F.groupby("type").fvr.max().sort_values(ascending=False).index
Y0, CH, GAP = 17.6, 2.3, 0.9
rows_, yy = [], Y0
for k, tp in enumerate(ORDER):
    g = F[F.type == tp].sort_values(["fvr", "patients"], ascending=[False, False])
    if k:
        yy += GAP
    for j, (_, r) in enumerate(g.iterrows()):
        rows_.append((tp, r, yy + CH / 2, j == 0))
        yy += CH
Y1 = yy
pg.panel(1, 3, "a", "False-validation rate (FVR) of every cohort")
X_ACC, X_RANK, X_END, X_PAT, X_EV, X_NUL, XD0, XD1, X_FVR, X_SIB = 8.0, 22.4, 25.8, 34.8, 41.0, 47.9, 49.8, 71.0, 76.0, 77.2
HY = Y0 - 0.45
for x_, s_, ha_ in ((1, "type", "left"), (X_ACC, "cohort", "left"), (X_RANK, "rank", "right"), (X_END, "end", "center"), (X_PAT, "patients", "right"), (X_EV, "events", "right"),
                    (X_NUL, "null\nmodels", "right"), (X_FVR, "FVR\n(%)", "right"), (X_SIB, "same\npatients as", "left")):
    pg.txt(x_, HY, s_, ha=ha_, va="bottom", weight="bold", linespacing=1.0)
axa = pg.ax(XD0, Y0, XD1 - XD0, Y1 - Y0)
axa.set_ylim(Y1, Y0); axa.set_xlim(0, 100); axa.set_yticks([]); axa.spines["left"].set_visible(False)
axa.set_xticks([0, 25, 50, 75, 100]); axa.set_xlabel("FVR (%)")
zero(axa, 2.5)
for tp, r, y_, first in rows_:
    if first:
        pg.txt(1, y_, tp, weight="bold")
        if tp != ORDER[0]:
            pg.line(1, y_ - CH / 2 - GAP / 2, 144.6, y_ - CH / 2 - GAP / 2, color=HAIR, lw=0.25)
    out_ = not r.independent_set
    pg.txt(X_ACC, y_, r.acc, color=MUTED if out_ else INK)
    pg.txt(X_RANK, y_, str(r["rank"]), ha="right")
    pg.txt(X_END, y_, r.endpoint, ha="center")
    pg.txt(X_PAT, y_, f"{r.patients:,.0f}", ha="right")
    pg.txt(X_EV, y_, f"{r.events:,.0f}", ha="right")
    pg.txt(X_NUL, y_, str(r.null_models), ha="right")
    pg.txt(X_FVR, y_, f"{100 * r.fvr:.1f}", ha="right")
    if out_:
        pg.txt(X_SIB, y_, r.same_patients_as, color=MUTED)
    point(axa, 100 * r.fvr, y_, NULLC if out_ else INK, ms=3.0)
key(pg, 1, 9.0, [(INK, "independent set (35)", "o"), (NULLC, "outside the set (8)", "o"), (INK, "2.5%, β = 0", "dot")])
pg.src("a", pd.DataFrame([{**{k: r[k] for k in ("type", "cohort", "rank", "endpoint", "patients", "events", "null_models", "fvr", "independent_set", "same_patients_as", "note")}}
                          for _, r, _, _ in rows_]))

# ================================================================== b  rates per cohort, aligned with a
XB = 86.9
pg.panel(XB, 3, "b", "Rates per cohort (%)")
CW, NW, GG = 4.5, 3.6, 0.6
XBC = {}
x_ = XB + 0.4
for c_, w_ in (("fvr_coex", CW), ("fvr_rand", CW), ("g1", GG), ("model_rate_not_naming", CW), ("models_not_naming", NW), ("model_rate_naming", CW), ("models_naming", NW), ("g2", GG),
               ("decoy_rate_not_naming", CW), ("rand_rate_not_naming", CW)):
    XBC[c_] = (x_, w_); x_ += w_
XBEND = x_
for c_, s_ in (("fvr_coex", "co-expr."), ("fvr_rand", "random"), ("model_rate_not_naming", "not named"), ("models_not_naming", "n"), ("model_rate_naming", "named"), ("models_naming", "n"),
               ("decoy_rate_not_naming", "co-expr."), ("rand_rate_not_naming", "random")):
    pg.txt(XBC[c_][0] + XBC[c_][1] / 2, Y0 - 0.5, s_, rotation=90, ha="center", va="bottom")
for a_, b_, s_ in (("fvr_coex", "fvr_rand", "null models"), ("model_rate_not_naming", "models_naming", "published models"), ("decoy_rate_not_naming", "rand_rate_not_naming", "decoys")):
    x0_, x1_ = XBC[a_][0], XBC[b_][0] + XBC[b_][1] - 0.15
    pg.txt((x0_ + x1_) / 2, 6.6, s_, ha="center", weight="bold")
    pg.line(x0_, 7.8, x1_, 7.8, color=INK, lw=0.4)
for tp, r, y_, first in rows_:
    for c_ in ("fvr_coex", "fvr_rand", "model_rate_not_naming", "model_rate_naming", "decoy_rate_not_naming", "rand_rate_not_naming"):
        x0_, w_ = XBC[c_]
        v_ = r[c_]
        if pd.isna(v_):
            tile(pg, x0_, y_ - CH / 2, w_, CH, EMPTY)
        else:
            tile(pg, x0_, y_ - CH / 2, w_, CH, BLUE(float(v_)), f"{100 * v_:.0f}")
    for c_ in ("models_not_naming", "models_naming"):
        x0_, w_ = XBC[c_]
        pg.txt(x0_ + w_ - 0.5, y_, f"{r[c_]:.0f}", ha="right")
pg.cbar(XB + 0.4, Y1 + 2.0, 19, 1.4, BLUE, 0, 100, [0, 50, 100], "rate (%)")
pg.src("b", pd.DataFrame([{k: r[k] for k in ("type", "cohort", "fvr_coex", "fvr_rand", "models_not_naming", "model_rate_not_naming", "models_naming", "model_rate_naming",
                                             "decoy_rate_not_naming", "rand_rate_not_naming")} for _, r, _, _ in rows_]))

# ================================================================== c  papers naming it, aligned with a
XC0, XC1 = XBEND + 1.3, XBEND + 15.5
pg.panel(XBEND + 0.4, 3, "c", "Papers")
pg.txt(XC1 + 5.0, HY, "n", ha="right", va="bottom", weight="bold")
axc = pg.ax(XC0, Y0, XC1 - XC0, Y1 - Y0)
axc.set_ylim(Y1, Y0); axc.set_xscale("log"); axc.set_xlim(0.8, 110); axc.set_yticks([]); axc.spines["left"].set_visible(False)
axc.set_xticks([1, 10, 100], ["1", "10", "100"]); axc.minorticks_off(); axc.set_xlabel("papers naming it")
for tp, r, y_, first in rows_:
    point(axc, r.records_naming, y_, NULLC if not r.independent_set else INK, ms=3.0)
    pg.txt(XC1 + 5.0, y_, str(r.records_naming), ha="right")
pg.src("c", pd.DataFrame([{k: r[k] for k in ("type", "cohort", "records_naming", "records_training_tcga_model", "independent_set")} for _, r, _, _ in rows_]))

# ================================================================== right column: d events, e literature, f and g models and decoys
XR, XRA, WR = 146.2, 153.2, 24.4


def tlab(c):
    tp, a = c.split("_", 1)
    return a if not a.startswith("GSE") else f"{tp} {a[3:]}"


pg.panel(XR, 3, "d", "FVR and events")
box = (XRA, 8.6, WR, 20.0)
ax = pg.ax(*box)
for out_, d_ in ((True, DIRECT[~DIRECT.independent_set]), (False, IND)):
    point(ax, d_.events, 100 * d_.fvr, NULLC if out_ else INK, ms=2.4)
ax.set_xscale("log"); ax.set_xlim(20, 1500); ax.set_ylim(0, 100); ax.set_xticks([30, 100, 300, 1000], ["30", "100", "300", "1,000"]); ax.minorticks_off(); ax.set_yticks([0, 50, 100])
ax.set_xlabel("events"); ax.set_ylabel("FVR (%)")
pg.txt(XRA + 0.5, 8.9, f"Spearman {SP.Q:.2f}, P = {SP.p:.2f}, n = {int(SP.cohorts)}")
lb = pd.concat([IND.nlargest(2, "fvr"), IND.nlargest(1, "events")]).drop_duplicates("cohort")
place(pg, box, ((20, 1500), (0, 100)), list(zip(lb.events, 100 * lb.fvr)), [tlab(c) for c in lb.cohort], scale=("log", "linear"), avoid=[(XRA, 7.8, XRA + 25, 10.1)],
      others=list(zip(DIRECT.events, 100 * DIRECT.fvr)))
pg.src("d", DIRECT[["cohort", "type", "events", "patients", "fvr", "independent_set"]])

pg.panel(XR, 35.6, "e", "Papers and FVR")
box = (XRA, 41.0, WR, 20.0)
ax = pg.ax(*box)
for out_, d_ in ((True, F[~F.independent_set]), (False, IND)):
    point(ax, d_.records_naming, 100 * d_.fvr, NULLC if out_ else INK, ms=2.4)
ax.set_xscale("log"); ax.set_xlim(0.8, 110); ax.set_ylim(0, 100); ax.set_xticks([1, 10, 100], ["1", "10", "100"]); ax.minorticks_off(); ax.set_yticks([0, 50, 100])
ax.set_xlabel("papers naming it"); ax.set_ylabel("FVR (%)")
rp = rho(IND.records_naming, IND.fvr)
pg.txt(XRA + 0.5, 41.3, f"Spearman {rp:.2f}, n = {len(IND)}")
lb = IND.nlargest(2, "fvr")
place(pg, box, ((0.8, 110), (0, 100)), list(zip(lb.records_naming, 100 * lb.fvr)), [tlab(c) for c in lb.cohort], scale=("log", "linear"), avoid=[(XRA, 40.2, XRA + 20, 42.5)],
      others=list(zip(F.records_naming, 100 * F.fvr)))
pg.src("e", F[["cohort", "type", "records_naming", "fvr", "independent_set"]])

pg.panel(XR, 68.0, "f", "Paper did not name it")
ax = pg.ax(XRA, 76.0, WR, 20.2)
ax.plot([0, 100], [0, 100], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * IND.fvr, 100 * IND.rand_rate_not_naming, RAND, m="D", ms=1.9)
point(ax, 100 * IND.fvr, 100 * IND.decoy_rate_not_naming, COEX, m="s", ms=2.0)
point(ax, 100 * IND.fvr, 100 * IND.model_rate_not_naming, SIG, ms=2.4, zorder=5)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks([0, 50, 100]); ax.set_yticks([0, 50, 100])
ax.set_xlabel("FVR (%)"); ax.set_ylabel("validated (%)")
key(pg, XR + 0.2, 71.5, [(SIG, f"models {rho(IND.fvr, IND.model_rate_not_naming):.2f}", "o"), (COEX, f"co-expr. {rho(IND.fvr, IND.decoy_rate_not_naming):.2f}", "s")], gap=1.6)
key(pg, XR + 0.2, 73.8, [(RAND, f"random {rho(IND.fvr, IND.rand_rate_not_naming):.2f}", "D")], gap=1.6)
pg.txt(160.6, 73.8, f"Spearman, n = {len(IND)}")
pg.src("f", IND[["cohort", "type", "fvr", "models_not_naming", "model_rate_not_naming", "decoy_rate_not_naming", "rand_rate_not_naming"]])

pg.panel(XR, 103.2, "g", "Paper named it")
d_ = IND[IND.model_rate_naming.notna() & (IND.models_naming > 0)]
ax = pg.ax(XRA, 108.6, WR, 17.4)
ax.plot([0, 100], [0, 100], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * d_.fvr, 100 * d_.model_rate_naming, SIG, ms=2.4)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks([0, 50, 100]); ax.set_yticks([0, 50, 100])
ax.set_xlabel("FVR (%)"); ax.set_ylabel("validated (%)")
pg.txt(XRA + WR, 124.7, f"Spearman {rho(d_.fvr, d_.model_rate_naming):.2f}, n = {len(d_)}", ha="right")
pg.src("g", d_[["cohort", "type", "fvr", "models_naming", "model_rate_naming"]])

# ================================================================== h, i  H4 Cochran's Q, shared rows
YA = 134.0
pg.panel(1, YA, "h", "Cochran's Q, cohorts within a type")
Q_ = H4[~H4.type.str.startswith(("all", "Spearman"))].sort_values("Q", ascending=False).reset_index(drop=True)
YQ, CQ = YA + 8.6, 2.6
XQ0, XQ1 = 21.0, 38.0
yq = lambda i: YQ + CQ * i + CQ / 2   # noqa: E731
HQ = YQ - 0.4
for x_, s_, ha_ in ((1, "type", "left"), (12.6, "cohorts", "right"), (19.6, "null\nmodels", "right"), (44.3, "Q", "right"), (48.2, "df", "right"), (59.4, "P", "right")):
    pg.txt(x_, HQ, s_, ha=ha_, va="bottom", weight="bold", linespacing=1.0)
axq = pg.ax(XQ0, YQ, XQ1 - XQ0, CQ * 9 + 0.6)
axq.set_ylim(YQ + CQ * 9 + 0.6, YQ); axq.set_xscale("log"); axq.set_xlim(0.8, 300); axq.set_yticks([]); axq.spines["left"].set_visible(False)
axq.set_xticks([1, 10, 100], ["1", "10", "100"]); axq.minorticks_off(); axq.set_xlabel("Q")
for i, r in Q_.iterrows():
    y_ = yq(i)
    pg.txt(1, y_, r.type)
    pg.txt(12.6, y_, f"{r.cohorts:.0f}", ha="right")
    pg.txt(19.6, y_, f"{r.null_models:.0f}", ha="right")
    axq.plot([r.df, r.df], [y_ - 0.8, y_ + 0.8], color=MUTED, lw=0.6, zorder=2, solid_capstyle="butt")
    point(axq, r.Q, y_, INK, ms=3.0, filled=bool(r.p < 0.05))
    pg.txt(44.3, y_, f"{r.Q:.1f}", ha="right")
    pg.txt(48.2, y_, f"{r.df:.0f}", ha="right")
    pow10(pg, 59.4, y_, r.p)
y_ = yq(8) + 0.6
pg.line(1, y_ - CQ / 2, 59.4, y_ - CQ / 2, color=INK, lw=0.3)
pg.txt(1, y_, "sum", weight="bold")
axq.plot([Qs.df, Qs.df], [y_ - 0.8, y_ + 0.8], color=MUTED, lw=0.6, zorder=2, solid_capstyle="butt")
point(axq, Qs.Q, y_, INK, ms=3.4)
pg.txt(44.3, y_, f"{Qs.Q:.1f}", ha="right", weight="bold")
pg.txt(48.2, y_, f"{Qs.df:.0f}", ha="right", weight="bold")
pow10(pg, 59.4, y_, Qs.p)
key(pg, 21.0, YA + 4.6, [(INK, "P < 0.05", "o"), (INK, "P ≥ 0.05", "open-o")])
pg.line(44.6, YA + 3.8, 44.6, YA + 5.4, color=MUTED, lw=0.6)
pg.txt(45.4, YA + 4.6, "df")
pg.src("h", H4[~H4.type.str.startswith("Spearman")][["type", "cohorts", "null_models", "Q", "df", "p"]])

XH = 61.6
pg.panel(XH, YA, "i", "FVR, shared models")
XH0, XH1 = XH + 1.2, XH + 18.4
axh = pg.ax(XH0, YQ, XH1 - XH0, CQ * 9 + 0.6)
axh.set_ylim(YQ + CQ * 9 + 0.6, YQ); axh.set_xlim(0, 100); axh.set_yticks([]); axh.spines["left"].set_visible(False)
axh.set_xticks([0, 50, 100]); axh.set_xlabel("FVR (%)")
zero(axh, 2.5)
pg.txt(XH1 + 7.6, HQ, "range", ha="right", va="bottom", weight="bold")
orows = []
for i, r in Q_.iterrows():
    v_ = {a.split("=")[0]: 100 * float(a.split("=")[1]) for a in r.rates.split(";")}
    vals = np.array(list(v_.values()))
    y_ = yq(i)
    axh.plot([vals.min(), vals.max()], [y_, y_], color=OTHER, lw=0.5, zorder=2, solid_capstyle="butt")
    point(axh, vals, np.full(len(vals), y_), INK, ms=2.3)
    pg.txt(XH1 + 7.6, y_, f"{vals.min():.0f}–{vals.max():.0f}", ha="right")
    orows += [{"type": r.type, "cohort": k, "rate_pct_common_null_models": v, "common_null_models": r.null_models} for k, v in v_.items()]
pg.src("i", pd.DataFrame(orows))

# ================================================================== j  H5
XJ = 90.4
pg.panel(XJ, YA, "j", "Replication in B, by ρ in A")
pg.txt(XJ, YA + 4.2, f"{H5.estimate['pairs']:,.0f} pairs, {H5.estimate['models']:.0f} models; P = {H5.estimate['H5 one-sided P (difference <= 0)']:.4f}")
rl = [("models, ρ ≥ 0.95", H5.estimate["replication rate, rho_A >= 0.95"], H5.estimate["pairs with rho_A >= 0.95"], SIG, "o", True),
      ("models, ρ < 0.95", H5.estimate["replication rate, rho_A < 0.95"], H5.estimate["pairs"] - H5.estimate["pairs with rho_A >= 0.95"], SIG, "o", False),
      ("null models", H5.estimate["replication rate in B of null models validated in A"], H5.n["replication rate in B of null models validated in A"], INK, "o", True),
      ("co-expr. decoys", H5.estimate["replication rate in B of co-expressed decoys validated in A"], H5.n["replication rate in B of co-expressed decoys validated in A"], COEX, "s", True)]
YL, CL = YA + 9.2, 3.0
XL0, XL1, XLV, XLN = XJ + 15.6, XJ + 32.6, XJ + 37.8, XJ + 45.2
pg.txt(XLV, YL - 0.4, "%", ha="right", va="bottom", weight="bold")
pg.txt(XLN, YL - 0.4, "pairs", ha="right", va="bottom", weight="bold")
axl = pg.ax(XL0, YL, XL1 - XL0, CL * 5 + 0.3)
axl.set_ylim(YL + CL * 5 + 0.3, YL); axl.set_xlim(0, 45); axl.set_yticks([]); axl.spines["left"].set_visible(False)
axl.set_xticks([0, 20, 40]); axl.set_xlabel("replication in B (%)")
zero(axl, 0)
lrows = []
for i, (s_, v_, n_, c_, m_, f_) in enumerate(rl):
    y_ = YL + CL * i + CL / 2
    pg.txt(XJ, y_, s_)
    point(axl, 100 * v_, y_, c_, m=m_, ms=3.4, filled=f_)
    pg.txt(XLV, y_, f"{100 * v_:.1f}", ha="right")
    pg.txt(XLN, y_, f"{n_:,.0f}", ha="right")
    lrows.append({"row": s_, "replication_rate": v_, "pairs": n_})
y_ = YL + CL * 4 + CL / 2
pg.txt(XJ, y_, "difference, models")
d_, lo_, hi_ = (100 * H5.estimate["difference"], 100 * H5.lo["difference"], 100 * H5.hi["difference"])
ci(axl, lo_, hi_, y_, INK, lw=0.7)
point(axl, d_, y_, INK, ms=3.0)
pg.txt(XLN, y_, f"{d_:.1f} ({lo_:.1f}–{hi_:.1f})", ha="right")
lrows.append({"row": "difference (points), 95% interval", "replication_rate": d_ / 100, "lo": lo_ / 100, "hi": hi_ / 100, "P_one_sided": H5.estimate["H5 one-sided P (difference <= 0)"]})
pg.src("j", pd.DataFrame(lrows))

# ================================================================== k  H5 by type (same weights)
XK = 138.2
pg.panel(XK, YA, "k", "Replication, by type (%)")


def wr(g):
    return np.average(g.replicated, weights=g.w) if len(g) else np.nan


BT = pd.DataFrame([{"type": t, "models": g.model.nunique(), "pairs": len(g), "hi": wr(g[g.pass_A]), "lo": wr(g[~g.pass_A]), "pairs_hi": int(g.pass_A.sum()), "pairs_lo": int((~g.pass_A).sum())}
                   for t, g in PAIRS.groupby("type")]).sort_values("pairs", ascending=False).reset_index(drop=True)
YT, CT = YA + 7.4, 3.2
for x_, s_, ha_ in ((XK, "type", "left"), (XK + 12.6, "models", "right"), (XK + 18.6, "pairs", "right"), (XK + 24.2, "ρ ≥ 0.95", "center"), (XK + 34.6, "ρ < 0.95", "center")):
    pg.txt(x_, YT - 0.4, s_, ha=ha_, va="bottom", weight="bold")
for i, r in BT.iterrows():
    y_ = YT + CT * i
    pg.txt(XK, y_ + CT / 2, r.type)
    pg.txt(XK + 12.6, y_ + CT / 2, str(r.models), ha="right")
    pg.txt(XK + 18.6, y_ + CT / 2, f"{r.pairs:,}", ha="right")
    for x_, c_ in ((XK + 19.2, "hi"), (XK + 29.6, "lo")):
        if pd.isna(r[c_]):
            tile(pg, x_, y_, 10.4, CT, EMPTY)
        else:
            tile(pg, x_, y_, 10.4, CT, BLUE(float(r[c_])), f"{100 * r[c_]:.0f}")
pg.src("k", BT.rename(columns={"hi": "replication_rho_ge_095", "lo": "replication_rho_lt_095", "pairs_hi": "pairs_rho_ge_095", "pairs_lo": "pairs_rho_lt_095"}))

# ================================================================== l  the law in sample, K = 10
YB = 174.0
TK = [0, 5, 10, 25, 50, 100]
SMALLTK = [0, 25, 100]
pg.panel(1, YB, "l", "Law and observed FVR, K = 10")
XS, WS, YS, HS = 9.0, 28.0, YB + 6.0, 23.0
ax = pg.ax(XS, YS, WS, HS)
C = LAW[10].merge(F[["cohort", "type", "fvr"]], on="cohort")
ax.plot([0, 100], [0, 100], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * C.pred, 100 * C.obs, INK, ms=2.2)
ax.set_xscale(SQRT[0], **SQRT[1]); ax.set_yscale(SQRT[0], **SQRT[1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks([0, 10, 25, 50, 100], ["0", "10", "25", "50", "100"]); ax.set_yticks([0, 10, 25, 50, 100], ["0", "10", "25", "50", "100"]); ax.minorticks_off()
ax.set_xlabel("law's FVR, in sample (%)"); ax.set_ylabel("observed FVR (%)")
pg.txt(XS + 0.6, YS + 0.3, f"Spearman {H6.spearman[10]:.2f}")
pg.txt(XS + 0.6, YS + 2.4, "P < 0.0001")
pg.txt(XS + 0.6, YS + 4.5, f"n = {int(H6.cohorts[10])}")
top = C.nlargest(1, "obs")
place(pg, (XS, YS, WS, HS), ((0, 100), (0, 100)), list(zip(100 * top.pred, 100 * top.obs)), [tlab(c) for c in top.cohort], scale=("sqrt", "sqrt"), avoid=[(XS, YS - 1, XS + 13, YS + 6)],
      others=list(zip(100 * C.pred, 100 * C.obs)), rings=(1.1, 1.5, 2.0, 2.8, 3.6, 4.6, 5.8, 7.2))
pg.src("l", C[["cohort", "type", "events", "n", "k", "pred", "obs", "pred_events"]].rename(columns={"pred": "law_in_sample", "pred_events": "events_only"}))

# ================================================================== m  K = 5 and K = 20 in sample, merged
YM = 210.2
pg.panel(1, YM, "m", "K = 5 and K = 20, in sample")
ax = pg.ax(XS, YM + 5.8, WS, 17.6)
ax.plot([0, 100], [0, 100], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * LAW[5].pred, 100 * LAW[5].obs, NULLC, ms=2.0, zorder=3)
point(ax, 100 * LAW[20].pred, 100 * LAW[20].obs, INK, ms=2.0, zorder=4)
ax.set_xscale(SQRT[0], **SQRT[1]); ax.set_yscale(SQRT[0], **SQRT[1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks([0, 10, 25, 50, 100], ["0", "10", "25", "50", "100"]); ax.set_yticks([0, 25, 100], ["0", "25", "100"]); ax.minorticks_off()
ax.set_xlabel("law's FVR, in sample (%)"); ax.set_ylabel("observed FVR (%)")
pg.txt(1, YM + 3.5, "Spearman")
key(pg, 9.0, YM + 3.5, [(NULLC, f"K = 5, {H6.spearman[5]:.2f}", "o"), (INK, f"K = 20, {H6.spearman[20]:.2f}", "o")], gap=1.4)
pg.src("m", pd.concat([LAW[k][["cohort", "events", "n", "k", "pred", "obs"]].assign(K=k) for k in (5, 20)]).rename(columns={"pred": "law_in_sample"}))

# ================================================================== n  out of sample, K = 10
pg.panel(41.0, YB, "n", "Predicted and observed FVR, half B, K = 10")
XN, WN, YNA, HN = 48.4, 46.0, YB + 6.0, 53.6
ax = pg.ax(XN, YNA, WN, HN)
C = OOS[10]
LIMO = 75
TKO = [0, 5, 10, 25, 50, 75]
ax.plot([0, LIMO], [0, LIMO], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * C.pred, 100 * C.obs, INK, ms=2.8)
ax.set_xscale(SQRT[0], **SQRT[1]); ax.set_yscale(SQRT[0], **SQRT[1])
ax.set_xlim(0, LIMO); ax.set_ylim(0, LIMO); ax.set_xticks(TKO, [str(t) for t in TKO]); ax.set_yticks(TKO, [str(t) for t in TKO]); ax.minorticks_off()
ax.set_xlabel("predicted FVR, half B (%)"); ax.set_ylabel("observed FVR, half B (%)")
tx = XN + 0.8
pg.txt(tx, YNA + 0.3, f"Spearman {HS8.spearman[10]:.2f}")
pg.txt(tx, YNA + 2.4, "P < 0.0001")
pg.txt(tx, YNA + 4.5, f"n = {int(HS8.cohorts[10])} cohorts")
pg.txt(tx, YNA + 6.6, f"events only, {HS8.spearman_events[10]:.2f}")
pg.txt(tx, YNA + 8.7, f"mean predicted {100 * HS8.pred_mean[10]:.1f}%")
pg.txt(tx, YNA + 10.8, f"mean observed {100 * HS8.obs_mean[10]:.1f}%")
place(pg, (XN, YNA, WN, HN), ((0, LIMO), (0, LIMO)), list(zip(100 * C.pred, 100 * C.obs)), [tlab(c) for c in C.cohort], scale=("sqrt", "sqrt"), avoid=[(XN, YNA - 1, XN + 20, YNA + 12)], tries=2500, seed=1, pad=0.12, rings=(1.1, 1.5, 2.0, 2.8, 3.6, 4.6, 5.8, 7.2, 9))
pg.src("n", C[["cohort", "type", "D_B", "splits", "n", "k", "pred", "obs", "pred_events"]])

# ================================================================== o  in sample against out of sample
XO = 96.6
pg.panel(XO, YB, "o", "Law in sample and out of sample")
XK_, XF_, XN_ = XO, XO + 4.2, XO + 18.9
# P and AUC columns 1.3 mm further right than in v19, to fit "< 0.0001" (no superscript) after the difference column
TL = [("sp", XO + 19.9, 7.9), ("law", XO + 27.8, 6.2), ("ev", XO + 34.0, 6.2), ("pool", XO + 70.5, 6.1), ("with", XO + 76.6, 6.1)]
XDF, XP = XO + 61.4, XO + 69.8
YH, CO_ = YB + 7.4, 2.72
pg.txt(XK_, YH - 0.4, "K", va="bottom", weight="bold")
pg.txt(XN_, YH - 0.4, "n", ha="right", va="bottom", weight="bold")
pg.txt(TL[0][1] + TL[0][2] / 2 - 0.08, YH - 0.4, "Spearman", ha="center", va="bottom", weight="bold")
pg.txt(TL[1][1] + TL[1][2] / 2 - 0.08, YH - 0.4, "law", ha="center", va="bottom", weight="bold")
pg.txt(TL[2][1] + TL[2][2] / 2 - 0.08, YH - 0.4, "events", ha="center", va="bottom", weight="bold")
pg.txt((XO + 40.6 + XDF) / 2, YH - 0.4, "difference (95% interval)", ha="center", va="bottom", weight="bold")
pg.txt(XP, YH - 0.4, "P", ha="right", va="bottom", weight="bold")
pg.txt(TL[3][1] + TL[3][2] / 2 - 0.08, YH - 0.4, "pooled", ha="center", va="bottom", weight="bold")
pg.txt(TL[4][1] + TL[4][2] / 2 - 0.08, YH - 0.4, "within", ha="center", va="bottom", weight="bold")
for x0_, x1_, s_ in ((TL[1][1], TL[2][1] + TL[2][2] - 0.15, "MAE"), (TL[3][1], TL[4][1] + TL[4][2] - 0.15, "AUC")):
    pg.txt((x0_ + x1_) / 2, YH - 3.9, s_, ha="center", weight="bold")
    pg.line(x0_, YH - 2.9, x1_, YH - 2.9, color=INK, lw=0.4)
trows = []
y_ = YH
for K in (10, 5, 20):
    s_in = dict(fit="in sample", n=H6.cohorts[K], sp=H6.spearman[K], law=H6.mae_law[K], ev=H6.mae_events[K], d=H6.mae_diff[K], lo=H6.mae_diff_lo[K], hi=H6.mae_diff_hi[K],
                p=H6.p_law_not_better[K], pool=H6.auc_null_models[K], within=AUCW[K][0])
    s_out = dict(fit="out of sample", n=HS8.cohorts[K], sp=HS8.spearman[K], law=HS8.mae_law[K], ev=HS8.mae_events[K], d=HS8.mae_diff[K], lo=HS8.mae_lo[K], hi=HS8.mae_hi[K],
                 p=HS8.p_law_not_better[K], pool=HS8.auc_pooled[K], within=HS8.auc_within_median[K])
    for j, q in enumerate((s_in, s_out)):
        yc = y_ + CO_ / 2
        if j == 0:
            pg.txt(XK_, yc, str(K), weight="bold")
        pg.txt(XF_, yc, q["fit"], weight="bold" if j else None)
        pg.txt(XN_, yc, f"{int(q['n'])}", ha="right")
        tile(pg, TL[0][1], y_, TL[0][2], CO_, BLUE(float(q["sp"])), f"{q['sp']:.2f}")
        tile(pg, TL[1][1], y_, TL[1][2], CO_, PURPLE(float(q["law"]) / 0.25), f"{q['law']:.3f}")
        tile(pg, TL[2][1], y_, TL[2][2], CO_, PURPLE(float(q["ev"]) / 0.25), f"{q['ev']:.3f}")
        pg.txt(XDF, yc, f"{minus(f'{q['d']:.3f}')} ({minus(f'{q['lo']:.3f}')} to {minus(f'{q['hi']:.3f}')})", ha="right")
        if q["p"] == 0:
            pg.txt(XP, yc, "< 0.0001", ha="right")
        else:
            pow10(pg, XP, yc, q["p"])
        tile(pg, TL[3][1], y_, TL[3][2], CO_, BLUE(float((q["pool"] - 0.5) / 0.5)), f"{q['pool']:.2f}")
        tile(pg, TL[4][1], y_, TL[4][2], CO_, BLUE(float((q["within"] - 0.5) / 0.5)), f"{q['within']:.2f}")
        trows.append({"K": K, **q})
        y_ += CO_
    y_ += 0.4
pg.src("o", pd.DataFrame(trows).rename(columns={"d": "mae_difference", "lo": "difference_lo", "hi": "difference_hi", "p": "p_law_not_better", "pool": "auc_pooled", "within": "auc_within_median",
                                                "sp": "spearman", "law": "mae_law", "ev": "mae_events_only"}))

# ================================================================== p  out of sample, K = 5
YC = 203.0
AXY, AXH = YC + 5.6, 25.0
pg.panel(XO, YC, "p", "K = 5, half B")
ax = pg.ax(XO + 7.0, AXY, 18.6, AXH)
C_ = OOS[5]
ax.plot([0, LIMO], [0, LIMO], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
point(ax, 100 * C_.pred, 100 * C_.obs, INK, ms=2.0)
ax.set_xscale(SQRT[0], **SQRT[1]); ax.set_yscale(SQRT[0], **SQRT[1])
ax.set_xlim(0, LIMO); ax.set_ylim(0, LIMO); ax.set_xticks([0, 10, 25, 75], ["0", "10", "25", "75"]); ax.set_yticks([0, 10, 25, 75], ["0", "10", "25", "75"]); ax.minorticks_off()
ax.set_xlabel("predicted, half B (%)"); ax.set_ylabel("observed, half B (%)")
pg.txt(XO + 7.4, AXY + 0.3, f"Spearman {HS8.spearman[5]:.2f}")
pg.txt(XO + 7.4, AXY + 2.4, f"n = {int(HS8.cohorts[5])}")
pg.src("p", C_[["cohort", "type", "D_B", "splits", "n", "k", "pred", "obs", "pred_events"]])

# ================================================================== q  ROC, K = 10, in sample and out of sample
XQ = 125.0
pg.panel(XQ, YC, "q", "ROC, K = 10")
ax = pg.ax(XQ + 7.0, AXY, 18.6, AXH)
ax.plot([0, 1], [0, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
L10 = L[L.K == 10].merge(OBS, on=["cohort", "id"])
R10 = SPLIT[SPLIT.K == 10]
assert abs(roc_auc_score(L10.validated, L10.p_pred) - H6.auc_null_models[10]) < 1e-9 and len(L10) == H6.pairs[10]
assert abs(roc_auc_score(R10.validated, R10.p_pred) - HS8.auc_pooled[10]) < 1e-9 and len(R10) == HS8.rows[10]
fi, ti, _ = roc_curve(L10.validated, L10.p_pred)
fo, to, _ = roc_curve(R10.validated, R10.p_pred)
ax.plot(fi, ti, color=NULLC, lw=0.9, zorder=2, solid_capstyle="butt")
ax.plot(fo, to, color=INK, lw=0.9, zorder=3, solid_capstyle="butt")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
ax.set_xlabel("false-positive rate"); ax.set_ylabel("true-positive rate")
pg.box(XQ + 12.0, AXY + AXH - 5.9, 13.4, 5.7, "white", z=1.5)
key(pg, XQ + 12.2, AXY + AXH - 4.4, [(NULLC, f"in sample, {H6.auc_null_models[10]:.2f}", "line")], sw=1.4)
key(pg, XQ + 12.2, AXY + AXH - 2.2, [(INK, f"half B, {HS8.auc_pooled[10]:.2f}", "line")], sw=1.4)
pg.src("q", pd.concat([pd.DataFrame({"fit": "in sample", "fpr": fi, "tpr": ti}), pd.DataFrame({"fit": "out of sample", "fpr": fo, "tpr": to})]))

# ================================================================== r  validated against the law's probability, K = 10, in sample and out of sample
XR2 = 153.0
pg.panel(XR2, YC, "r", "Law and validated")
ax = pg.ax(XR2 + 7.0, AXY, 17.6, AXH)
bins = np.linspace(0, 1, 11)
ax.plot([0, 100], [0, 100], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
brows = []
for name, X_, col, z_ in (("in sample", L10, NULLC, 3), ("out of sample", R10, INK, 4)):
    X_ = X_.assign(bin=np.clip(np.digitize(X_.p_pred, bins) - 1, 0, 9))
    B_ = X_.groupby("bin").agg(law=("p_pred", "mean"), share=("validated", "mean"), n=("validated", "size")).reset_index().assign(fit=name)
    ax.plot(100 * B_.law, 100 * B_.share, color=col if col != NULLC else OTHER, lw=0.4, zorder=2)
    point(ax, 100 * B_.law, 100 * B_.share, col, ms=2.2, zorder=z_)
    brows.append(B_)
    print(name, "bin n", B_.n.min(), B_.n.max(), "validated", int(X_.validated.sum()), "of", len(X_))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks([0, 50, 100]); ax.set_yticks([0, 50, 100])
ax.set_xlabel("law's probability (%)"); ax.set_ylabel("validated (%)")
key(pg, XR2 + 7.4, AXY + 0.4, [(NULLC, "in sample", "o")])
key(pg, XR2 + 7.4, AXY + 2.6, [(INK, "half B", "o")])
pg.src("r", pd.concat(brows).assign(bin_lo=lambda d: d.bin / 10, bin_hi=lambda d: (d.bin + 1) / 10))

pg.save("Fig10")
