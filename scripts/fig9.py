"""Figure 9 (NC v16): published prognostic risk models and the pipeline that builds them, against matched decoys
(validation/PLAN_round13_published_risk_models.md with amendments 1-2; icb/v19_audit from applications/v19_*.R and
decoys/v19_audit.py).

a  the screen: records found, screened, eligible, scored, with decoys, independent of the external cohort
b  validated share in TCGA (training) and in the independent cohort: models, co-expressed and random-gene decoys
c  the decoy estimate of the false discovery rate of "this model validates", by population and decoy class
d  each independent model among its 32 co-expressed decoys in the independent cohort (directed C-index)
e  rho, the share of a model's decoys whose C-index it exceeds
f  the pipeline run on theme pools and on their decoy pools, by cancer type
g  the pipeline by theme pool: validated theme units against the number their decoys predict
"""
import sys, glob, json
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
K.apply()
A = K.D + "/icb/v19_audit/"
V = K.D + "/signature-decoy-fdr/validation/"
U = pd.read_csv(A + "audit_models_units.csv")
U = U[U.n_coex.notna()]
R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "models_[A-Z]*.csv.gz"))])
DEC = R[R.kind == "decoy"]
AM = pd.read_csv(A + "audit_models.csv")
PU = pd.read_csv(A + "audit_pipeline_units.csv")
PA = pd.read_csv(A + "audit_pipeline.csv")
PD = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "pipeline_decoys_*.csv.gz"))])
PT = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "pipeline_targets_*.csv.gz"))])
SEARCH = json.load(open(V + "round13_search/search_summary.json"))
PM = pd.read_csv(V + "published_models.csv")
SCR = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(V + "round13_search/screen_*.csv"))])
IND = U[~U.uses_external.astype(bool)]
MOD, COEX, RAND = K.INK, K.TEAL, K.LGRAY


def am(pop, q, dec="coex"):
    x = AM[AM.population.str.startswith(pop) & (AM.quantity == q) & (AM.decoys == dec)]
    return x.iloc[0]


C = K.Canvas(236)

# ================================================================== a  the screen
ax = C.panel("a", 24, 12, 40, 38, head="From search to scored models")
steps = [("found by the search", sum(v["count"] for v in SEARCH.values())), ("screened", len(SCR)),
         ("eligible", int(SCR.eligible.astype(str).str.upper().eq("TRUE").sum())), ("scored", int(PM.scored.sum())),
         ("with decoys (≥ 3 genes)", len(U)), ("independent of the cohort", len(IND))]
for i, (lab_, n) in enumerate(steps):
    y = len(steps) - 1 - i
    ax.barh(y, np.log10(n), 0.62, color=K.tint(MOD, 0.25 + 0.12 * i), edgecolor="white", lw=0.6, zorder=3)
    ax.text(np.log10(n) + 0.05, y, f"{n:,}", va="center", ha="left", fontsize=K.T_SMALL, color=K.INK)
ax.set_yticks(range(len(steps))); ax.set_yticklabels([s for s, _ in steps][::-1], fontsize=K.T_SMALL)
ax.set_xlim(0, 3.9); ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["1", "10", "100", "1,000"])
ax.set_xlabel("records or models (log scale)"); K.frame(ax, left=False)

# ================================================================== b  validated share, training and independent
ax = C.panel("b", 84, 12, 40, 38, head="Validated, in TCGA and outside it")
rT = DEC[(DEC["class"] == "rand") & DEC.model.isin(IND.model)]
vals = {"TCGA (training)": (IND.tcga_validated.mean(), IND.d_tcga.mean(),
                            ((rT.tcga_p < 0.05) & (rT.tcga_worse == 1)).groupby(rT.model).mean().mean()),
        "independent cohort": (IND.validated.mean(), IND.d_coex.mean(), IND.d_rand.mean())}
for i, (k, v) in enumerate(vals.items()):
    for j, (val, col) in enumerate(zip(v, (MOD, COEX, RAND))):
        ax.bar(i + (j - 1) * 0.26, 100 * val, 0.24, color=col, edgecolor="white", lw=0.5, zorder=3)
        ax.text(i + (j - 1) * 0.26, 100 * val + 1.5, f"{100 * val:.0f}", ha="center", va="bottom",
                fontsize=K.T_SMALL - 0.3, color=K.INK)
ax.set_xticks([0, 1]); ax.set_xticklabels(list(vals)); ax.set_xlim(-0.55, 1.55)
ax.set_ylim(0, 90); ax.set_yticks([0, 25, 50, 75]); ax.set_ylabel("validated (%; log-rank P < 0.05)")
K.legend(ax, [(MOD, "published models", "s"), (COEX, "co-expressed decoys", "s"), (RAND, "random-gene decoys", "s")],
         x=0.48, y=0.95, dy=0.085, s=8, size=K.T_SMALL)

# ================================================================== c  FDR of "this model validates"
ax = C.panel("c", 146, 12, 32, 38, head="False discovery rate")
rows = [("independent, co-expressed", am("primary", "fdr")), ("independent, random genes", am("primary", "fdr", "rand")),
        ("all models", am("all models", "fdr")), ("paper used the cohort", am("models whose", "fdr"))]
for i, (lab_, r) in enumerate(rows):
    y = len(rows) - 1 - i
    K.hbar_interval(ax, y, 100 * r.lo, 100 * r.hi, K.INK2 if i else K.INK, lw=0.9)
    K.dots(ax, [100 * r.estimate], [y], COEX if "random" not in lab_ else RAND, s=16, z=5)
    ax.text(103, y, f"{100 * r.estimate:.0f}%", va="center", ha="left", fontsize=K.T_SMALL, color=K.INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=K.T_SMALL)
ax.set_xlim(0, 100); ax.set_xticks([0, 50, 100]); ax.set_xlabel("estimated FDR (%; 95% interval)")
ax.set_ylim(-0.6, len(rows) - 0.4)

# ================================================================== d  each model among its decoys
ax = C.panel("d", 24, 70, 64, 104, head="Each independent model among its decoys")
o = IND.sort_values("ext_cdir").reset_index(drop=True)
for y, r in o.iterrows():
    dc = DEC[(DEC.model == r.model) & (DEC["class"] == "coex")].ext_cdir
    q = np.quantile(dc, [0.05, 0.25, 0.5, 0.75, 0.95])
    ax.plot([q[0], q[4]], [y, y], color=K.tint(COEX, 0.45), lw=1.0, solid_capstyle="butt", zorder=2)
    ax.plot([q[1], q[3]], [y, y], color=COEX, lw=2.2, solid_capstyle="butt", zorder=3)
    K.dots(ax, [r.ext_cdir], [y], MOD, s=10, filled=bool(r.validated), z=5)
    ax.text(0.302, y, r.type, va="center", ha="left", fontsize=K.T_SMALL - 1.2, color=K.INK2)
K.ref(ax, 0.5, axis="v")
ax.set_xlim(0.3, 0.8); ax.set_ylim(-1, len(o)); ax.set_yticks([]); K.frame(ax, left=False)
ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]); ax.set_xlabel("directed C-index in the independent cohort")
K.legend(ax, [(MOD, "model, validated", "o", True), (MOD, "model, not validated", "o", False),
              (COEX, "co-expressed decoys,", "line")], x=0.70, y=0.17, dy=0.035, s=8, size=K.T_SMALL)
ax.text(0.70 + 0.035, 0.17 - 3 * 0.035, "25–75% and 5–95%", transform=ax.transAxes, va="center", ha="left",
        fontsize=K.T_SMALL, color=K.INK)

# ================================================================== e  rho
ax = C.panel("e", 110, 70, 68, 30, head="Where a model sits among its decoys")
bins = np.linspace(0, 1, 11)
for X, col, lab_, off in ((IND, MOD, f"independent (n = {len(IND)})", -0.022),
                          (U[U.uses_external.astype(bool)], K.SLATE, "paper used the cohort "
                           f"(n = {int(U.uses_external.astype(bool).sum())})", 0.022)):
    h, _ = np.histogram(X.rho, bins)
    ax.bar(bins[:-1] + 0.05 + off, 100 * h / len(X), 0.042, color=col if col == MOD else K.tint(col, 0.55),
           edgecolor="white", lw=0.4, zorder=3, label=lab_)
ax.axhline(10, color=K.RULE, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_ylabel("models (%)")
ax.set_xlabel("ρ, share of co-expressed decoys the model exceeds")
K.legend(ax, [(MOD, f"independent (n = {len(IND)})", "s"),
              (K.tint(K.SLATE, 0.55), f"paper used the cohort (n = {int(U.uses_external.astype(bool).sum())})", "s")],
         x=0.04, y=0.95, dy=0.12, s=8, size=K.T_SMALL)
ax.text(0.04, 0.62, f"mean ρ {IND.rho.mean():.2f} (independent)", transform=ax.transAxes, ha="left", va="top",
        fontsize=K.T_SMALL, color=K.INK)

# ================================================================== f  the pipeline by type
ax = C.panel("f", 110, 118, 68, 56, head="The pipeline on theme and decoy pools")
by = PU.groupby("type").agg(t=("validated", "mean"), d=("d_coex", "mean"), r=("d_rand", "mean")).sort_values("d")
for y, (tp, r) in enumerate(by.iterrows()):
    ax.plot([100 * r.d, 100 * r.t], [y, y], color=K.RULE, lw=0.9, zorder=2)
    K.dots(ax, [100 * r.d], [y], COEX, s=12, z=4)
    K.dots(ax, [100 * r.t], [y], MOD, s=12, z=5)
ax.set_yticks(range(len(by))); ax.set_yticklabels(by.index, fontsize=K.T_SMALL)
ax.set_xlim(-2, 100); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("pools yielding a validated model (%)")
fb = PA[(PA.population == "theme pools") & (PA.quantity == "fdr") & (PA.decoys == "coex")].iloc[0]
ax.text(0.98, 0.05, f"FDR {100 * fb.estimate:.0f}% ({100 * fb.lo:.0f}–{100 * fb.hi:.0f}%)", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=K.T_SMALL, color=K.INK)
K.legend(ax, [(MOD, "15 theme pools", "o"), (COEX, "20 co-expressed decoy pools each", "o")], x=0.50, y=0.66,
         dy=0.075, s=8, size=K.T_SMALL)

# ================================================================== g  the pipeline by theme pool
ax = C.panel("g", 24, 194, 154, 30, head="By theme: validated theme units against what their decoys predict")
bp = PU.groupby("pool").agg(v=("validated", "sum"), e=("d_coex", "sum")).sort_values("e")
x = np.arange(len(bp))
ax.bar(x - 0.18, bp.v, 0.34, color=MOD, edgecolor="white", lw=0.5, zorder=3)
ax.bar(x + 0.18, bp.e, 0.34, color=COEX, edgecolor="white", lw=0.5, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels([{"emt": "EMT", "m6a": "m6A"}.get(p_, p_) for p_ in bp.index], fontsize=K.T_SMALL, rotation=30,
                   ha="right")
ax.set_ylabel("of 11 cancer types"); ax.set_ylim(0, 7); ax.set_yticks([0, 2, 4, 6]); ax.set_xlim(-0.6, len(bp) - 0.4)
K.legend(ax, [(MOD, "validated theme models", "s"), (COEX, "expected from decoy pools", "s")], x=0.02, y=0.95,
         dy=0.13, s=8, size=K.T_SMALL)

C.save("Figure9_NC16")
print("Figure9_NC16")
