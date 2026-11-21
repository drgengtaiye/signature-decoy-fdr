"""Estimates of validation/PLAN_round13_published_risk_models.md (A9-A10, B4) from the outputs of
applications/v19_published_models.R and v19_pipeline.R.

    python v19_audit.py   -> icb/v19_audit/audit_*.csv
"""
import glob, numpy as np, pandas as pd

G = "/path/to/analysis_root"
A = f"{G}/icb/v19_audit"
TYPES = ["STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA"]
B = 10000


def strat_boot(units, fn, seed):
    """units: DataFrame with a 'type' column; fn(DataFrame) -> dict of numbers; resampled within type."""
    rng = np.random.default_rng(seed)
    idx = {t: np.flatnonzero(units.type.values == t) for t in units.type.unique()}
    out = []
    for _ in range(B):
        pick = np.concatenate([rng.choice(ix, len(ix), replace=True) for ix in idx.values()])
        out.append(fn(units.iloc[pick]))
    return pd.DataFrame(out)


def summarize(est, boot, label):
    rows = []
    for k, v in est.items():
        b = boot[k].replace([np.inf, -np.inf], np.nan).dropna()
        rows.append(dict(population=label, quantity=k, estimate=v, lo=b.quantile(0.025) if len(b) else np.nan,
                         hi=b.quantile(0.975) if len(b) else np.nan, n_boot=len(b)))
    return rows


def fdr_parts(u, dcol):
    V = u.validated.sum(); E = u[dcol].sum(); M = len(u)
    return dict(M=M, validated=V, validated_share=V / M, decoy_share=E / M, excess=(V - E) / M,
                fdr=(E / V) if V > 0 else np.nan)


# ------------------------------------------------------------------ part B: the pipeline on theme pools
def part_b():
    T = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/pipeline_targets_*.csv.gz"))])
    D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/pipeline_decoys_*.csv.gz"))])
    for X in (T, D):
        X["made"] = (X.n_model >= 2).astype(int)
    agg = D.groupby(["type", "pool", "class"]).agg(d=("validated", "mean"), d_made=("made", "mean")).unstack("class")
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    U = T.set_index(["type", "pool"]).join(agg).reset_index()
    U.to_csv(f"{A}/audit_pipeline_units.csv", index=False)
    rows = []
    for dc, lab in (("d_coex", "coex"), ("d_rand", "rand")):
        est = fdr_parts(U, dc); bt = strat_boot(U, lambda x: fdr_parts(x, dc), 20261119)
        rows += [dict(r, decoys=lab) for r in summarize(est, bt, "theme pools")]
        if lab == "coex":
            p_h2 = float((bt.excess <= 0).mean())
            rows.append(dict(population="theme pools", quantity="H2 one-sided P (excess <= 0)", estimate=p_h2,
                             decoys=lab))
    rows.append(dict(population="decoy pools", quantity="share validated (all co-expressed decoy pools)",
                     estimate=D[D["class"] == "coex"].validated.mean(), decoys="coex"))
    rows.append(dict(population="decoy pools", quantity="share validated (all random-gene decoy pools)",
                     estimate=D[D["class"] == "rand"].validated.mean(), decoys="rand"))
    rows.append(dict(population="theme pools", quantity="share making a model", estimate=T.made.mean()))
    rows.append(dict(population="decoy pools", quantity="share making a model (co-expressed)",
                     estimate=D[D["class"] == "coex"].made.mean(), decoys="coex"))
    pd.DataFrame(rows).to_csv(f"{A}/audit_pipeline.csv", index=False)
    by = U.groupby("type").agg(units=("validated", "size"), validated=("validated", "sum"), decoy_exp=("d_coex", "sum"))
    by["fdr"] = by.decoy_exp / by.validated.replace(0, np.nan)
    by.to_csv(f"{A}/audit_pipeline_by_type.csv")
    byp = U.groupby("pool").agg(units=("validated", "size"), validated=("validated", "sum"), decoy_exp=("d_coex", "sum"))
    byp["fdr"] = byp.decoy_exp / byp.validated.replace(0, np.nan)
    byp.to_csv(f"{A}/audit_pipeline_by_pool.csv")
    return U


# ------------------------------------------------------------------ part A: published models
# theme family (A10): the first of the plan's theme terms found in the paper's recorded theme, else "other"
FAMILY_TERMS = [("ferroptosis", "ferroptosis"), ("cuproptosis", "cuproptosis"), ("pyroptosis", "pyroptosis"),
                ("necroptosis", "necroptosis"), ("anoikis", "anoikis"), ("disulfidptosis", "disulfidptosis"),
                ("mitophagy", "mitophagy"), ("autophagy", "autophagy"), ("m6a", "m6A"),
                ("n6-methyladenosine", "m6A"), ("hypoxia", "hypoxia"), ("glycolysis", "glycolysis"),
                ("senescence", "senescence"), ("immune", "immune-related"), ("epithelial-mesenchymal", "EMT")]


def theme_family(t):
    t = str(t).lower()
    hits = [(t.find(k), f) for k, f in FAMILY_TERMS if k in t]
    return min(hits)[1] if hits else "other"


def part_a():
    R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/models_*.csv.gz")) if "decoy_info" not in f])
    PM = pd.read_csv(f"{G}/signature-decoy-fdr/validation/published_models.csv")
    PM = PM[PM.scored.astype(bool)].set_index("model_id")
    T = R[R.kind == "target"].set_index("model")
    Dc = R[(R.kind == "decoy") & (R["class"] == "coex")]
    Dr = R[(R.kind == "decoy") & (R["class"] == "rand")]
    T["d_coex"] = Dc.groupby("model").validated.mean()
    T["d_rand"] = Dr.groupby("model").validated.mean()
    T["n_coex"] = Dc.groupby("model").size()
    m = Dc.merge(T.ext_cdir.rename("t_c"), left_on="model", right_index=True)
    m["below"] = (m.ext_cdir < m.t_c) + 0.5 * (m.ext_cdir == m.t_c)
    T["rho"] = m.groupby("model").below.mean()
    T["p95"] = Dc.groupby("model").ext_cdir.quantile(0.95)
    T["call_C"] = (T.ext_cdir > T.p95).astype(int)
    T["tcga_validated"] = ((T.tcga_p < 0.05) & (T.tcga_worse == 1)).astype(int)
    T["d_tcga"] = Dc.assign(v=((Dc.tcga_p < 0.05) & (Dc.tcga_worse == 1)).astype(int)).groupby("model").v.mean()
    T = T.join(PM[["theme", "uses_external", "pmcid", "year", "n_genes"]])
    T["family"] = T.theme.map(theme_family)
    T.reset_index().to_csv(f"{A}/audit_models_units.csv", index=False)
    U = T[T.n_coex.notna()].reset_index()                       # models with decoys (three or more genes)
    indep = U[~U.uses_external.astype(bool)]
    primary = indep if len(indep) >= 40 else U
    plabel = "independent of the external cohort" if len(indep) >= 40 else "all models (fewer than 40 independent)"
    rows = []
    for pop, X in ((f"primary: {plabel}", primary), ("all models with decoys", U),
                   ("models whose paper used the external cohort", U[U.uses_external.astype(bool)])):
        if not len(X):
            continue
        for dc, lab in (("d_coex", "coex"), ("d_rand", "rand")):
            est = fdr_parts(X, dc); bt = strat_boot(X, lambda x: fdr_parts(x, dc), 20261116)
            rows += [dict(r, decoys=lab) for r in summarize(est, bt, pop)]
            if lab == "coex":
                rows.append(dict(population=pop, quantity="H1 one-sided P (excess <= 0)",
                                 estimate=float((bt.excess <= 0).mean()), decoys=lab))
        bt = strat_boot(X, lambda x: dict(rho=x.rho.mean()), 20261116)
        rows += [dict(r, decoys="coex") for r in summarize(dict(rho=X.rho.mean()), bt, pop)]
        # (C): directed C-index above the co-expressed decoys' 95th percentile; decoys clear it 5% by construction
        vC = X.call_C.sum()
        rows.append(dict(population=pop, quantity="calls by C above decoy p95", estimate=vC, decoys="coex"))
        rows.append(dict(population=pop, quantity="fdr by C (0.05 M / calls)",
                         estimate=0.05 * len(X) / vC if vC else np.nan, decoys="coex"))
        rows.append(dict(population=pop, quantity="TCGA (in-sample) validated share", estimate=X.tcga_validated.mean()))
        rows.append(dict(population=pop, quantity="TCGA (in-sample) decoy validated share", estimate=X.d_tcga.mean()))
    pd.DataFrame(rows).to_csv(f"{A}/audit_models.csv", index=False)
    for key in ("family", "type"):
        g = U.groupby(key).agg(models=("validated", "size"), validated=("validated", "sum"),
                               decoy_exp=("d_coex", "sum"), rho=("rho", "mean"))
        g["fdr"] = g.decoy_exp / g.validated.replace(0, np.nan)
        g.to_csv(f"{A}/audit_models_by_{key}.csv")
    return U


if __name__ == "__main__":
    if glob.glob(f"{A}/pipeline_decoys_*.csv.gz"):
        u = part_b(); print("part B", len(u), "units")
        print(pd.read_csv(f"{A}/audit_pipeline.csv").to_string())
    if glob.glob(f"{A}/models_[A-Z]*.csv.gz"):
        u = part_a(); print("part A", len(u), "models")
        print(pd.read_csv(f"{A}/audit_models.csv").to_string())
