"""Estimates of validation/PLAN_round15_literature_cohorts_controls.md A4 (round 13 A9-A10 on the whole search; amendments
1-2) from the outputs of applications/v20_published_models.R.

    python3 v20_audit.py   -> icb/v20_round15/audit_models*.csv
"""
import glob, numpy as np, pandas as pd
import v19_audit as V19

G = "/path/to/analysis_root"
A = f"{G}/icb/v20_round15"
SEED = 20261203


def boot_rows(X, pop, extra=None):
    rows = []
    for dc, lab in (("d_coex", "coex"), ("d_rand", "rand")):
        est = V19.fdr_parts(X, dc); bt = V19.strat_boot(X, lambda x: V19.fdr_parts(x, dc), SEED)
        rows += [dict(r, decoys=lab) for r in V19.summarize(est, bt, pop)]
        if lab == "coex":
            rows.append(dict(population=pop, quantity="H1 one-sided P (excess <= 0)",
                             estimate=float((bt.excess <= 0).mean()), decoys=lab))
    bt = V19.strat_boot(X, lambda x: dict(rho=x.rho.mean()), SEED)
    rows += [dict(r, decoys="coex") for r in V19.summarize(dict(rho=X.rho.mean()), bt, pop)]
    return rows


def main():
    R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/models_[A-Z]*.csv.gz"))])
    PM = pd.read_csv(f"{G}/signature-decoy-fdr/validation/published_models_r15.csv")
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
    T["tcga_validated"] = ((T.tcga_p < 0.05) & (T.tcga_worse == 1)).astype(int)
    T = T.join(PM[["theme", "uses_external", "pmcid", "year", "n_genes", "round", "hr_weights"]], rsuffix="_pm")
    T["family"] = T.theme.map(V19.theme_family)
    T["years"] = np.where(T.year <= 2021, "2019-2021", "2022-2025")
    T["size"] = pd.cut(T.n_genes, [0, 5, 10, 1000], labels=["2-5", "6-10", "11+"]).astype(str)
    T.reset_index().to_csv(f"{A}/audit_models_units.csv", index=False)
    U = T[T.n_coex.notna()].reset_index()
    indep = U[~U.uses_external.astype(bool)]
    rows = []
    pops = [("primary: independent of the external cohort", indep), ("all models with decoys", U),
            ("models whose paper used the external cohort", U[U.uses_external.astype(bool)]),
            ("sensitivity: independent, without hazard-ratio weights", indep[~indep.hr_weights.astype(bool)]),
            ("independent, round-13 models (new decoys)", indep[indep["round"] == 13]),
            ("independent, round-15 models", indep[indep["round"] == 15])]
    for pop, X in pops:
        if len(X):
            rows += boot_rows(X, pop)
    rows.append(dict(population="by training-set size", quantity="not run: training-set sizes were not recorded "
                                                                  "in screening", estimate=np.nan))
    pd.DataFrame(rows).to_csv(f"{A}/audit_models.csv", index=False)
    for key in ("type", "family", "years", "size"):
        for lab, X in (("independent", indep), ("all", U)):
            g = X.groupby(key).agg(models=("validated", "size"), validated=("validated", "sum"),
                                   decoy_exp=("d_coex", "sum"), rho=("rho", "mean"))
            g["validated_share"] = g.validated / g.models
            g["decoy_share"] = g.decoy_exp / g.models
            g["fdr"] = g.decoy_exp / g.validated.replace(0, np.nan)
            g.to_csv(f"{A}/audit_models_by_{key}_{lab}.csv")
    return U, indep


if __name__ == "__main__":
    U, indep = main()
    print("models with decoys", len(U), "independent", len(indep))
    print(pd.read_csv(f"{A}/audit_models.csv").to_string())
    for key in ("type", "family", "years", "size"):
        print(pd.read_csv(f"{A}/audit_models_by_{key}_independent.csv").to_string())
