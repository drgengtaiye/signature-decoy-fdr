"""Round 15 part D (plan D3, amendment 4): H7, the assays' decoy FDR of "validates" (round-13 estimator, pooled over assays and cohorts) against the independent
published models of A4; one-sided P from 10,000 bootstrap resamples, assays and published models resampled separately (seed 20261210). Primary cohorts: TCGA
and the independent set. Sensitivities (amendment 4 items 2-3), same estimator and bootstrap: every built cohort (independent set plus the overlapping ones),
and ER-positive patients of the annotated cohorts. Secondary: each assay's rho per cohort.
-> icb/v20_round15/assay_units.csv (primary), assay_units_all.csv (every set), assay_h7.csv (primary), assay_h7_sensitivity.csv,
   assay_by_assay.csv (primary), assay_by_assay_sensitivity.csv"""
import numpy as np, pandas as pd
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v20_round15"
S = pd.read_csv(f"{O}/assay_scores.csv.gz")
S = S[S.scorable.astype(bool)]
CO = pd.read_csv("/groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/validation/round15_cohorts/built_cohorts_r15.csv")
IND = set(CO[CO.built.astype(bool) & CO.independent_set.astype(bool)].cohort) | {"TCGA"}     # amendment 4: primary on the independent set and TCGA
S["set"] = ["erpos" if c.endswith("_ERpos") else ("primary" if c in IND else "overlap") for c in S.cohort]
M = pd.read_csv(f"{O}/audit_models_units.csv"); M = M[M.n_coex.notna() & ~M.uses_external.astype(bool)]
fdr = lambda x, col="d_coex": x[col].sum() / x.validated.sum() if x.validated.sum() else np.nan  # noqa: E731


def units(S):
    T = S[S.kind == "target"][["assay", "type", "cohort", "set", "measured", "k", "cdir", "p", "worse", "validated"]]
    D = S[S.kind == "decoy"]
    g = D.groupby(["assay", "cohort", "class"]).validated.mean().unstack("class").rename(columns={"coex": "d_coex", "rand": "d_rand"}).reset_index()
    c = D[D["class"] == "coex"].merge(T[["assay", "cohort", "cdir"]], on=["assay", "cohort"], suffixes=("", "_t"))
    rho = c.assign(x=c.cdir < c.cdir_t).groupby(["assay", "cohort"]).x.mean().rename("rho").reset_index()
    n = D[D["class"] == "coex"].groupby(["assay", "cohort"]).size().rename("n_coex").reset_index()
    return T.merge(g, on=["assay", "cohort"]).merge(rho, on=["assay", "cohort"]).merge(n, on=["assay", "cohort"])


def h7(U, label):
    rng = np.random.default_rng(20261210)                         # each set restarts the registered seed
    AS = sorted(U.assay.unique()); B = 10000
    by_assay = {a: (U[U.assay == a].validated.sum(), U[U.assay == a].d_coex.sum()) for a in AS}
    types = {t: x[["validated", "d_coex"]].to_numpy(float) for t, x in M.groupby("type")}
    fa, fm = [], []
    for _ in range(B):
        pick = rng.choice(len(AS), len(AS))
        v = sum(by_assay[AS[i]][0] for i in pick); d = sum(by_assay[AS[i]][1] for i in pick)
        fa.append(d / v if v else np.inf)
        vm = dm = 0.0
        for t, arr in types.items():                              # models resampled within type, as in A4
            ii = rng.integers(0, len(arr), len(arr)); vm += arr[ii, 0].sum(); dm += arr[ii, 1].sum()
        fm.append(dm / vm)
    fa, fm = np.array(fa), np.array(fm); finite = np.isfinite(fa)
    return dict(set=label, assays=len(AS), cohorts=U.cohort.nunique(), units=len(U), validated=int(U.validated.sum()), decoy_expected=U.d_coex.sum(),
                fdr_assays=fdr(U), fdr_assays_rand=fdr(U, "d_rand"), fdr_lo=np.quantile(fa[finite], .025), fdr_hi=np.quantile(fa[finite], .975),
                resamples_no_validation=int((~finite).sum()), mean_rho=U.rho.mean(), fdr_models=fdr(M), fdr_models_lo=np.quantile(fm, .025),
                fdr_models_hi=np.quantile(fm, .975), difference=fdr(U) - fdr(M), diff_lo=np.quantile((fa - fm)[finite], .025),
                diff_hi=np.quantile((fa - fm)[finite], .975), p_one_sided=float(np.mean(fa >= fm)))


def by_assay(U, label):
    A = U.groupby("assay").agg(type=("type", "first"), cohorts=("cohort", "size"), validated=("validated", "sum"), decoy_expected=("d_coex", "sum"),
                               decoy_rand=("d_rand", "sum"), mean_rho=("rho", "mean"), median_cdir=("cdir", "median")).reset_index()
    A["fdr"] = A.decoy_expected / A.validated.replace(0, np.nan)
    return A.assign(set=label)


ALL = units(S); ALL.to_csv(f"{O}/assay_units_all.csv", index=False)
SETS = {"primary": ALL[ALL.set == "primary"], "all_built": ALL[ALL.set.isin(["primary", "overlap"])], "erpos": ALL[ALL.set == "erpos"]}
SETS["primary"].to_csv(f"{O}/assay_units.csv", index=False)
R = pd.DataFrame([h7(u, k) for k, u in SETS.items()])
R[R.set == "primary"].drop(columns="set").to_csv(f"{O}/assay_h7.csv", index=False)
R[R.set != "primary"].to_csv(f"{O}/assay_h7_sensitivity.csv", index=False)
print(R.set_index("set").T.round(4).to_string())
BA = pd.concat([by_assay(u, k) for k, u in SETS.items()])
BA[BA.set == "primary"].drop(columns="set").to_csv(f"{O}/assay_by_assay.csv", index=False)
BA[BA.set != "primary"].to_csv(f"{O}/assay_by_assay_sensitivity.csv", index=False)
print(BA.round(3).to_string(index=False))
