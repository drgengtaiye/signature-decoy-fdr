"""Round 15 parts B4-B6 and C2 (validation/PLAN_round15_literature_cohorts_controls.md; amendment 3) from
applications/v20_cohort_scores.R and v20_law.R.

  B4  false-validation rate of every built cohort (null models validated under (L)), with the published models' and
      their decoys' validation rates there
  H4  within a type, cohorts' rates differ: Cochran's Q over null models scorable in every independent cohort, summed
  H5  a decoy check in the cohort where a model validated predicts replication in another (bootstrap, seed 20261206)
  H6  the law's predicted rate tracks the observed one (permutation, seed 20261207) and beats an events-only prediction
      (paired bootstrap over cohorts, seed 20261208)

A paper names a cohort when its text names the cohort's accession or one of the series that share its patients
(amendment 3, decision 6; amendment 2 of round 13, item 5).

    python3 v20_cohorts_audit.py   -> icb/v20_round15/cohorts_*.csv
"""
import glob, re
import numpy as np, pandas as pd
from scipy import stats

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
A = f"{G}/icb/v20_round15"
B = 10000
CO = pd.read_csv(f"{G}/signature-decoy-fdr/validation/round15_cohorts/built_cohorts_r15.csv")
PM = pd.read_csv(f"{G}/signature-decoy-fdr/validation/published_models_r15.csv").set_index("model_id")
GROUPS = [{"GSE84437", "GSE84433", "GSE84426", "GSE147163"}, {"GSE39582", "GSE40967"},
          {"GSE17536", "GSE17537", "GSE17538", "GSE29621", "GSE14333"}, {"GSE58812", "GSE103091", "GSE83937"},
          {"GSE20711", "GSE20713"}, {"GSE96058", "GSE81540"}, {"GSE157011", "GSE157009", "GSE157010"},
          {"GSE41271", "GSE42127"}, {"GSE13213", "GSE11969"}, {"GSE26939", "GSE36471"}, {"GSE21050", "GSE71118"},
          {"OV-AU", "ICGC"}, {"HTMCP-CC", "CGCI-HTMCP-CC", "HTMCP", "CGCI"}, {"METABRIC"}]


def acc_of(cohort):
    return cohort.split("_", 1)[1]


def names_of(cohort):
    a = acc_of(cohort)
    for g in GROUPS:
        if a in g:
            return g
    return {a}


def paper_names(model, cohort):
    txt = str(PM.at[model, "validation_cohorts"]) if model in PM.index else ""
    return any(re.search(r"(?<![A-Za-z0-9])" + re.escape(n) + r"(?![A-Za-z0-9])", txt) for n in names_of(cohort))


S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/cohort_scores_*.csv.gz"))], ignore_index=True)
S["kind"] = S["kind"].fillna("null")          # pandas reads the label "null" as missing
S = S[S.scorable.astype(bool)].copy()
S["validated"] = S.validated.astype(float)
built = CO[CO.built.astype(bool)]
indep = set(built.cohort[built.independent_set.astype(bool)])

# ------------------------------------------------------------------ B4
rows = []
for _, c in built.iterrows():
    x = S[S.cohort == c.cohort]
    nul = x[x.kind == "null"]; tg = x[x.kind == "target"]; dc = x[(x.kind == "decoy") & (x["class"] == "coex")]
    tg_unnamed = tg[[not paper_names(m, c.cohort) for m in tg.id]]
    dc_unnamed = dc[dc.model.isin(tg_unnamed.id)]
    rows.append(dict(type=c.type, cohort=c.cohort, endpoint=c.endpoint, patients=c.patients, events=c.events,
                     independent_set=c.cohort in indep, null_models=len(nul), fvr=nul.validated.mean(),
                     fvr_coex=nul[nul["class"] == "coex"].validated.mean(),
                     fvr_rand=nul[nul["class"] == "rand"].validated.mean(),
                     models_not_naming=len(tg_unnamed), model_rate_not_naming=tg_unnamed.validated.mean(),
                     decoy_rate_not_naming=dc_unnamed.validated.mean(),
                     models_naming=len(tg) - len(tg_unnamed),
                     model_rate_naming=tg[~tg.id.isin(tg_unnamed.id)].validated.mean()))
for _, c in CO[~CO.built.astype(bool) & CO.note.str.startswith("same patients")].iterrows():
    twin = re.search(r"as (\S+):", c.note).group(1)
    r = dict(next(r for r in rows if r["cohort"] == twin)); r.update(cohort=c.cohort, independent_set=False,
                                                                       note=f"rate of {twin}")
    rows.append(r)
F = pd.DataFrame(rows).sort_values(["type", "fvr"], ascending=[True, False])
F.to_csv(f"{A}/cohorts_fvr.csv", index=False)

# ------------------------------------------------------------------ H4: Cochran's Q within types
qrows, Qsum, dfsum = [], 0.0, 0
for tp, cs in built[built.cohort.isin(indep)].groupby("type"):
    cs = list(cs.cohort)
    if len(cs) < 2:
        continue
    x = S[(S.kind == "null") & S.cohort.isin(cs)].pivot_table(index="id", columns="cohort", values="validated")
    x = x.dropna()
    k, n = x.shape[1], x.shape[0]
    Cj, Ri, N = x.sum(0).values, x.sum(1).values, x.values.sum()
    den = k * N - (Ri ** 2).sum()
    Q = (k - 1) * (k * (Cj ** 2).sum() - N ** 2) / den if den > 0 else 0.0
    qrows.append(dict(type=tp, cohorts=k, null_models=n, Q=Q, df=k - 1, p=stats.chi2.sf(Q, k - 1),
                      rates=";".join(f"{c}={v:.3f}" for c, v in zip(x.columns, x.mean(0)))))
    Qsum += Q; dfsum += k - 1
qrows.append(dict(type="all (sum)", cohorts=np.nan, null_models=np.nan, Q=Qsum, df=dfsum, p=stats.chi2.sf(Qsum, dfsum)))
FI = F[F.independent_set]
rho_e = stats.spearmanr(FI.fvr, FI.events, alternative="greater")
qrows.append(dict(type="Spearman FVR vs events (independent set)", cohorts=len(FI), Q=rho_e.statistic, p=rho_e.pvalue))
pd.DataFrame(qrows).to_csv(f"{A}/cohorts_h4.csv", index=False)

# ------------------------------------------------------------------ H5: decoy check in A predicts replication in B
T = S[(S.kind == "target") & S.cohort.isin(indep)][["type", "cohort", "id", "cdir", "validated"]]
D = S[(S.kind == "decoy") & (S["class"] == "coex") & S.cohort.isin(indep)][["cohort", "model", "cdir"]]
m = D.merge(T.rename(columns={"id": "model", "cdir": "t_c"})[["cohort", "model", "t_c"]], on=["cohort", "model"])
m["below"] = (m.cdir < m.t_c) + 0.5 * (m.cdir == m.t_c)
rho = m.groupby(["model", "cohort"]).below.mean().rename("rho_A")
T = T.join(rho, on=["id", "cohort"])
pairs = []
for (tp, mod), x in T.groupby(["type", "id"]):
    va = x[(x.validated == 1) & x.rho_A.notna()]
    for _, a in va.iterrows():
        for _, b in x[x.cohort != a.cohort].iterrows():
            if not paper_names(mod, b.cohort):
                pairs.append(dict(type=tp, model=mod, A=a.cohort, B=b.cohort, rho_A=a.rho_A, pass_A=a.rho_A >= 0.95,
                                  replicated=b.validated))
P = pd.DataFrame(pairs)
P["w"] = 1.0 / P.groupby("model").model.transform("size")
P.to_csv(f"{A}/cohorts_h5_pairs.csv", index=False)


def h5_stat(p):
    g1, g0 = p[p.pass_A], p[~p.pass_A]
    r1 = np.average(g1.replicated, weights=g1.w) if len(g1) else np.nan
    r0 = np.average(g0.replicated, weights=g0.w) if len(g0) else np.nan
    return r1, r0, r1 - r0


r1, r0, d = h5_stat(P)
rng = np.random.default_rng(20261206)
mods = P[["type", "model"]].drop_duplicates()
by_t = {t: g.model.values for t, g in mods.groupby("type")}
grp = {k: g for k, g in P.groupby("model")}
boot = []
for _ in range(B):
    pick = np.concatenate([rng.choice(v, len(v), replace=True) for v in by_t.values()])
    bp = pd.concat([grp[k] for k in pick])
    boot.append(h5_stat(bp)[2])
boot = np.array(boot)
h5 = [dict(quantity="pairs", estimate=len(P)), dict(quantity="models", estimate=P.model.nunique()),
      dict(quantity="pairs with rho_A >= 0.95", estimate=int(P.pass_A.sum())),
      dict(quantity="replication rate, rho_A >= 0.95", estimate=r1),
      dict(quantity="replication rate, rho_A < 0.95", estimate=r0),
      dict(quantity="difference", estimate=d, lo=np.nanquantile(boot, 0.025), hi=np.nanquantile(boot, 0.975)),
      dict(quantity="H5 one-sided P (difference <= 0)", estimate=float(np.nanmean(boot <= 0)))]
# secondary: replication in B of null models and of co-expressed decoys validated in A
for kind, lab in (("null", "null models"), ("decoy", "co-expressed decoys")):
    X = S[(S.kind == kind) & S.cohort.isin(indep)]
    if kind == "decoy":
        X = X[X["class"] == "coex"]
    X = X.merge(built[["cohort", "type"]], on=["cohort", "type"])
    piv = X.pivot_table(index="id", columns="cohort", values="validated")
    num = den = 0.0
    for cA in piv.columns:
        va = piv[piv[cA] == 1]
        for cB in piv.columns:
            if cB != cA and cA.split("_")[0] == cB.split("_")[0]:
                v = va[cB].dropna(); num += v.sum(); den += len(v)
    h5.append(dict(quantity=f"replication rate in B of {lab} validated in A", estimate=num / den if den else np.nan,
                   n=den))
pd.DataFrame(h5).to_csv(f"{A}/cohorts_h5.csv", index=False)

# ------------------------------------------------------------------ H6: the law
L = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A}/law_predictions*.csv.gz"))]).drop_duplicates(["cohort", "id", "K"])
obs = S[(S.kind == "null")][["cohort", "id", "validated"]]
h6 = []
for K in (10, 5, 20):
    LK = L[L.K == K].merge(obs, on=["cohort", "id"])
    C = LK.groupby("cohort").agg(pred=("p_pred", "mean"), obs=("validated", "mean"), n=("validated", "size"),
                                  k=("validated", "sum"), events=("events", "first")).reset_index()
    # events-only: leave-one-cohort-out binomial regression of the observed rate on sqrt(events)
    loo = []
    for i in range(len(C)):
        tr = C.drop(i)
        Xtr = np.column_stack([np.ones(len(tr)), np.sqrt(tr.events)])
        import statsmodels.api as sm
        f = sm.GLM(np.column_stack([tr.k, tr.n - tr.k]), Xtr, family=sm.families.Binomial()).fit()
        loo.append(float(f.predict(np.array([[1.0, np.sqrt(C.events[i])]]))[0]))
    C["pred_events"] = loo
    C.to_csv(f"{A}/cohorts_law_K{K}.csv", index=False)
    rho_l = stats.spearmanr(C.pred, C.obs).statistic
    rng = np.random.default_rng(20261207)
    perm = np.array([stats.spearmanr(C.pred, rng.permutation(C.obs.values)).statistic for _ in range(B)])
    e_law, e_ev = np.abs(C.pred - C.obs).values, np.abs(C.pred_events - C.obs).values
    rng = np.random.default_rng(20261208)
    bd = []
    for _ in range(B):
        ix = rng.choice(len(C), len(C), replace=True); bd.append(e_law[ix].mean() - e_ev[ix].mean())
    bd = np.array(bd)
    LK["y"] = LK.validated
    try:
        from sklearn.metrics import roc_auc_score
        auc = roc_auc_score(LK.y, LK.p_pred)
    except Exception:
        auc = stats.mannwhitneyu(LK.p_pred[LK.y == 1], LK.p_pred[LK.y == 0]).statistic / ((LK.y == 1).sum() * (LK.y == 0).sum())
    h6.append(dict(K=K, cohorts=len(C), spearman=rho_l, p_perm_one_sided=float((perm >= rho_l).mean()),
                   mae_law=e_law.mean(), mae_events=e_ev.mean(), mae_diff=e_law.mean() - e_ev.mean(),
                   mae_diff_lo=np.quantile(bd, 0.025), mae_diff_hi=np.quantile(bd, 0.975),
                   p_law_not_better=float((bd >= 0).mean()), auc_null_models=auc, pairs=len(LK)))
pd.DataFrame(h6).to_csv(f"{A}/cohorts_h6.csv", index=False)

if __name__ == "__main__":
    pd.set_option("display.width", 250)
    print(F[["type", "cohort", "patients", "events", "independent_set", "null_models", "fvr", "models_not_naming",
             "model_rate_not_naming", "decoy_rate_not_naming"]].to_string())
    print(pd.DataFrame(qrows).to_string()); print(pd.DataFrame(h5).to_string()); print(pd.DataFrame(h6).to_string())
