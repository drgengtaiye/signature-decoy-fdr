"""Round 18 (validation/PLAN_round18_law_out_of_sample.md): H6b and H6c on halves. -> icb/v23_law_split/law_split_h6.csv, law_split_K<K>.csv"""
import glob
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
G = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"; O = f"{G}/v23_law_split"; B = 10000
L = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{O}/law_split*.csv.gz"))]).drop_duplicates(["cohort", "split", "K", "id"])
print(L.shape, L.cohort.nunique(), "cohorts")
rows = []
for K in (10, 5, 20):
    X = L[L.K == K]
    C = X.groupby("cohort").agg(pred=("p_pred", "mean"), obs=("validated", "mean"), n=("validated", "size"), k=("validated", "sum"), D_B=("D_B", "mean"), splits=("split", "nunique"), type=("type", "first")).reset_index()
    loo = []
    for i in range(len(C)):
        tr = C.drop(i)
        f = sm.GLM(np.column_stack([tr.k, tr.n - tr.k]), np.column_stack([np.ones(len(tr)), np.sqrt(tr.D_B)]), family=sm.families.Binomial()).fit()
        loo.append(float(f.predict(np.array([[1.0, np.sqrt(C.D_B[i])]]))[0]))
    C["pred_events"] = loo
    C.to_csv(f"{O}/law_split_K{K}.csv", index=False)
    rho = stats.spearmanr(C.pred, C.obs).statistic
    rng = np.random.default_rng(20261302)
    perm = np.array([stats.spearmanr(C.pred, rng.permutation(C.obs.values)).statistic for _ in range(B)])
    e_l, e_e = np.abs(C.pred - C.obs).values, np.abs(C.pred_events - C.obs).values
    rng = np.random.default_rng(20261303)
    bd = np.array([(lambda ix: e_l[ix].mean() - e_e[ix].mean())(rng.choice(len(C), len(C))) for _ in range(B)])
    within = []
    for _, g in X.groupby("cohort"):
        if 0 < g.validated.sum() < len(g):
            within.append(roc_auc_score(g.validated, g.p_pred))
    rows.append(dict(K=K, cohorts=len(C), mean_splits=C.splits.mean(), rows=len(X), spearman=rho, p_perm=float((perm >= rho).mean()), mae_law=e_l.mean(), mae_events=e_e.mean(),
                     mae_diff=e_l.mean() - e_e.mean(), mae_lo=np.quantile(bd, .025), mae_hi=np.quantile(bd, .975), p_law_not_better=float((bd >= 0).mean()),
                     spearman_events=stats.spearmanr(C.pred_events, C.obs).statistic, auc_pooled=roc_auc_score(X.validated, X.p_pred), auc_within_median=np.median(within), auc_within_q25=np.quantile(within, .25), auc_within_q75=np.quantile(within, .75), obs_mean=C.obs.mean(), pred_mean=C.pred.mean()))
R = pd.DataFrame(rows); R.to_csv(f"{O}/law_split_h6.csv", index=False); print(R.T.round(4).to_string())
# for comparison: the in-sample law of round 15, within-cohort AUC
I = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{G}/v20_round15/law_predictions*.csv.gz"))]).drop_duplicates(["cohort", "id", "K"])
S = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{G}/v20_round15/cohort_scores_*.csv.gz"))]); S["kind"] = S["kind"].fillna("null")
S = S[S.kind == "null"][["cohort", "id", "validated"]]
Y = I[I.K == 10].merge(S, on=["cohort", "id"]); w = [roc_auc_score(g.validated, g.p_pred) for _, g in Y.groupby("cohort") if 0 < g.validated.sum() < len(g)]
print("in sample, K = 10: pooled AUC", round(roc_auc_score(Y.validated, Y.p_pred), 3), "within-cohort median", round(float(np.median(w)), 3), "n", len(w))
