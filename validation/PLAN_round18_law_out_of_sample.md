# Plan, round 18: the principal-component law out of sample (2026-11-06)

Written and deposited before the split runs. Round 15 H6 is registered and stands as reported, but it fits gamma_c, the survival loadings of the K components, on the survival of the same cohort
whose validation rate it then predicts, while its events-only comparator is leave-one-cohort-out. So its Spearman 0.89 and its mean absolute error of 0.073 against 0.169 are in-sample against
out-of-sample, and they show that the law accounts for the rates, not that it predicts them. This round tests prediction.

## Design
For each of the 35 cohorts of the independent set (round 15 amendment 3), K = 10 (primary), 5, 20:
1. The components of the z-scored genes and each null model's loadings a_mc are computed on all patients, without outcome, as in applications/v20_law.R.
2. 20 random splits of the patients into halves A and B, stratified by event (seed 20261301 + cohort index). gamma is fitted by Cox regression on the K components in A only.
3. beta_mc = a_mc . gamma, and p_mc = Phi(0.80 sqrt(D_B) beta_mc - 1.96), D_B the events in B.
4. Observed in B alone: rule (L) of round 13, score split at its median in B, log-rank P < 0.05, high group worse. Null models with direction as in round 15 (higher score, shorter survival).
5. Per cohort and K, predicted FVR and observed FVR are means over null models and splits. A split whose Cox fit does not converge or with D_B < K + 1 is dropped and counted.

## Tests
H6b (primary): Spearman of predicted and observed FVR over the 35 cohorts, K = 10, one-sided permutation P (10,000, seed 20261302).
H6c (primary): mean absolute error of the law against an events-only prediction, the same leave-one-cohort-out binomial regression of observed FVR on sqrt(mean D_B) as round 15, here on the half-cohort rates; one-sided
paired bootstrap over cohorts (10,000, seed 20261303).
Secondary: K = 5 and 20; pooled AUC of p_mc for validation in B, and its within-cohort median; both also on the full cohorts in sample, for comparison.

## Reading the result, fixed now
- If H6b is positive with Spearman at least 0.6 and H6c shows the law's error below events-only: the text may say the law predicts a cohort's rate out of sample, and gives the out-of-sample numbers beside the in-sample ones.
- Otherwise: the text says the law accounts for the rates in sample, that it does not predict them from half of the patients, and it does not use the word predict for the in-sample result. The cohort table, measured directly, is then the usable output.
Whichever holds is reported. Rates in halves are lower than in full cohorts because halves have fewer events, and the law allows for that through D_B.

Output: icb/v23_law_split/. Scripts: applications/v23_law_split.R, decoys/v23_law_split_audit.py.
