# Analysis plan, round 16: the false discovery rate of the integrative machine-learning pipeline

Written 2026-11-03, with the plan of round 15 and before either was started. Its SHA-256 is recorded in
PLAN_round16_ml_pipeline.sha256 and the file is deposited before the first step. Changes are made only by dated
amendments deposited before the step they change.

## Question

A widely reused pipeline (Liu et al., Nat Commun 2022;13:816, and the many papers that follow it) fits 101
combinations of 10 survival-learning algorithms (random survival forest, elastic net, lasso, ridge, stepwise Cox,
CoxBoost, partial least squares Cox, supervised principal components, generalized boosted regression, survival support
vector machine) to a gene set in a training cohort, keeps the combination with the highest mean C-index across the
validation cohorts, and reports the resulting score as validated. Round 13 asked of the LASSO-Cox pipeline how often
it validates a model from a theme's genes against how often it does from decoy pools of other genes; this round asks
the same of the integrative pipeline, whose choice among 101 fits is made on the validation cohorts themselves.

## What is known when this is written

Everything listed in PLAN_round15_literature_cohorts_controls.md. The integrative pipeline has not been run in this
project; no package implementing it is installed (of its learners only glmnet is).

## Design

1. Implementation: the published implementation, Mime1 (github.com/l-magnificence/Mime, the commit current on the day
   of installation, recorded), function ML.Dev.Prog.Sig with all 101 combinations and its default settings. If it
   cannot be installed or run, the 101 combinations are implemented from the Methods of Liu et al. with the same
   learners and default settings, stated. Learner packages are installed once from CRAN or their published sources
   and their versions recorded.

2. Units: the 11 types of round 13; per type, the 15 theme pools of round 13 and, of each pool's round-13 decoy pools
   (20 co-expressed, 5 random-gene, in their stored order), the first 6 co-expressed and the first 2 random-gene: 9
   pools per theme pool, 1,485 runs. Every pool is restricted to genes measured in TCGA and in every validation cohort
   of its type.

3. Cohorts: training, the TCGA primary tumors of the type with TCGA-CDR overall survival; validation, the eligible
   cohorts of the type from round 15 part B (the round-11 cohort always among them).

4. Pipeline, per pool: univariate Cox in TCGA, genes with P < 0.05 kept (fewer than two: no model); the 101
   combinations fitted in TCGA; the combination with the highest mean C-index over the validation cohorts kept, as
   published; seeds 20261211 + run index. The kept score is judged in each validation cohort under the round-13 rule
   (L) (median split, log-rank P < 0.05, high score faring worse).

5. Claim and estimate. Primary claim: the kept model validates (L) in the round-11 cohort of the type. Secondary: it
   validates in every validation cohort of the type. FDR_ML = sum over theme units of their decoy pools' validation
   rate divided by the number of theme units validated (the estimator of round 13 B4), 95% interval from 10,000
   bootstrap resamples of units stratified by type (seed 20261212).

6. Comparison with LASSO-Cox: the round-13 pipeline (univariate Cox, LASSO-Cox, lambda.min) rerun on the same
   restricted pools, fold ids seed 20261213 + run index, judged the same way.

## Hypotheses

* H8 (primary): FDR_ML with its interval, and the share of decoy pools that yield a validated model under the
  integrative pipeline.
* H9: theme pools validate more often than their decoy pools under the integrative pipeline (one-sided bootstrap P).
* H10: decoy pools validate more often under the integrative pipeline than under LASSO-Cox on the same pools (one-sided
  paired bootstrap P over pools stratified by type, seed 20261214).
* Secondary: the mean C-index of the kept model in the validation cohorts for theme and decoy pools; which learner
  combinations are kept, for theme and decoy pools.

## Reporting

Every result is reported whatever it shows. The pipeline is assessed, not any paper or its authors; the original paper
is cited for the pipeline's description only. Runs that fail are counted and reported with the error. Scripts:
applications/v21_ml_pipeline.R; outputs icb/v21_ml/. Numbers already reported do not change.
