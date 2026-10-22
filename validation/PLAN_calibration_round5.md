# Analysis plan: weighted-mean null, known-truth calibration of the target-decoy FDR, dependence-aware intervals, origin-cohort exclusion and the randomized interaction null

Written 2026-10-21, after an external review of manuscript v4 and before any analysis below was run.
Its SHA-256 is recorded in PLAN_calibration_round5.sha256 and the file is deposited before the first
analysis starts.

## What is known when this is written

Every result in manuscript v4. In addition, while reading the code to write this plan:

* The curated weights of the 12 weighted-mean signatures were tabulated. Eight carry unit weights
  (CTL_Jiang, GeparSixto_Karn, IFN_Ayers, IRG_Ayers, proliferation_Pabla, PTEN-MITF_Cabrita,
  TGFB_Mariathasan, Tcell-inflamed_Ayers), for which the unweighted mean of the null and the weighted
  mean of the signature are the same score. Four carry weights of +1 and -1 (EMT_Thompson,
  IMPRES_Auslander, IRG_Yang, MPS_Guijarro); for them the response, prognosis and survival-under-blockade
  benchmarks scored the null without the signs.
* In the response benchmark those four have Wilcoxon P of 0.27, 0.42, 0.21 and 0.34 against the
  published null and 0.15, 0.38, 0.23 and 0.15 against the exact-size N0 bank of the target-decoy
  calibration, which carries the signs. Neither benchmark has been recomputed.
* The one signature called by the N0 bank and not by the benchmark is CYT_Davoli (GSVA, 7 genes, null
  rung 7; P = 0.052 against 0.043). The benchmark draws 200 sets per rung; the manuscript's "300 draws"
  is to be checked.
* The 623 signature x cohort tests are 52 x 12 less CRMA_Shukla in ICB_Miao1.
* No provenance information about where each signature was developed has been compiled; the rule
  for using it is fixed below before the matrix is read.

## 1. Signed null for weighted-mean signatures

For each affected signature (the four above) and each benchmark — response in the 13 checkpoint
cohorts (native_response_matched_v2.csv), prognosis in 20 TCGA types x 3 endpoints
(native_prognosis_matched.csv), survival under checkpoint blockade (native_icbsurv_matched.csv) — 200
random sets are drawn at the signature's own rung K from the same gene pool. Of the K genes, round(K m/k)
chosen at random receive weight -1 and the rest +1 (m the number of negatively weighted genes among the k
the signature has in that cohort), and the set is scored as crossprod(Z, w)/K, the signature's own
formula. The statistic is the benchmark's own. rand_med, rand_p95, rand_mean and rand_sd of those rows are
replaced; every other row, and every observed score, is left as it is. Seeds: 20261021 + the sum of the
character codes of the cohort and the signature names.

Everything downstream is recomputed from the patched tables: the response count and its Wilcoxon P
values, the label-permutation inference (the signed sets added to the stored score matrices),
v3_dissociation.csv, the standardized effects, the prognostic and survival-under-blockade results, and
Figs. 1, 3 and 7. The corrected values replace the published ones whatever their direction, and Note S2
lists each change.

## 2. Why Fig. 3 and Fig. 4 call 20 and 21

Four constructions of the folded-AUC null are compared on the 12 cohorts, all with the same observed
AUCs: (i) ladder rung, unsigned, 200 draws (Fig. 3 as published); (ii) ladder rung, signed (analysis 1),
200 draws; (iii) exact size, signed, 1,000 draws (the N0 bank of Fig. 4); (iv) exact size, signed, 200
draws, as 500 random subsamples of the 1,000-draw bank (seed 20261021). For each the number of signatures
called (two-sided Wilcoxon P < 0.05, as in Fig. 3; and P < 0.05 with a positive median, as in Fig. 4) is
reported, and under (iv) the distribution of the count and the call frequency of every signature whose
call differs between constructions. The text then states which null each figure uses and attributes the
difference to whichever factor moves the count.

## 3. Calibration under the global null (all targets false)

A per-cohort score store is built: the 52 targets, the 2,912 decoys and the N0 1,000-draw bank, the
bank redrawn by the code, seed and draw order of null_matched_v3_coexpr.R and the decoys scored by the
code and seed of calib_score_v3.R. Reproduction is required before use: every stored N0 null AUC and every
decoy AUC in calib_v3 within 1e-10.

Response labels are permuted within each cohort, 1,000 times (seed 20261021), and the whole target-decoy
pipeline is rerun on each permutation: AUCs of targets, decoys and bank sets; the bank median transformed
per statistic (folded, directed); the Wilcoxon signed-rank test over the 12 cohorts; a call when P < 0.05
with a positive median. Every call is false. For each permutation b: R_b, the number of targets called,
and E_b, the sum over targets of the call rate among their decoys (pool-A co-expressed decoys, primary;
pool-B co-expressed and pool-A random-gene decoys, secondary; for the directed statistic the decoys
declared in their target's direction).

Primary metric: rho = mean(E_b) / mean(R_b), with a 95% interval from resampling permutations, for each
statistic. Secondary: the per-test call rate of targets, co-expressed decoys and random-gene decoys under
the label null, against the nominal 2.5%; the per-target agreement between its call frequency and its
decoys' mean call rate.

Decision: the estimator is described as calibrated under the global null if rho lies within 0.8–1.25 for
both statistics; otherwise the direction (rho > 1 conservative, rho < 1 anti-conservative) and size are
stated wherever the FDR is reported.

## 4. Known truth with signal (semi-synthetic)

The score store above plus each cohort's gene-wise z-scores. In each cohort the number of responders is
kept and responders are drawn without replacement with probability proportional to exp(beta u), u the
standardized mean z-score of a causal gene set G*.

Scenarios, 200 replicates each (seed 20261021 + scenario index):
* P1, inflamed program: G* the genes of 3 targets drawn at random from the inflamed/IFN family.
* P2, non-inflamed program: G* the genes of 3 targets drawn at random from the confounder and other
  families.
* P3, program outside every signature: G* the genes of one pool-B co-expressed decoy drawn at random;
  all 52 targets are then false, and estimation uses only the pool-A decoys sharing no gene with G*.
* beta = 0.5 and 1.0 for each of P1–P3, on all 12 cohorts; and P1, P2 at beta = 1.0 on 6 cohorts drawn
  at random per replicate.

Truth: a target is true when at least 10% of its measured genes lie in G* (primary); at least one gene,
and at least 25%, are sensitivity definitions. In P3 every target is false.

Estimators: pool-A co-expressed decoys (primary); pool-A random-gene decoys (co-expression not matched);
pool-A co-expressed decoys whose internal r lies within 0.05 of the target's (closest matched).

Metrics per scenario and statistic (folded, directed): mean realized FDP (V/R, 0 when R = 0) against the
mean estimated FDR; mean V against mean E; coverage of the realized FDP by the target bootstrap 95%
interval and by the cluster bootstrap interval of analysis 5.

Decision: the estimator is described as approximately unbiased in a scenario when the mean estimated
FDR lies within 5 percentage points of the mean realized FDP. If it underestimates by more than that in
any scenario, the title, abstract and conclusions describe the FDR as an empirical estimate referenced to
matched decoys, with the scenario stated, and no property of control is claimed.

## 5. Dependence-aware intervals for the observed FDR

For N0 folded and N0 directed (primary), and the closest-matched and natural-pool N2 directed estimates
(secondary): (a) the published bootstrap over the 52 targets; (b) a bootstrap over the 13 correlation
clusters of Fig. 2d, each carried with its targets and their decoys; (c) Monte Carlo error of the decoys
alone, targets fixed and decoys resampled within each target; (d) a bootstrap over the 12 cohorts, every
Wilcoxon call of targets and decoys recomputed; (e) clusters and cohorts resampled together. 6,000
replicates for (a)–(c), 2,000 for (d) and (e). The stability of (b) is reported with clusters cut at
1 - |r| = 0.4, 0.5 and 0.6.

Decision: the main text reports (e) as the interval for each FDR and the others in a table.

## 6. Leaving out each signature's origin cohorts

A provenance matrix is compiled from the original publications, classifying each of the 12 cohorts for
each signature as D (its outcome data were used to select genes, fit weights or train the signature),
C (the gene set was predefined but the published association with checkpoint outcome was established in
that cohort), V (used only to test the finished signature) or unused. Primary: for each signature the D
and C cohorts are removed from the target's test and from the tests of all its decoys, and the calls and
the FDR (N0 folded and directed) are recomputed. Secondary: D, C and V removed. The count and the FDR are
reported against the published ones whatever their direction; signatures left with fewer than 6 cohorts
are reported as not testable.

## 7. Randomized interaction

(a) No-interaction simulation for the pooled tests. In each survival trial PFS is drawn from the fitted
model PFS ~ score + arm (v7_lib.reduced_simulator: Breslow baseline, observed follow-up as censoring), and
in each pathological-response comparison pCR from the fitted logistic pCR ~ score + arm; the pooled models
are refitted on each of 2,000 replicates and the calibrated one-sided P is the share of replicates whose
interaction z is at least as extreme in the benefit direction. Applied to: the first-four plan (IMmotion150
+ PRINCE stratified Cox, and the I-SPY2 logistic), the three-trial primary, and the seven-comparison
combination (Stouffer). This P replaces the arm-label permutation as the reported check on the Wald P;
the permutation is kept and described as testing the sharp null of no treatment effect.

(b) Heterogeneity: the per-trial interaction (log HR, log OR) with Cochran's Q, I^2 and the
DerSimonian–Laird tau^2 and random-effects estimate, separately for the four survival trials and the three
pCR comparisons; leave-one-trial-out for the seven-comparison combination and for the primary.

(c) IMmotion150 entered as atezolizumab alone against sunitinib only, in the four-survival-trial model and
the seven-comparison combination (sensitivity).

## 8. Reporting

Wording of the randomized results follows what the design supports: an interaction with
checkpoint-containing regimens, not a property of checkpoint blockade alone; the first four comparisons
are exploratory evidence (their plan was deposited after the per-trial results were known), the three
further trials a replication under a plan deposited before analysis, and the seven-comparison P a
combination, not a confirmation. The TCGA-STAD patient counts (389 and 383) are reconciled.
