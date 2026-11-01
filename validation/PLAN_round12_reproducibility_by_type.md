# Analysis plan, round 12: how reproducible is a cancer type's content estimate?

Written 2026-11-01, after manuscript v14 and before any computation below. Its SHA-256 is recorded in
PLAN_round12_reproducibility_by_type.sha256 and the file is deposited before the first computation.

## Question

Round 11 (PLAN_round11_external_replication.md) failed its primary test: over 11 cancer types, the inflamed family's
mean rho (a signature's place among its co-expressed decoys, for survival) in independent cohorts did not track the
TCGA value of the type (Spearman -0.28, one-sided P = 0.81). Two explanations differ in what they imply. (i) The value
of a type is not reproducible even within one cohort at the sample sizes involved, so no replication was to be
expected and a type-level claim from one cohort of this size says little. (ii) It is reproducible within TCGA but does
not carry over to other cohorts (platform, endpoint, population, decoy pool). This plan separates the two by splitting
each cohort at random into halves, which measures the patient-sampling noise of every estimate, and asks whether the
TCGA-external discordance of round 11 exceeds that noise.

## What is known when this is written

* Every result of manuscript v14 and rounds 1 to 11, including the per-type values of the inflamed, confounder and
  other families in TCGA (icb/v15_content/atlas_by_type.csv) and in the 11 external survival cohorts
  (icb/v17_external/per_cohort.csv).
* The TCGA overall-survival cohorts and their sizes (patients / events): BLCA 405/177, BRCA 1081/151, CESC 290/71,
  COAD 279/69, ESCA 184/77, GBM 153/121, HNSC 519/220, KIRC 531/175, KIRP 287/44, LGG 511/125, LIHC 365/130,
  LUAD 502/182, LUSC 494/212, MESO 85/73, OV 303/183, PAAD 177/93, SARC 259/98, SKCM 103/29, STAD 389/156,
  UCEC 174/32. Several types have few events (SKCM 29, UCEC 32, KIRP 44), so large sampling noise in those types is
  foreseeable; how large, and whether it accounts for round 11, is not known.
* The structure of the stored round-7 TCGA fit (icb/v12_pkg/tcga_OS_fit.rds): it holds each cohort's per-patient
  scores of the 52 targets and all 20,800 decoys. Nothing has been computed from them on any subset of patients.

## A. Split halves in TCGA

For each of the 20 TCGA types, 500 random splits into two halves, stratified by event status (within the patients with
an event and within the censored patients, a random permutation; the first floor(n/2) to half A, the rest to half B;
seed 20261105). Scores are the stored per-patient scores of the round-7 fit (targets by their engines, the ten sets of
32 co-expressed decoys; scores use no outcome, so they are not recomputed per half). In each half: Harrell's C of every
target and co-expressed decoy (decoyfdr::cindex_rows), folded; rho = share of the target's 320 co-expressed decoys whose
folded C lies below the target's (ties count one half), as in round 10; family means over the inflamed (20), confounder
(8) and other (24) signatures. Check first: on the whole cohort these computations must reproduce the stored units
(icb/v12_pkg/tcga_OS_units.csv.gz, fold_obs) and atlas_by_type.csv (rho by family) to 1e-10.

1. Primary. Split-half reliability of the inflamed family's type values: for each split, the Pearson correlation r over
   the 20 types between halves A and B; the median over splits, carried to full size by Spearman-Brown,
   R = 2r / (1 + r). Interval: 2,000 resamples of the 20 types (seed 916100), the median over splits recomputed in each.
   Pre-stated reading: R of 0.7 or more, the TCGA type pattern is a stable property of TCGA and explanation (ii) is
   favored; R below 0.5, the type values are dominated by sampling noise at TCGA's sizes and explanation (i) is favored;
   in between, both contribute. The Spearman form is reported beside it.
2. Sampling SD of each type's full-cohort value: SD_full = SD over splits of (A - B) / 2. With it, a second estimate
   of reliability, 1 - mean(SD_full^2) / var(type values) (variance over the 20 types of the full-cohort values).
3. The same two quantities for the confounder and other families (descriptive).

## B. Split halves in the external cohorts

The same in each of the 11 external survival cohorts of round 11 (seed 20261106), with per-patient scores recomputed
from the stored round-11 decoys (icb/v16_external/external_decoys.rds) by the steps of decoyfdr::fit_units, with the
random draws of the null bank made but not scored so that the weights given to the weighted-mean decoys are those of
round 11. Check: on the whole cohort the folded C of every target and decoy must equal round 11
(icb/v16_external/external_surv_units.csv.gz) to 1e-10; if it does not, fit_units is rerun unchanged with the scores
kept, and the check repeated. External reliability across the 11 types (as A1, with 11 types) and SD_full per cohort.

## C. Does the round-11 discordance exceed sampling noise?

For the 11 types of round 11, D_t = external value - TCGA value of the inflamed family, with v_t = SD_full,TCGA^2 +
SD_full,ext^2. Q = sum over types of (D_t - weighted mean D)^2 / v_t (weights 1 / v_t), against chi-square with 10
degrees of freedom; the weighted mean absorbs an overall offset between TCGA and the external cohorts (different decoy
pools and platforms). P < 0.05: the discordance exceeds what sampling produces (heterogeneity between cohorts); P of
0.05 or more: it does not. I^2 = max(0, (Q - 10) / Q) is reported. The same for the confounder family (descriptive).

## D. What the round-11 test could detect (descriptive)

Simulation under transportability (10,000 draws, seed 916101): true type values theta_t = m + R_TCGA (TCGA value - m)
for the 11 types (m their mean; R_TCGA from A2, the variance-based reliability, bounded to [0, 1]); a TCGA value drawn
as theta_t + N(0, SD_full,TCGA^2) and an external value as theta_t + N(0, SD_full,ext^2); the Spearman correlation over
the 11 types. Reported: its median, the share of draws at or below -0.28 (round 11's value), and power, the share
above the one-sided 5% critical value of the Spearman with 11 untied values (the 95th percentile of 100,000
permutations, seed 916102).

## E. Sample size and single-signature claims (descriptive)

1. log SD_full against log events over the 31 cohorts (20 TCGA, 11 external), least squares; the slope and the number
   of events at which the fitted SD of a type value is 0.10 and 0.05.
2. Signature level: for each signature and cohort, rho in half A and half B. The share of (signature, cohort) cells
   with rho of 0.95 or more in half A that also have it in half B, averaged over splits and both directions, against
   the share of all cells with rho of 0.95 or more (the rate expected if the halves were unrelated); TCGA and external
   separately. The Pearson correlation over cells between the halves, median over splits.

Every result is reported whatever it shows. Halves have half the patients of the full cohorts, which the Spearman-Brown
step and the SD_full formula address for the type values but not for E2, which is stated as a half-cohort rate.

## F. Scripts, outputs, displays

applications/v18_split_half.R (A, B: per-split rho tables) -> icb/v18_split/; decoys/v18_reproducibility.py (A to E)
-> icb/v18_split/. A supplementary figure (Additional file 1: Fig. S9), a Results paragraph after the external
cohorts, a Methods paragraph, Additional file 2 sheets, and at most one clause in the abstract. Numbers already
reported do not change.
