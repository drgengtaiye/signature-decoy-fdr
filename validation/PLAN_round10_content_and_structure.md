# Analysis plan, round 10: gene content and co-expression structure, by endpoint and cancer type

Written 2026-10-31, after manuscript v12 and before any analysis below was run. Its SHA-256 is recorded in
PLAN_round10_content_and_structure.sha256 and the file is deposited before the first analysis starts.

## Question

The decoys split what a signature achieves over its null into what any module of non-signature genes with its
structure achieves (structure) and what its own genes add (content). The manuscript reports this split only as library
false discovery rates. This plan asks, from stored outputs, whether the split differs between response to checkpoint
blockade and survival in untreated patients, between signature families, and between cancer types.

## What is known when this is written

* Every result of manuscript v12, in particular: calls by family (response, folded: 16 of the 20 inflamed signatures
  and 0 of the 8 confounder signatures among 21 calls; TCGA overall survival: 5 of 20 and 5 of 8 among 13 calls, 5 of
  the 13 also response calls); the FDR estimates and pi0 (response folded 18%, pi0 0.23; TCGA overall survival 55%,
  pi0 0.42); the per-signature place among its decoys for response (Fig. 4h, signature level); CIN70 (Carter) above
  its decoys' 95th percentile in 9 of the 20 cancer types (Fig. 6l); 59% of signature x cancer-type gains positive
  (Fig. 6j); decoys' co-expression per cancer type (Fig. 6c, d).
* The pi0 inputs, icb/v12_pkg/*_fold_target_vs_decoys.csv (each target's signature-level P among its decoys), have
  been used for pi0 and Fig. 4h but not tabulated by family or contrasted between endpoints.
* Inputs inspected for their format only: icb/v12_pkg/{response,tcga_OS,tcga_DSS,tcga_PFI,breast_pc1,breast_mean,
  hallmark_OS}_units.csv.gz (per cohort and unit: fold_obs, fold_med, fold_p95, dir_obs, dir_med, dir_p95; a decoy
  row carries its target's null; decoy classes coex and rand, sets 1-10); icb/v5_signature_clusters.csv (family,
  cluster_fixed); the cohort cancer types of figstyle.ICB_TUM; icb/v12_pkg/breast_origin_excluded_tests.csv.
* Given the call counts by family, the expected direction of the primary contrast (inflamed signatures gaining
  content for response relative to survival, confounder signatures not) is anticipated. The contrast is therefore a
  pre-specified description, not a blind test, and is reported whatever it shows. Nothing below has been computed.

## A. Definitions

For unit u (a target or one of its decoys) in cohort or cancer type c: Delta(u, c) = fold_obs - fold_med, read
against the target's null as stored. Decoys are the target's co-expressed pool-A decoys over all ten sets (coex);
random-gene decoys (rand) are reported beside them.

* Structure S(s, c): the median Delta over the co-expressed decoys of target s in c.
* Content G(s, c) = Delta(s, c) - S(s, c).
* Place among decoys rho(s, c): the share of the co-expressed decoys of s whose Delta in c lies below the target's,
  ties counted one half. If the signature's genes add nothing to a module of its structure, rho is near uniform
  with mean 0.5. rho is free of the endpoint's scale (AUC against C-index) and of cohort size.
* Signature level, per setting E: rho_bar(s, E), the mean of rho over the cohorts of E; content call: the Wilcoxon
  signed-rank test of G(s, c) over the cohorts of E, two-sided P < 0.05 with a positive median (the manuscript's call
  rule applied to G; nominal one-sided rate 2.5%).
* Structure per setting: the share of co-expressed decoy x cohort tests with fold_obs above fold_p95 (5% nominal),
  overall and by the target's family; the median over targets and cohorts of S.

Settings: response (12 cohorts) and TCGA overall survival (20 types), primary; TCGA disease-specific survival and the
progression-free interval, secondary; breast PC1 (every cohort, and with the origin tests removed), breast mean z and
Hallmark overall survival, library-level only (the families do not apply). Folded throughout; directed for response
as a secondary analysis in Additional file 2.

## B. Analyses

1. Content by endpoint and family (primary). For each of the 52 signatures, C(s) = rho_bar(s, response) -
   rho_bar(s, TCGA OS). Primary contrast: mean C over the 20 inflamed signatures minus mean C over the 8 confounder
   signatures. P value: 100,000 random reassignments of the two family labels among the 28 signatures (seed 915000),
   two-sided. 95% interval: 10,000 resamples of the 13 correlation clusters (cluster_fixed) with replacement,
   recomputing the contrast from the signatures drawn (seed 915001). Reported beside it: the family means of rho_bar
   per setting (inflamed, other, confounder, all 52) with cluster intervals; the number of content calls per family
   and setting; for all 52, rho_bar(response) against rho_bar(TCGA OS) by the paired Wilcoxon test over signatures
   (descriptive, as signatures are correlated).

2. The same cancer types. TCGA restricted to the types of the response cohorts (SKCM for melanoma, KIRC for kidney,
   LUAD and LUSC for lung, BLCA for bladder, STAD for gastric): rho_bar(s, TCGA OS matched) and the primary contrast
   of 1 recomputed with it; per tumor type, the inflamed and confounder family means of rho in the response cohorts
   and in the matched TCGA types.

3. Pan-cancer atlas. The 52 x 20 matrix rho(s, t) for TCGA overall survival and the 52 x 12 matrix for response.
   Per cancer type t: the inflamed and confounder family means of rho; the claim "these signatures are prognostic in
   t", with calls R_t = targets with fold_obs above fold_p95 in t, and its single-type decoy estimate FDR_t = sum over
   targets of the share of their co-expressed decoys with fold_obs above fold_p95 in t, divided by R_t (reported as
   computed, above 100% included); the same per response cohort. Descriptive; no claim about any single type is
   tested.

4. Libraries. For each setting of A: mean rho_bar over all targets and over called targets, the median S, the decoy
   clearance share and the reported pi0; side by side, as the library-level summary of content against structure.

Script: signature-decoy-fdr/decoys/v15_content_structure.py; outputs icb/v15_content/. No package code changes;
nothing reported elsewhere changes.

## C. Displays and text

* A new main figure (Figure 8): the split of each signature's gain into structure and content for response and
  overall survival; rho_bar by signature between the two endpoints, by family; the family means with intervals for
  every setting of A and for the matched types; the atlas of 3 with the per-type estimates; the library summary of 4.
* A Results subsection placed last, a Methods subsection, a Discussion paragraph, Additional file 2 sheets and, if
  space allows within 250 words, one abstract sentence. Titles and wording state what the analyses show, whatever it
  is; numbers already reported do not change.
