# Analysis plan, round 9: the decoy estimate beside the standard false discovery rate estimates

Written 2026-10-30, after an internal review of manuscript v10 and before any analysis below was run. Its SHA-256 is
recorded in PLAN_round9_standard_fdr_comparison.sha256 and the file is deposited before the first analysis starts.

## What is known when this is written

* Every result of manuscript v10. The review's point: the paper never sets its estimate beside the estimates an
  analyst would compute anyway from the signatures' own P values (Benjamini-Hochberg, Storey) or from permuted
  outcomes (as in SAM), so a reader cannot tell what the decoys add.
* Inputs already on disk, inspected for their format only: the per-target calls with two-sided signed-rank P values
  (icb/v12_pkg/{response,tcga_OS,breast_pc1,hallmark_OS}_calls.csv.gz); the global-null permutations
  (icb/v12_truth/null_pi0.csv, icb/v12_pkg/*_globalnull_per_perm.csv); the planted-signal replicates
  (icb/v12_truth/signal_*.csv, tcga_planted_*.csv).
* Numbers already seen: the mean number of calls per library under the global null (response 1.926 folded and 1.032
  directed; TCGA overall survival 2.655 and 1.305; breast PC1 1.22; Hallmark 1.525) and the mean number of calls per
  planted setting (round-7 summaries). The nominal estimate below is m x t divided by such numbers, so its global-null
  ratio and its value in each planted setting follow arithmetically from known numbers. These comparisons are
  therefore descriptive, not tests, and are reported whatever they show. Nothing below has been computed.

## A. Analyses

Notation. m targets; the call of target s is P_s < 0.05 (two-sided signed-rank over cohorts) with a positive median
Delta, i.e. the one-sided P value p_s = P_s / 2 if the median Delta > 0 and 1 - P_s / 2 otherwise is below t = 0.025.
R = number of calls. The decoy estimates are those already reported and are not recomputed.

1. Real data. For the response analysis (folded and directed, null N0), TCGA overall survival (folded), breast
   cancer PC1 with every cohort (folded; the analysis on which its global null was run) and the Hallmark sets for
   overall survival (folded):
   a. nominal estimate F_nom = m t / R (Storey's FDR(t) with pi0 = 1, the rate the Benjamini-Hochberg procedure
      implies at threshold t);
   b. Storey estimate F_St = pi0_St m t / R, pi0_St = #{p_s > 0.5} / (0.5 m), capped at 1;
   c. permutation estimate F_perm = Rbar_perm / R, Rbar_perm the mean number of calls per library under the stored
      global-null permutations of the same analysis (as SAM estimates the FDR); and pi0_St x F_perm;
   d. the number of targets with Benjamini-Hochberg adjusted p_s at most 0.05 and at most 0.10, and the BH level at
      which exactly the R calls would be made (m p_(R) / R);
   e. beside them, the decoy estimate and its pi0-corrected form as reported.
   Output: icb/v14_standard/standard_real.csv; Table 3 of the manuscript.

2. Global null. For each analysis of 1, the ratio of the nominal expectation m t to the realized mean number of
   calls under the stored permutations, with a 95% interval from 4,000 bootstrap resamples of the permutations
   (seed 914000), beside the decoy estimate's reported ratio. The permutation estimate is exact here by construction
   and is so stated. Output: icb/v14_standard/standard_globalnull.csv.

3. Planted signal (16 response settings, 6 TCGA settings of round 7). Per setting, as the round-7 summaries do
   (ratios of totals over replicates): realized FDR = sum V / sum R under the call-probability truth (response also
   under the gene-content truth); nominal estimate = reps x m t / sum R; permutation estimate = reps x Rbar_perm /
   sum R where a global null of the same cohorts and statistic was run (response, 12 cohorts, folded and directed;
   TCGA overall survival, folded), not available for the four 6-cohort settings; the decoy estimates as reported.
   Reported: per estimate, the number of settings in which it falls below the realized FDR, under each truth; per
   replicate with calls (12 cohorts, folded, beta = 1, as Fig. 5j), the share in which each estimate is at or above
   the realized proportion. Output: icb/v14_standard/standard_planted.csv.

Script: signature-decoy-fdr/decoys/v14_standard_fdr.py. No package code changes; nothing reported elsewhere changes.

## B. Displays

* Figure 5 is reorganized. Kept: a (with the nominal expectation m t added to each histogram), b, d, e, the TCGA
  planted panel, the decomposition, the replicate-by-replicate scatter and the estimate-minus-realized panel. New:
  (i) the global-null ratio of analysis 2 for the nominal and decoy estimates, by analysis; (ii) the estimates of
  analysis 3 against the realized FDR, setting by setting. Moved to a new Additional file 1 Fig. S7, unchanged in
  content: the per-unit call rates, the number of true signatures, the interval coverage and pi0 under the global
  null.
* Table 3 of the manuscript: analysis 1.

## C. Text

Analysis results enter the Results (one paragraph), Table 1 (one row), Table 3, the Discussion and, if space allows
within 250 words, the abstract. The main text is shortened for reading without changing any number; panel citations
follow the reorganized Figure 5 and Fig. S7.
