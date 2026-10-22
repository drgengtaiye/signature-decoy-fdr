# Plan, round 8: descriptive displays for the redrawn figures

Written 2026-10-23. Displays 1-3 and the data of display 12 were drawn earlier the same day, before this plan was
written; the others are computed after its deposit. A first version deposited minutes earlier said that none had
been computed; it is kept as PLAN_round8_figure_displays_superseded.md. The SHA-256 of this version is in
PLAN_round8_figure_displays.sha256.

## Scope

The round-8 figures redraw the round-7 results. No estimate, interval, call, P value or count reported in the text
changes, and no display below enters the abstract, the conclusions or any test. Every display is drawn from outputs
that already exist (score matrices, unit tables, call tables, replicate tables, fits) and is labelled in its legend
as descriptive.

## Displays

1. ROC curve of Chemokine (Messina) in the Gide cohort with its 32 co-expressed decoys and the pointwise 5th-95th
   percentile of its 1,000-draw null bank (first decoy set; folded direction).
2. Standard deviation of each matched null against cohort size, with the Mann-Whitney standard deviation of a folded
   AUC as reference.
3. Per cohort, the share of tests clearing a fixed 48-gene null against the matched null.
4. Eigenvalues of the 52-signature correlation matrix in each checkpoint cohort (the quantity from which the Li-Ji
   effective number, already reported, is computed).
5. Similarity of the correlation structure between cohorts: the correlation, over the 1,326 signature pairs, of the
   pairwise score correlations in two cohorts.
6. Gene overlap (Jaccard index) against score correlation for every signature pair.
7. Genes shared by the most signatures.
8. Signed-rank P against median gain for every signature (response; first decoy set, N0).
9. The number of signatures called and the estimated FDR as the per-unit P threshold varies from 0.001 to 0.2
   (folded, ten sets pooled); the pre-specified threshold, 0.05, stays the only one reported.
10. Quantile-quantile plots of the per-unit signed-rank P of targets and of decoys against the uniform.
11. Replicate-level realized against estimated FDR in the planted-signal simulations.
12. Kaplan-Meier curves, split at the median, of a called signature and of its median co-expressed decoy in one
    cohort (TCGA KIRP, CIN70 (Carter), overall survival; METABRIC, SOTIRIOU-GGI, its endpoint). A log-rank P is
    printed for orientation only.
13. Heatmaps of the median gain of every signature in every TCGA cancer type and every breast cohort.

## Rules

The examples (1, 12) were chosen by looking at effect sizes: Chemokine (Messina) in Gide in round 7, CIN70
(Carter) in KIRP as the called signature with the largest overall-survival gain in any cancer type, SOTIRIOU-GGI
in METABRIC among the largest gains in the largest breast cohort. They illustrate the procedure and are not
representative; the legends say so. They are fixed here and not replaced by others after drawing. If a display shows something at odds with
a reported result, the display is kept and the discrepancy is reported to the authors.
