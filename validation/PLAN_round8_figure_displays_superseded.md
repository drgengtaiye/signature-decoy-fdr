# Plan, round 8: descriptive displays for the redrawn figures

Written 2026-10-23, before any of the displays below is computed. Its SHA-256 is recorded in
PLAN_round8_figure_displays.sha256 and the file is deposited before the first display is drawn.

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

Examples (1, 12) are fixed here and not replaced by others after drawing. If a display shows something at odds with
a reported result, the display is kept and the discrepancy is reported to the authors.
