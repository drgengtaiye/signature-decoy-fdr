# Amendment 2 to PLAN_round15_literature_cohorts_controls.md: a gene removed when a cohort is loaded

Written 2026-11-03, after the scoring of A3 ran for ten types and stopped for BRCA, and before any of its outputs was
summarized or read beyond the runs' completion lines and the decoy-matching lines of their logs. Its SHA-256 is recorded
in AMENDMENT_2_round15_literature_cohorts_controls.sha256 and it is deposited before BRCA is scored.

## What was seen

applications/v20_published_models.R stops when a frozen model holds a gene absent from the analysis matrices. The frozen
table (amendment 1) maps genes against the row names of the cohorts' files, as round 13 did; the loaders then remove
genes with non-finite or constant values. In BRCA one model, BRCA_PMC11332736 (8 genes), holds DUX4, which is in the
TCGA-BRCA file but removed by the TCGA loader. No model of the other ten types holds such a gene (their runs passed the
same check), so their outputs do not depend on this amendment.

## Change

A gene in a cohort's file but removed by its loader counts as unmeasured. A model keeps its other genes, with their
coefficients, if they are still at least 80% of its printed genes and at least two (the rule of round 13 A5);
otherwise it is not scored. BRCA_PMC11332736 is scored on 7 of its 8 genes. The frozen table is unchanged; the change is
in the scoring script, which prints every such case.
