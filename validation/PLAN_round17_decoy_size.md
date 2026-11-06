# Plan, round 17: does the size of a decoy differ from its target's, and does that move the estimate? (2026-11-06)

Written and deposited before any sensitivity run. Counts of round 17 S1 for the response library were read first, while
checking a note from the figure redraw: 167 of 623 cohort tests have a target measured at a size other than its decoys'.

Why it can happen: a decoy has one size per target, the median over cohorts of the target's measured genes (decoyfdr R/decoys.R),
and is drawn from genes measured in every cohort, so its size is fixed. The target loses genes on some platforms, so its measured
size moves by cohort. Each decoy is read against its target's bank in that cohort, which has the target's measured size. The text says
decoys are grown "to the target's size".

## Analyses
S1. Per library (response; TCGA OS, DSS, PFI; breast PC1): cohort tests where the sizes differ, by direction, and the quantiles of target
size over decoy size.
S2. Primary. Response, folded and directed, co-expressed and random-gene decoys. `decoyfdr::call_units(units, exclude = ...)` drops a target's
and its decoys' tests in every cohort where the target's measured size differs from its decoys'. Then `decoy_fdr` with its defaults
(B = 6000, seed 20261019). A unit left with fewer than 6 cohorts gets P = 1, the package's rule. Compared with the published estimates.
S3. The same with a tolerance: drop only tests where the sizes differ by more than 10% of the decoys' size.
S4. S2 and S3 for TCGA OS, DSS, PFI and breast PC1.
Control. With nothing excluded the script must return the published calls and estimates exactly. If it does not, stop.

Dropping tests lowers power for targets and decoys alike, so calls can fall without any bias having been present. The comparison that matters
is the estimate and the calls that remain.

## Decision rule, fixed now
- If the primary folded and directed estimates under S2 and S3 each move by 5 percentage points or less, or stay inside the published 95% interval
  over signatures, the design stands: the wording is corrected and the sensitivity is reported in Methods.
- Otherwise it goes to the authors before anything else, because decoys would then need regrowing per cohort.

Nothing else changes. Output: icb/v22_size/. Script: applications/v22_size_check.R. Uses the frozen package install, no new decoys.
