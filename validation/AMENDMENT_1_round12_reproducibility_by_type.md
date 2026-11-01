# Amendment 1 to PLAN_round12_reproducibility_by_type.md

Written 2026-11-01, after the TCGA part of section A had been computed once and before any external split was
analyzed. Its SHA-256 is recorded in AMENDMENT_1_round12_reproducibility_by_type.sha256 and it is deposited before the
rerun.

## What happened

The check of section A failed narrowly. The folded C of every target and decoy, recomputed on whole TCGA cohorts from
the stored scores, matched the stored units to 6e-16, but the inflamed family's whole-cohort rho differed from round 10
(icb/v15_content/atlas_by_type.csv) by up to 0.0004, above the 1e-10 required. The cause is ties. Folding maps C to
max(C, 1 - C); when a target's C lies below one half and a decoy's lies the same distance above it, the two folded
values are mathematically equal but can differ in the last binary place. Round 10 read them from the stored units
(15 significant digits), where they are equal and count one half; the first run of applications/v18_split_half.R
compared the unrounded values and counted them as one below the other. 24 of the 1,040 (signature, type) cells
changed, by 1/640 to 5/640.

Before the cause was found, the family values of that first run were summarized once: for the inflamed family, the
median split-half Pearson r over the 20 types was 0.305 (Spearman-Brown 0.47) and the variance-based reliability 0.67;
confounder 0.55 and 0.62; other 0.62 and 0.70. The external run was cancelled during its splits, before any split
output was written; only its per-patient scores had been saved (change 3).

## Changes

1. Ties. In rho, a decoy's folded C counts as tied with the target's when the two differ by less than 1e-12 (distinct
   values of C differ by at least one half over the number of comparable pairs, more than 1e-8 in these cohorts).
   The check of section A is unchanged and must now pass to 1e-10; the same rule applies in section B.
2. Computation. Splits are processed in blocks of 25 per cohort, in parallel, for speed. The splits themselves are
   unchanged (drawn serially in cohort order from seeds 20261105 and 20261106 before any computation).
3. External scores. The per-patient scores computed by the first external run (icb/v18_split/external_scores.rds,
   written before the run was cancelled) are reused; the check of section B against round 11 applies to them as
   planned.

Nothing else in the plan changes. Both parts are rerun from the start with these rules, and the results reported are
those of the rerun.
