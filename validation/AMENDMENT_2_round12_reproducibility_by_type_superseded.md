# Amendment 2 to PLAN_round12_reproducibility_by_type.md

Written 2026-11-01, after the rerun under Amendment 1 and before any of its results beyond the checks was read (the
analysis script stopped at the check). Its SHA-256 is recorded in AMENDMENT_2_round12_reproducibility_by_type.sha256 and
it is deposited before the analysis is run.

## What happened

In the rerun, the folded C of every target and decoy matched the stored units in all 31 cohorts (to 6e-16), and the
whole-cohort rho matched rounds 10 and 11 to 1e-10 in every (signature, cohort) cell but one: APM (Wang) in TCGA COAD,
where the rerun gives 1/640 more than round 10. There the target's C equals 39/77 and one decoy's 38/77, so their
folded values are equal; the stored units, written to 15 significant digits, hold them one unit apart
in the 15th digit (0.506493506493506 against 0.506493506493507), and round 10 counted the decoy above the target. The
rerun, under the tie rule of Amendment 1, counts them tied, which is the correct reading. The inflamed family's COAD
value moves by 0.00008 (0.2891 against 0.2890); no reported number changes at the precision reported.

## Change

The check is accepted when every mismatch is of this kind: a cell whose stored target and decoy folded values differ by
less than 1e-13 without being equal. The script verifies that for each mismatched cell and stops otherwise. Round 12
uses its own whole-cohort values throughout; the round-10 numbers in the manuscript are left as reported, and the
discrepancy is stated in Additional file 1, Note S2.
