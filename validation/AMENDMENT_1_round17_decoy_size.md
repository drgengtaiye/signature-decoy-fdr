# Amendment 1 to round 17 (2026-11-06): the control

The plan's control asked that the script, with nothing excluded, return the published estimates exactly. It held for response, TCGA OS and TCGA DSS and
failed once: TCGA PFI, folded, random-gene decoys, 10.885% against the published 10.833%. The run stopped, as the plan says.

Cause: the published calls were computed from the fit in memory; this script reads the stored units table, written to 15 significant digits. A Delta that
is about 1e-16 in memory is exactly 0 after that round trip, and an exact zero is dropped by the signed-rank test where a tiny value takes the
lowest rank. Near-ties between cohorts behave the same way. Discrete AUC values make both common enough to see. In PFI all 37 units whose P differs have such a
Delta, and one of 4,160 random-gene calls changes.

Change: the baseline for S2-S4 is the run with nothing excluded from the stored units. The control now passes if each baseline estimate is within 0.1 percentage
point of the published one, and any difference is printed. No published number changes. The decision rule stands, read against the published estimates.
