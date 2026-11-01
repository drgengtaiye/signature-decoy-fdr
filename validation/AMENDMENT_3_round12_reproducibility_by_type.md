# Amendment 3 to PLAN_round12_reproducibility_by_type.md: a sensitivity analysis added after the results

Written 2026-11-01, after every pre-specified result of sections A to E had been computed and read
(icb/v18_split/reliability.csv, transport_test.csv, simulation_round11_power.csv, sd_against_events.csv,
signature_level.csv). It adds a post hoc sensitivity analysis and changes no pre-specified analysis; its SHA-256 is
recorded in AMENDMENT_3_round12_reproducibility_by_type.sha256 and it is deposited before the analysis is run.

## Why

Section A2 estimates a type value's sampling SD at full size as SD over splits of (A - B) / 2, which is exact for a
statistic that is an average over patients, where the two halves move in opposite directions from the whole. The
family value of rho is not: it ranks folded C values, and folding maps a C below one half and one above it to the same
side. Where the inflamed signatures' C lies near one half, the two halves of a split therefore move together rather
than apart, and (A - B) / 2 understates the noise. The results show it: across the 500 splits the halves' values
correlated positively in TCGA STAD (0.75), LUSC (0.69), COAD (0.19) and PAAD (0.13) and in the external cervical
(0.54), head and neck (0.69) and lung adenocarcinoma (0.70) cohorts, all types whose value lies well below 0.5 (not
every such type: TCGA KIRC -0.64); they were near zero in TCGA GBM (-0.08) and the external sarcoma cohort (0.09), and
-0.41 to -0.95 elsewhere; in TCGA STAD SD(A - B) / 2 was 0.053 while either half varied over splits with SD 0.15. The
quantities built on SD_full (A2, C, D, E1) may therefore be optimistic: reliabilities too high, Q too large, power too
high. The split-half correlations of A1 and B and the single-signature rates of E2 do not use SD_full.

## Sensitivity analysis (post hoc)

1. For each of the 31 cohorts, 200 nonparametric bootstrap resamples of its patients (seeds 20261107 for TCGA and
   20261108 for the external cohorts, drawn serially in cohort order before any computation), with the same scores,
   decoys and rho (tie rule of Amendment 1); SD_boot of each family value over the resamples.
2. With SD_boot in place of SD_full: the variance-based reliability of A2 and B, the test of C (Q, P, I2), the
   simulation of D (with R_TCGA recomputed from SD_boot) and the fit of E1.
3. Reported beside the pre-specified results, labeled post hoc, with the reason above. The pre-specified results are
   not replaced, and no conclusion is drawn from the sensitivity analysis that the pre-specified results contradict
   without saying so.

Scripts: applications/v18_split_half.R (apps tcga_boot and external_boot) -> icb/v18_split/rho_boot_*.csv.gz;
decoys/v18_reproducibility.py (section S) -> icb/v18_split/*_boot.csv.
