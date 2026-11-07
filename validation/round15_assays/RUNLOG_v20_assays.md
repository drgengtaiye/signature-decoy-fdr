# Run log, applications/v20_assays.R (round 15 part D)

- 2026-11-06, run 1 (12 cores): BCI finished in all 11 cohorts; the run then stalled at the start of EndoPredict for 6.5 hours and was stopped. Cause: the script held
  the TCGA pair in a variable `P`, which is also the TCGA data path that `applications/common.R` sets and `tcga_cohorts()` reads. The second `load_pair()` call
  therefore passed a list of expression matrices to `file.path()`. No result was written. Log kept as `icb/v20_round15/logs/v20_assays_run1_stalled.log`.
- 2026-11-06, run 2: the pair renamed `PR`; nothing else changed (same tables, seeds, cores). The other round-13 and round-15 scripts call `load_pair()` once per
  process and were not affected.
