# Acceptance test of `check_models()` (decoyfdr 0.3.0)

Written 7 November 2026, before the function was run on real data. This is a reproduction check of software, not a new
analysis: it asks whether the packaged function returns the numbers of round 15 part A (and the scores of part B5) from
the same inputs and seeds. It tests no hypothesis and adds no estimate to the paper.

## What the function does

`check_models()` packages the per-model audit of rounds 13 and 15 (applications/v20_published_models.R,
v20_cohort_scores.R; decoys/v20_audit.py, the rho of v20_cohorts_audit.py):

- a model's genes are kept if measured in every cohort the decoys are matched in; a model short of genes there keeps
  the rest when they are at least 80% of its printed genes and at least two (round 15 amendment 2);
- its linear score (coefficient x within-cohort z-score, summed) is declared to mark shorter survival when high;
- 32 co-expressed and 8 random-gene decoys per model of three or more genes, from `grow_decoys()` with the given seed;
  each decoy carries its model's coefficients in a random order (Mersenne-Twister, the same seed, decoy by decoy in the
  order of `grow_decoys()`) and the direction of its own association with survival in the training cohort;
- rule (L) in every validation cohort: median split, log-rank P < 0.05, high group worse;
- in a validation cohort outside the matching set, a model or decoy is scored when at least 80% of its genes are
  measured there, on those genes (round 13 A5, round 15 B5);
- per model and cohort: the share of co-expressed and of random-gene decoys validated, and rho, the share of
  co-expressed decoys whose directed C-index is below the model's (ties one half).

## Run

applications/v24_check_models_repro.R, for each of the 11 types of round 15: the frozen models of the type
(validation/published_models_r15.csv, `scored`), training cohort TCGA, decoys matched in TCGA and the round-11 external
cohort of the type, seed 20261202 + type index, exclude = the theme pools; validation cohorts = every built cohort of the
type (built_cohorts_r15.csv). The package is installed into a new library (icb/v24_tool/Rlib); the frozen install
icb/v12_pkg/Rlib6 is not touched.

## Criteria (all must hold)

1. Decoys: the same decoy names, in the same order, with the same genes (models_decoy_info_<TYPE>.csv.gz) and the same
   shuffled coefficients to 8 significant digits and directions (models_<TYPE>.csv.gz).
2. Part A, the external cohort and TCGA: for every model and decoy, `validated` identical; directed C-index and log-rank
   P equal within 1e-10 (absolute) and 1e-8 (relative).
3. Per model (audit_models_units.csv): d_coex, d_rand, n_coex and rho identical.
4. Estimates from the function's output equal those of audit_models.csv to the last digit: validated count, decoy
   share and FDR for the independent models (427) and all models with decoys (678). Bootstrap intervals are not
   compared (the package resamples with R's generator, the paper with numpy's).
5. Part B5 (cohort_scores_<TYPE>.csv.gz, published models and decoys): `scorable` and `validated` identical in every
   cohort; directed C-index within 1e-6 (part B read the decoys' coefficients rounded to 8 significant digits). One
   known exception, stated before the run: in BRCA_GSE96058, model BRCA_PMC11332736 was scored on 7 genes in part A
   (DUX4 not measured in TCGA after loading) and on 8 in part B; the function follows part A, so its C-index there
   differs from part B's (validated in both).

Any other difference fails the test and is reported, with its cause, before the function is released.
