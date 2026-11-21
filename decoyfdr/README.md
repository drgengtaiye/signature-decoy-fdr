# decoyfdr

Matched decoy modules for the false discovery rate of a claim about a library of gene-expression signatures.

A claim such as "k of these 52 signatures predict response" is a claim about a library. `decoyfdr` estimates
how many of the k calls are false with the target-decoy design of peptide identification: every signature (the
target) receives decoy modules of genes that are in no signature of the library, matched to it on size, scoring
engine, expression-decile composition and internal co-expression in every cohort; targets and decoys go through
the identical test against random gene sets of the target's size scored by its engine; and the share of each
target's decoys that is called estimates the probability that a structure-matched module of non-signature genes
would be called. The estimated FDR is the sum of those shares over targets divided by the number of targets
called. `calibrate_null()` checks the estimate where the truth is known, by permuting outcomes within cohorts.

For one prognostic risk model or many, `check_models()` asks the same question of each validation: how often does a
model of the same structure, built from genes unrelated to it, pass the same test in the same cohort? See
"Checking prognostic risk models" below.

## Install

```r
install.packages("remotes")
remotes::install_github("drgengtaiye/signature-decoy-fdr", subdir = "decoyfdr")
```

Rcpp, matrixStats and survival are required (a C++ compiler is needed to build the package). GSVA (Bioconductor)
is needed only for the `"gsva"` and `"ssgsea"` engines, RSpectra speeds up `"pc1"` for large sets, and testthat,
knitr and rmarkdown are used by the tests and the vignettes (building the vignettes needs pandoc).

## Use

```r
library(decoyfdr)
# one dfdr_cohort per cohort: expression (genes x samples) and a 0/1 outcome or a data frame of time and event
cohorts <- list(dfdr_cohort(expr1, response1, "cohort1"), dfdr_cohort(expr2, response2, "cohort2"), ...)
lib <- dfdr_library(sets,                      # named list of gene vectors
                    engine = "gsva",           # "mean", "wmean", "gsva", "ssgsea" or "pc1", per set
                    direction = c(sigA = 1))   # +1: high score declared favorable; -1 unfavorable
res <- run_decoy_fdr(cohorts, lib, n_coex = 32, n_rand = 8, n_sets = 10, n_draw = 1000, keep_draws = 200,
                     cores = 8)
res$fdr                                         # calls and estimated FDR, folded and directed
cal <- calibrate_null(res$fit, B = 200)         # estimated / realized false calls under permuted outcomes
```

`res$fdr` gives, per statistic and decoy class, the signatures called, the estimate pooled over the decoy sets
(`fdr`, with a bootstrap interval over targets or, with `clusters =`, over clusters of correlated signatures), the
range of the per-set estimates (`fdr_set_min`, `fdr_set_max`) and the pi0-corrected estimate (`pi0`, `fdr_pi0`).
Results do not depend on `cores`: every target and cohort draws from its own L'Ecuyer-CMRG stream.

A unit (target or decoy) is called when the two-sided Wilcoxon signed-rank test over cohorts of its statistic
minus the median of its target's null bank gives P < 0.05 with a positive median; at least six cohorts are
needed. The statistic is the AUC for a binary outcome and Harrell's C for survival, folded (direction chosen
after the fact) or in the declared direction.

## Checking prognostic risk models

A prognostic risk model is a set of genes with Cox coefficients, fitted in one cohort (often TCGA) and "validated" in
another when its high-score group fares worse. `check_models()` gives every model of three or more genes 32
co-expressed and 8 random-gene decoys (gene sets of its size and expression level, the co-expressed ones also of
its internal co-expression, from genes in no model), each carrying the model's coefficients in a random order and
the direction of its own association with survival in the training cohort, and tests model and decoys in each
validation cohort under one rule: the score split at its median, log-rank P < 0.05, high group worse.

```r
library(decoyfdr)
tcga  <- dfdr_cohort(expr_tcga, data.frame(time = os_time, event = os_event, row.names = samples), "TCGA_COAD")
valid <- dfdr_cohort(expr_gse, clin_gse, "GSE39582")            # expression: genes x samples, same gene symbols
model <- c(GENE1 = 0.21, GENE2 = -0.35, GENE3 = 0.12, GENE4 = 0.40)  # Cox coefficients as printed
chk <- check_models(model, train = tcga, validate = valid, seed = 1, cores = 4)
chk                     # validated?, share of the model's decoys validated, rho
cohort_fvr("GSE39582")  # how often null models fitted in TCGA validated in this cohort in the paper
```

For each model and validation cohort, `chk$units` gives the model's directed C-index and log-rank P, whether it
validated, the share of its co-expressed (`d_coex`) and random-gene (`d_rand`) decoys that validated, and `rho`, the
share of co-expressed decoys it beats. Over a set of models, `validation_fdr(chk)` sums `d_coex` and divides by the
number validated: the estimated share of the validations that are false. For one model, `d_coex` says how often its
structure alone validates there; `rho` ranks it among its decoys but does not sort true from false (in the paper,
validated models with rho >= 0.95 replicated in another independent cohort in 30.5% of pairs, the others in 21.9%).
Decoys are matched in the training and validation cohorts by default (`match =`); give `exclude =` the gene set the
model was selected from (for instance a ferroptosis list) to keep it out of the decoys, as the paper did.

`cohort_fvr()` returns the table of the paper's 41 validation cohorts (43 rows, two of which repeat the patients of
others; 11 cancer types): for each, the share of null models validated there, where a null model is the usual
pipeline (univariate Cox screen, then LASSO Cox) run in TCGA on a pool of genes outside every published theme and
model. The rates run from 0% to 85%.

## What the estimate estimates

A signature counts as a false call when it is called no more often than a module of non-signature genes with its
structure. The plain estimate sums the decoys' call rate over every target, true ones included, so it is the
estimate with pi0 = 1 and an upper bound for that target when the decoys are matched; under permuted outcomes,
where every target is null, it is unbiased. The pi0-corrected estimate multiplies it by the share of targets whose
statistic does not stand out from their own decoys' (Storey's estimator with the decoys as the null reference).
Averaging over several independent decoy sets (`n_sets`, ten by default) reduces the estimate's dependence on one
draw of decoys.

## Reproducing the paper's numbers

`tests/testthat/test-paper.R` reproduces the 21 folded and 27 directed calls and the FDRs of 17.3% and 7.4% of
the response benchmark (its first decoy set) from the per-cohort statistics in
`inst/extdata/response_units_N0.csv.gz`. The scripts in `../applications/` run the package on the response cohorts,
on 20 TCGA cancer types, on breast cancer cohorts and on the MSigDB Hallmark library, with the seeds given there.
`check_models()` reproduces the paper's audit of its whole literature search (728 published models, 11 cancer types) from the same inputs
and seeds: `../applications/v24_check_models_repro.R` runs it, `../decoys/v24_check_models_repro.py` compares the
output with the paper's tables, and `../validation/TOOL_check_models_acceptance.md` states the criteria.
