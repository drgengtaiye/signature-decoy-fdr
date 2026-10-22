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

## Install

```r
install.packages("remotes")
remotes::install_github("drgengtaiye/signature-decoy-fdr", subdir = "decoyfdr")
```

Only Rcpp and matrixStats are required (a C++ compiler is needed to build the package). GSVA (Bioconductor) is
needed only for the `"gsva"` and `"ssgsea"` engines, RSpectra speeds up `"pc1"` for large sets, and survival,
testthat, knitr and rmarkdown are used by the tests and the vignette (building the vignette needs pandoc).

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
