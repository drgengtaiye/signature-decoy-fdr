# signature-decoy-fdr

Code for **Matched decoy modules estimate the false discovery rate of claims about gene-expression signatures**.

A claim that so many signatures of a library predict an outcome can carry a false discovery rate. Each signature
receives decoy modules of genes in no signature, matched to it on size, scoring engine, expression level and
internal co-expression in every cohort; signatures and decoys pass through the identical test against random gene
sets of the signature's size and engine, and the decoys' call rate estimates the false calls.

## The package

`decoyfdr/` is an installable R package implementing the procedure: cohort preparation, five scoring engines, matched
null banks, decoy growth in several independent sets, the call rule, the estimate (and its pi0-corrected form) with
target or cluster bootstraps, and the calibration of the estimate under permuted outcomes.

```r
install.packages("remotes")
remotes::install_github("drgengtaiye/signature-decoy-fdr", subdir = "decoyfdr")
```

It needs Rcpp, matrixStats and survival (and a C++ compiler); GSVA from Bioconductor is needed for the GSVA and ssGSEA
engines.

For prognostic risk models, `check_models()` tests any model given as Cox coefficients named by gene, with 32
co-expressed and 8 random-gene decoys matched to it, in each validation cohort under one rule (median split, log-rank
P < 0.05, high group worse), and reports the share of its decoys that validate and how many it beats;
`validation_fdr()` turns that into the estimated share of false validations over a set of models. `cohort_fvr()`
gives the false-validation rates of the 43 validation cohorts of the paper, also as the searchable page
`docs/validation_cohorts.html`. The function reproduces the paper's audit of 707 published models from the same
inputs and seeds (`validation/TOOL_check_models_result.md`).

See `decoyfdr/README.md` and the vignette. Its tests reproduce the response-benchmark calls and estimates of the
paper from the stored per-cohort statistics.

## The analyses

| Folder | Contents |
|---|---|
| `decoyfdr/` | the R package |
| `applications/` | the package run on the 12 response cohorts, on 20 TCGA cancer types, on breast cancer cohorts and on the MSigDB Hallmark library |
| `breast/` | download and preparation of the breast cancer cohorts and the Venet et al. signature library |
| `scripts/` | the response and prognostic benchmarks, the dependence between signatures, figures |
| `decoys/` | the target-decoy estimate on the response benchmark and its calibration (global null, planted signal, intervals) |
| `validation/` | the analysis plans with their SHA-256 hashes, when each was deposited (`PLAN_DEPOSITS.md`), and the provenance of the immune and breast signatures |
| `docs/` | the searchable page of validation-cohort false-validation rates |

`README.txt` gives the order in which the scripts run. The analysis scripts carry the absolute paths of the system
they were run on: lay out the inputs listed in `INPUTS.csv` under a local root and run
`python scripts/set_root.py /path/to/root` to rewrite them. `ENVIRONMENT.txt` lists the software versions used. All
input data are public.

License: MIT (see `LICENSE`). Citation metadata: `CITATION.cff`.
