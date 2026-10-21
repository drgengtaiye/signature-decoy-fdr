# Analysis plan: the treatment interaction pooled across trials

Written before the pooled model was fitted and deposited to a private cloud object store, whose object
timestamp is the external record.

## Why

The within-arm contrast and the per-trial interaction tests are limited by something this paper has
measured: the 17 immune signatures inside one trial are so correlated that they behave as one or two
independent tests, so no single trial can drive its P value low however extreme the count looks. The
way out is not another statistic on the same 17 but one pre-specified score tested across all the
randomized evidence at once. That removes the multiplicity entirely and raises the patient count from
263 to roughly 700.

## The score, already defined

The 12-gene inflamed score used in Figure 4d of the manuscript, defined there and not chosen here.
It is standardised within each trial before fitting, so a coefficient is per standard deviation within
that trial.

## Models

1. **Survival trials.** IMmotion150 and PRINCE pooled, a Cox model stratified by trial:
   `PFS ~ score + arm + score:arm`, strata(trial). The estimand is the score x arm coefficient, whose
   negative sign means the score predicts more benefit from the checkpoint arm.
2. **Pathological complete response trials.** The two I-SPY2 contrasts pooled, logistic regression
   with a trial indicator: `pCR ~ score + arm + score:arm + trial`. A positive coefficient means the
   score predicts more benefit from the checkpoint arm.
3. **All four.** The two model coefficients combined by Stouffer's method, weighted by the square root
   of the number of patients, in the benefit direction.

## Tests

- **Primary.** The score x arm coefficient in model 1, one-sided in the benefit direction, Wald test.
  Significant at P < 0.05.
- **Secondary, stated in advance.** Model 2; the Stouffer combination of the two; a permutation of the
  arm label within trial (2,000 draws) as a check on the Wald P of model 1; the same model with the
  three confounder-control scores in place of the inflamed score, where no interaction is expected;
  and each trial separately, which the manuscript already reports.

**What counts as a failure.** If the primary coefficient does not reach one-sided P < 0.05 the pooled
analysis is reported as null and the manuscript's per-trial results stand unchanged. No alternative
score, endpoint or subset will be substituted.

## What this cannot show

Two trials in two diseases with different comparators are pooled on the assumption that a treatment
interaction on the log hazard scale is common to them, which is an assumption and not a finding. The
patients are those already analysed here. A positive result means the score carries information about
who benefits from the checkpoint arm in these trials, not that it is ready to select patients.
