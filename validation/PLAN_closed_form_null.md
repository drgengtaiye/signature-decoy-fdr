# Analysis plan: a closed-form matched null

Written before the approximation was fitted. Deposited to a private cloud object store, whose object
timestamp is the external record.

## The question

Every claim in this paper is made against a null built by simulation: 200 random gene sets per size
rung, per cohort, per scoring engine. That is the expensive part of the protocol and the main reason a
reader cannot apply it quickly to their own cohort. This plan asks whether the null can be predicted
instead of drawn, from statistics of the cohort that cost one pass over the expression matrix, and how
much accuracy is lost.

## Data

The null rungs already computed for the benchmark (`icb/native_response_matched_v2.csv`): for each
cohort, engine and size rung, the median, mean, 95th percentile and standard deviation of the folded
AUC over 200 random gene sets. No outcome-based quantity is fitted here that is not already public in
that file.

## The approximation

For a score that averages gene-wise z-scores, the unfolded AUC of a random set of size k has a
standard deviation that falls with k but not to zero, because genes are correlated and a large random
set converges on the cohort's dominant axis of variation. The approximation is

    sigma(k) = sigma_inf + (sigma_1 - sigma_inf) / sqrt(k_eff),    k_eff = k / (1 + (k - 1) * rho),

with sigma_1 the standard deviation of single-gene folded AUC in that cohort, rho the mean absolute
correlation between random gene pairs, and sigma_inf a cohort-level floor. Folded AUC follows as
median and mean of a half-normal at 0.5. Two free parameters are fitted across cohorts: how sigma_inf
scales with sigma_1, and one multiplicative constant per engine.

## Test

Leave-one-cohort-out: the parameters are fitted on the other cohorts and the held-out cohort's null is
predicted at every rung from its own sigma_1 and rho alone. Reported are the relative error in the
predicted standard deviation and 95th percentile, and, for the 52 signatures, how often the predicted
null changes whether a signature clears its null at the 95th percentile. The approximation is useful if
the median relative error is below 10% and fewer than 5% of the 52 signature calls change; it is
reported as inadequate otherwise, and the package keeps drawing the null.

## What this cannot show

The approximation is fitted on checkpoint cohorts of 24 to 348 patients with bulk expression; nothing
here says how it behaves in much larger cohorts, in single-cell data, or for scoring engines other
than the three used here.
