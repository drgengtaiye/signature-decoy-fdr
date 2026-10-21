# Analysis plan: survival under checkpoint blockade, against the same matched null

Written before the analysis was run. The question is whether immune signatures carry information about
survival *under* checkpoint blockade, measured with the same instrument already used for survival without it
(20 TCGA cancer types) and for response (12 checkpoint cohorts). The existing prognostic benchmark compares a
response endpoint under treatment with a survival endpoint without it; this analysis removes that mismatch.

## Cohorts

Every cohort of the ORCESTRA `clinical_icb` collection held locally (21 series) with, for the endpoint in
question, at least 25 patients and at least 10 events. ICB_VanDenEnde is excluded for the reason already fixed
in this study: its matrix is pseudogene-dominated and lacks CD8A, IFNG, CXCL9, GZMB and PTPRC. No cohort is
added or dropped after the results are seen.

## Endpoints

Overall survival is primary; progression-free survival is secondary. Both come from the collection's own
`survival_time_*` and `event_occurred_*` fields.

## Measurement

Signatures are scored by each publication's own engine, as everywhere else in this study. Performance is the
folded concordance index, as in the TCGA benchmark. The null is 200 random gene sets per size rung (4, 7, 12,
20, 32, 50, 80, 130, 260, 440 genes) per engine per cohort, drawn from that cohort's own measured genes, with
each signature compared against the rung nearest its own gene count under its own engine. Effects are reported
as the gain over the null median and as z = (observed − null median) / null standard deviation.

## Pre-specified analyses

1. Primary: the median gain in C-index for overall survival, per family, and the count of signature × cohort
   tests clearing the null's 95th percentile, against the 5% expected.
2. Primary contrast: inflamed/interferon against confounder controls, on the z scale, tested by a bootstrap
   over the 13 signature correlation clusters (one-sided) and by the within-cohort permutation used elsewhere.
3. Comparison with the two existing settings: the inflamed family's z under checkpoint blockade against its z
   for survival without checkpoint blockade (+0.27) and for response (+1.54).

## What would change the paper

The claim that immune signatures carry treatment-specific information is strengthened only if the inflamed
family's gain under checkpoint blockade exceeds the confounder controls' at a one-sided cluster-bootstrap
P < 0.05. If it does not, the result is reported as it comes out, beside the existing benchmarks, and the
paper's claims are left as they are.
