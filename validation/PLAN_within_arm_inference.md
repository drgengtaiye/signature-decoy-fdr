# Analysis plan: inference for the within-arm contrast

Written before the test was run and deposited to a private cloud object store, whose object timestamp
is the external record.

## The gap this fills

The paper's central claim is that what immune signatures measure is specific to the treatment given,
and its cleanest evidence is the within-arm comparison: immune signatures clear a matched null inside
the checkpoint arms of randomized trials and not inside the comparator arms of the same trials. That
comparison is currently reported as descriptive. The statistic is already defined and already
published in the manuscript; what is missing is a test of it against a null in which the arm label
carries no information. This plan adds that test and nothing else. No statistic, cohort, family or
threshold is changed.

## The statistic, as already defined in the paper

For one randomized comparison, the difference between arms in the fraction of the 17 immune signatures
that clear their own within-arm null in the benefit direction. Each signature is scored as elsewhere in
this paper, fitted inside one arm (a Cox coefficient per standard deviation for progression-free
survival in IMmotion150 and PRINCE, a directional AUC for pathological complete response in the two
I-SPY2 contrasts), and compared with 200 random gene sets of the size rung nearest its own, drawn from
that trial's own genes and scored in the same patients. Clearing means exceeding the 95th percentile of
that null in the benefit direction. The observed values are those already in the manuscript: 9 of 17
against 0 of 17 in IMmotion150, 6 of 17 against 0 of 17 in PRINCE, 14 against 15 in the I-SPY2
pembrolizumab contrast and 6 against 4 in the durvalumab contrast.

## Test

The arm label is permuted within each trial, 2,000 times. Everything else is recomputed inside the
permuted arms, including the random-set null, so a permutation is a complete repetition of the
procedure on patients who have been relabelled. A permutation preserves each arm's size, so the
imbalance between arms and the correlation between signatures both enter the null.

- **Primary.** The four one-sided permutation P values combined by Fisher's method. Significant at
  P < 0.05.
- **Secondary, stated in advance.** Each trial's own P value; the same test for the three confounder
  controls, where no difference is expected; and the two trials with survival endpoints separately,
  since in the I-SPY2 contrasts both arms receive chemotherapy and the manuscript already reports that
  the inflamed state predicts response to chemotherapy there.

**What counts as a failure.** If the combination does not reach P < 0.05 the within-arm contrast stays
descriptive in the manuscript and is reported as untestable at this sample size. No other statistic
will be substituted afterwards.

## What this cannot show

It tests the arm contrast within these four comparisons only. It does not make the comparison
randomized with respect to the signature, it adds no patient, and a positive result does not establish
that any individual signature predicts benefit.
