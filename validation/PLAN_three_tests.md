# Analysis plan: three further tests

Written before any of the three quantities below was computed, and deposited to a private cloud object
store whose object timestamp is the external record. Three hypotheses are stated, each with one primary
endpoint and its own decision rule. They are independent questions and no correction is applied across
them; each is reported as one of three, whatever the outcome.

## Test A. What makes a cohort able to show the effect

Signature effects differ more between cohorts than between signatures. The hypothesis is mechanical: a
signature can only separate responders where the thing it measures actually varies. If every patient in
a cohort is equally inflamed, no immune signature can discriminate there.

- **Predictor, computed without outcome.** In each cohort, every gene is standardised across patients,
  the inflamed/interferon family's genes are averaged per patient, and the predictor is the standard
  deviation of that score across the cohort's patients.
- **Outcome.** The cohort's mean effect over the 52 signatures, each measured as its AUC minus the
  median of its own size- and engine-matched null, as already computed for the benchmark.
- **Cohorts.** The 12 development cohorts and the 7 registered validation cohorts, 19 in all.
- **Primary.** One-sided Spearman correlation between predictor and outcome, positive direction.
  Significant at P < 0.05.
- **Secondary.** The same with the renal cohorts removed; the partial correlation given log cohort size;
  and the predictor against the inflamed family's own mean effect.

## Test B. Whether the library is one axis

Fifty-two signatures are published as if they were alternatives. If they are one axis, a score that
takes no position on which of them is right should do as well as the best of them.

- **The axis, computed without outcome.** In each held-out cohort the 52 signature scores are
  standardised across patients and their first principal component is taken, its sign fixed so that it
  correlates positively with the inflamed family's mean score.
- **The comparator, fixed in advance.** The single signature with the largest mean effect across the 12
  development cohorts, chosen before any validation cohort is scored.
- **Primary.** In the 7 registered validation cohorts, the paired difference in folded AUC between the
  axis and that signature, with a cluster bootstrap over cohorts (4,000 resamples). The axis is declared
  non-inferior if the lower bound of the one-sided 95% interval exceeds -0.02 AUC.
- **Secondary.** The axis's rank among the 52 in each cohort, and its effect against the matched null.

## Test C. The same treatment specificity, one gene at a time

The paper shows that signatures clear their matched null inside checkpoint arms and not inside the
comparator arms of the same trials. If that is a property of what immune genes measure and not of how
signatures are assembled, it should hold gene by gene.

- **Genes.** The union of the genes of the 17 inflamed/interferon signatures present in a trial, against
  the union of the confounder-control signatures' genes as the control set.
- **Within each arm.** IMmotion150 and PRINCE: the Cox coefficient per standard deviation for
  progression-free survival. I-SPY2's two contrasts: the AUC for pathological complete response. Each
  gene is compared against 200 random genes of the same trial in the same patients, one-sided in the
  benefit direction, exactly as the signature-level within-arm analysis was done.
- **Arm sizes.** The larger arm is subsampled to the smaller arm's size, 40 draws, and the matched
  counts are the primary ones.
- **Primary.** Across the four randomised comparisons, the difference between checkpoint and comparator
  arms in the fraction of immune genes clearing their within-arm null, against 2,000 permutations of the
  arm label within trial. One-sided. Significant at P < 0.05.
- **Secondary.** The same for the control gene set, where no difference is expected.

## What none of these can show

All three reuse cohorts already analysed here; none adds a patient. A, if positive, describes which
cohorts can show an effect, not why a patient responds. B, if positive, says the published library is
redundant, not that the axis is a useful test. C shares its patients with the signature-level result and
is therefore a decomposition of it, not an independent confirmation.
