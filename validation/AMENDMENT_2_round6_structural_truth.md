# Amendment to PLAN_round6_applications.md, analysis 4

Written 2026-10-21, before analysis 4 was run.

The plan estimates each signature's expected directed delta as its mean over the 200 replicates of a
scenario. In scenarios P1 and P2 of round 5 the causal gene set is redrawn in every replicate (and in P3
the causal decoy), so a mean over replicates averages over different truths and defines none of them.

The expected delta is instead estimated within each replicate: for that replicate's causal set and beta,
M = 20 further label sets are drawn by the same generative model (responders drawn without replacement
with weights exp(beta u), their number kept), from a separate random stream (seed 20261025 + 1000 x
scenario index + replicate), so that the replicate's own labels, calls and estimate are exactly those of
round 5. For every target and every pool-A co-expressed decoy the delta is averaged over the M draws in each
cohort, and the median over the replicate's cohorts is its expected delta. A target is structurally true in
that replicate when its expected delta exceeds the 95th percentile of the expected deltas of its eligible
decoys (those the estimator uses: in the target's declared direction for the directed statistic; in P3,
those sharing no causal gene). The realized false discovery proportion under this truth is reported, as a
ratio of totals over replicates, beside the decoy estimate and the gene-content proportion of round 5.
The folded statistic is reported beside the directed one. Round 5's calls are to be reproduced exactly
before the new quantity is computed.
