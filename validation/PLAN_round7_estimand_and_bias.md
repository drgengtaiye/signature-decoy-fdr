# Analysis plan, round 7: the estimand made operational, the direction of the bias, and decoy-draw variability

Written 2026-10-22, after an external review of manuscript v6 and before any analysis below was run. Its SHA-256 is
recorded in PLAN_round7_estimand_and_bias.sha256 and the file is deposited before the first analysis starts.

## What is known when this is written

* Every result of manuscript v6: the response benchmark (21 folded / 27 directed calls; FDR 17% and 7% with the
  first decoy set, 19.9% and 8.0% with the package's decoys, seed 20261022); the structural-truth reading of the
  round-5 planted signal (estimate at or above the realized proportion in all 16 settings, truth = expected delta
  above the 95th percentile of the decoys'); TCGA overall survival 13 folded / 16 directed calls, FDR 52% (target
  bootstrap 34-98%, cluster bootstrap 22-123%) and 39%, DSS 56%/56%, PFI 23%/34%, global-null ratio 1.01 and 0.89;
  breast cancer (45 Venet sets, 13 cohorts) first principal component 31 calls, FDR 13%, mean z 23 calls, 44%.
* The review's points: the structural truth was defined after the gene-content result and does not match the
  formal null (call probability at most that of the decoys); the plain estimate counts every target's decoy rate
  and so is an upper bound when some targets are true (pi0 = 1); signal outside the library makes the estimate
  high while under-matched decoy coherence makes it low, so the net bias in TCGA is unknown; breast cohorts include
  some signatures' development cohorts; one decoy draw moves the headline (17% vs 20%); the target-bootstrap
  interval covered the realized proportion in only 10-45% of planted replicates (gene-content truth); the decoy
  pool shrinks as a library grows.
* No analysis below has been run. The provenance of the breast signatures is being compiled from the original
  publications by a separate reader who sees no results (analysis 5).

Seeds below are arbitrary labels.

## 0. Package changes (decoyfdr 0.2.0), before any analysis that uses the package

* grow_decoys() gains n_sets: n_sets independent sets of n_coex co-expressed and n_rand random-gene decoys per
  target. Set 1 is drawn from the target's own L'Ecuyer-CMRG stream exactly as in version 0.1.0, set s > 1 from
  substream s - 1 of that stream, so set 1 reproduces the decoys of a 0.1.0 run with the same seed. Decoys are
  ordered set by set, so the weight permutations fit_units() draws for weighted-mean decoys in set 1 are those of
  0.1.0 too. Each decoy carries its set.
* decoy_fdr() reports, beside the estimate pooled over all sets (equal to the mean of the per-set estimates), the
  smallest and largest per-set estimate, and the pi0-corrected estimate of analysis 3 with its interval.
* run_decoy_fdr() defaults to n_sets = 10.
* Check: rerunning the round-6 response, TCGA and breast scripts with n_sets = 10 reproduces, in set 1, the round-6
  decoys, calls and estimates exactly (analysis 4). R CMD check with no warning; tests and vignette pass in a fresh
  R library with only CRAN and Bioconductor installs.

## 1. Structural truth by call probability (response, planted signal of round 5)

The formal null of a target s is P(s called) <= pi_s, pi_s being the call probability of its eligible decoys. It
is made operational within each replicate of the 16 round-5 planted settings (8 scenarios x folded/directed; 200
replicates), each replicate's labels, calls and estimate regenerated exactly as in v10_calib_sim.py signal and
checked against v10_calib_signal.csv. For that replicate's causal set and beta, M = 200 further label sets are
drawn by the same generative model from a separate stream (seed 710000 + 1000 x scenario index + replicate); the
full call procedure (with the replicate's cohorts) is run on each. f_t = target t's call frequency over the M
sets; pi_t = mean call frequency over the M sets of its eligible decoys (pool-A co-expressed; in the target's
declared direction for the directed statistic; in P3 sharing no causal gene, as the estimator uses them).
Primary truth: t is null when f_t <= pi_t. Sensitivity: t is true only when a one-sided exact binomial test of
f_t against pi_t (M trials) gives P < 0.05.

Reported per setting, as ratios of totals over replicates: realized FDR = sum V / sum R (V = calls of null
targets in the replicate's own labels), the plain estimate sum E / sum R, the pi0-corrected estimate (analysis 3);
and the decomposition of the plain estimate's expectation into sum over null targets of pi_t, sum over true
targets of pi_t (the part that is conservative by construction) and sum over null targets of (pi_t - f_t) (the
decoys' excess over null targets). The round-6 truth (expected delta above the decoys' 95th percentile) is
reported beside it; this truth replaces it as the primary reading.

## 2. Prospective planted signal in TCGA overall survival, with new causal gene sets

Data: icb/v11_pkg/tcga_OS_fit.rds (20 types; 52 targets; 32 co-expressed and 8 random decoys per target, set 1;
200 kept null draws per target and type; kept scores). Folded statistic only (analysis 7). Six settings, 100
replicates each:

* T1 library: causal genes = the union of the genes of 2 targets drawn uniformly from the 52.
* T2 outside: causal genes = the genes of 1 MSigDB Hallmark set (msigdbr, the version installed) drawn uniformly
  among those with at least 15 genes in the decoy pool, restricted to its pool genes.
* T3 mixed: the union of the genes of 1 target drawn uniformly and of a T2 set.
* beta = 0.3 and 0.6 (log hazard ratio per standard deviation of the causal score).

Outcome model, in each type: u = mean z over the causal genes measured (at least 3; otherwise beta = 0 in that
type), standardized; event times exponential with hazard lambda exp(-beta u) / mean(exp(-beta u)), lambda = events
/ total follow-up time of the type; censoring times drawn with replacement from the reverse Kaplan-Meier
distribution of the type's follow-up; observed time and event = min and indicator. A higher causal score means
longer survival. The replicate's outcome and the M = 100 truth draws come from streams seeded 720000 + 1000 x
setting + replicate and 730000 + 1000 x setting + replicate. Every stored row (targets, decoys, kept null draws) is
read by Harrell's C, deltas against the kept bank's folded median, signed-rank calls over the 20 types (at least
6). Decoys sharing causal genes are kept, as they would be in real data. Reported as in analysis 1 (realized FDR
under the call-probability truth, plain and pi0-corrected estimates, decomposition, coverage of analysis 6).

## 3. pi0-corrected estimate

For each target t with eligible decoys d (as the estimate uses them), s = sign(median delta) x (-log10 P) of the
signed-rank test; q_t = (1 + #{d: s_d >= s_t}) / (1 + n_t). pi0 = min(1, #{t: q_t > lambda} / (m (1 - lambda))),
lambda = 0.5 (sensitivity 0.25 and 0.75). Corrected FDR = pi0 x E / R. Intervals: bootstrap over targets and over
clusters, pi0 recomputed in every resample. Reported for every application beside the plain estimate, which is
described as the estimate with pi0 = 1. Its calibration is read in analyses 1 and 2 and in the 1,000 label
permutations of round 5 (v10_calib_sim.py null).

## 4. Ten independent decoy sets

The headline estimate of every application becomes the estimate pooled over 10 decoy sets (32 co-expressed and 8
random-gene decoys each), the package default, with the range of the 10 per-set estimates. Reruns with the seeds of
round 6 (decoys 20261022, scoring 20261023, global null 20261024) and decoyfdr 0.2.0: response (folded and
directed), TCGA overall survival with its global null, DSS and PFI, breast first principal component with its
global null, breast mean z, and the breast absolute-coherence decoys (seed of round 6). The closest-matched
variants as in round 6. The first decoy set of the response analysis (calib_decoys_v3.R, 17%), on which the
comparisons between nulls (Table 2) and the calibration experiments were planned, stays the basis of those
analyses and is reported as one draw beside the ten.

## 5. Breast cancer: origin cohorts left out

A provenance matrix of the 45 measured Venet sets x 13 cohorts, coded from each original publication by the rule of
signature_provenance_notes.md (D: outcome used to pick genes or fit weights, cut-offs or a model; C: the cohort's
own outcome prompted or shaped the signature; V: a finished signature only tested; blank: not used), through the
patient overlaps between series (Uppsala patients in GSE4922 and Loi; GSE2990 within Loi; NKI-295 containing 61 of
van 't Veer's 78). Frozen (SHA-256 recorded) before it is applied. D and C tests are left out for the target and
its decoys (call_units(exclude =)), and the analysis rerun on the 10-set fit; a conservative variant also leaves
out the borderline cases the reader lists. The left-out estimate becomes the breast headline; the full estimate is
reported beside it.

## 6. Intervals

In analyses 1 and 2, every replicate gets a 95% interval from 1,000 bootstrap resamples over targets and over the
13 correlation clusters (v5_signature_clusters.csv, cluster_fixed), for the plain and the pi0-corrected estimate.
Coverage = share of replicates (R >= 1) whose interval contains (a) the setting's realized FDR (sum V / sum R) and
(b) the replicate's own realized proportion. From this version the cluster-bootstrap interval is the one given in
the text and tables, whatever the coverage; the target bootstrap is given in Additional file 2. Breast clusters:
average linkage on 1 - |median over cohorts of the Spearman correlation between first-principal-component scores|,
cut at 0.5. Hallmark clusters (analysis 8) the same on mean-z scores over the 20 TCGA types.

## 7. TCGA: folded only

The declared direction of an immune signature concerns response to checkpoint blockade; immune infiltration is
prognostic in opposite directions in different cancer types. The directed statistic therefore does not test a
declared claim in TCGA, and the TCGA results in the text are folded only; the directed numbers move to Additional
file 2 with this reason.

## 8. The decoy pool at library scale

(a) Pool size (genes measured in all 20 TCGA types and in no set) for: the 52 immune signatures with the curated
exclusion set (as used), MSigDB Hallmark, Hallmark + C2 canonical pathways, Hallmark + C2 CP + C5 GO biological
process. (b) For each of the 50 Hallmark sets as a target, 4 co-expressed decoys (mean engine, tol 0.02, seed
740000) grown across the 20 types from each of the three MSigDB pools; reported: median ratio of achieved to
target internal correlation and the share of targets with all 4 decoys within 0.05 of the target. (c) The 50
Hallmark sets, scored by mean z, on TCGA overall survival with the package (n_sets = 10, 1,000 null draws, pool =
genes in no Hallmark set, seeds 741000/741001/741002): folded calls, plain and pi0-corrected FDR with cluster
intervals, global null with 200 permutations. A demonstration; no decision rests on it.

## 9. Text

The TCGA estimate is reported as the pi0 = 1 estimate with its cluster interval and the pi0-corrected estimate;
"about half ... are false" is not used. Both opposing biases (under-matched decoy coherence; signal outside the
library) are stated with the result of analysis 2. The Discussion limits "any library" by analysis 8.

## Decision rules

* All numbers are reported whatever they are. Numbers of v6 that change because of analysis 4 are given beside the
  new ones in Additional file 2.
* If set 1 of analysis 4 does not reproduce round 6, the package is corrected before anything else is reported.
* If the pi0-corrected estimate falls below the realized FDR in any setting of analyses 1-2 by more than a third
  of the realized value, the text says so wherever the corrected estimate is given and the plain estimate stays
  the headline; otherwise both are given with equal weight.
