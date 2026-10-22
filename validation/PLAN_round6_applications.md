# Analysis plan, round 6: software, two further applications, and the estimand

Written 2026-10-21, after an external review of manuscript v5 and before any analysis below was run.
Its SHA-256 is recorded in PLAN_round6_applications.sha256 and the file is deposited before the first
analysis starts.

## What is known when this is written

* Every result in manuscript v5, including the prognostic benchmark of the 52 signatures in 20 TCGA
  cancer types against a ladder null (Fig. 7g–j of v5: median gain +0.009 in C-index; six signatures
  exceeding their null more often than chance, four of them confounder controls; the confounder family
  clearing its null in 18.1% of cancer types against 10.5% for the inflamed family). The prognostic
  application below is therefore not blind to how the targets behave; it is new only in its decoys.
* The published findings of Venet et al. (PLoS Comput Biol 2011) in breast cancer: most random gene sets
  are associated with outcome, through proliferation, and many published signatures are not better than
  random sets of the same size. Their supplementary bundle (Dataset S1) has been downloaded and its
  signature object parsed without R: 48 sets (47 published signatures and the meta-PCNA metagene) with
  Entrez identifiers, symbols and PMIDs. No expression or outcome in it has been read.
* The software review asked for (an installable tool) does not exist yet.

## 1. Software

An R package, decoyfdr, in the repository subdirectory decoyfdr/, implementing the procedure of the paper:
cohort preparation (finite, non-constant genes; gene-wise z-scores; expression deciles), scoring by mean
z, signed weighted mean, GSVA, ssGSEA or first principal component; the signed statistic (AUC for a
binary outcome; Harrell's C for survival, larger score with longer survival); exact-size matched null
banks (signed for weights of both signs); co-expression-matched decoys grown by the multi-cohort
seed-and-grow of calib_decoys_v3.R (size, expression-decile composition and internal correlation in every
cohort; near-best band 0.02), plus uniform random-gene decoys; the folded and directed deltas; the
Wilcoxon signed-rank call over cohorts (P < 0.05 with a positive median, at least 6 cohorts), with the
exact tie-aware P of v10_calib_lib.signed_rank_p; the decoy FDR with a bootstrap over targets or over
clusters; and the global-null calibration by within-cohort outcome permutation.

Validation, required before the package is used below:
* V1. From the per-cohort statistics stored in calib_v3/ (N0), the package reproduces the 21 folded and 27
  directed calls and the FDRs of 17.3% and 7.4% exactly (point estimates; intervals to Monte Carlo error).
* V2. Its AUC equals calib_score_v3.R's auc_s and its C equals survival::concordance on test data, to 1e-12.
* V3. Rerun end to end on the 12 response cohorts with decoys and null banks drawn afresh by the package
  (seed 20261022): calls and FDRs are reported beside the published ones whatever they are. No decision
  rests on V3; it measures how much the estimate moves when every random element is redrawn.

## 2. Application A: prognosis in 20 TCGA cancer types

Cohorts, samples and endpoint exactly as in null_matched_prog.R: primary tumors, overall survival
(primary), disease-specific survival and the progression-free interval (secondary), the 20 types with at
least 80 patients and 25 events. Library: the 52 signatures, each scored by its publication's engine.
Null: N0, 1,000 draws of the signature's own size per cancer type from all measured genes (signed for the
four weighted-mean signatures with weights of both signs). Decoys: grown by the package across the 20
types, 32 co-expressed from pool A (genes in no curated signature) and 8 random-gene decoys per target,
seed 20261022. Statistics: folded C (primary, the benchmark's statistic) and directed C (secondary; a
signature declared to mark sensitivity is taken to mark longer survival, one declared to mark resistance
shorter survival). Unit of analysis: Wilcoxon signed-rank over the 20 types on delta = C - median of the
null bank under the same transformation; a call when P < 0.05 with a positive median. FDR = mean
per-target decoy call rate x 52 / signatures called, with bootstraps over targets and over the 13
correlation clusters of Fig. 2d (6,000 replicates). Global null: (time, event) pairs permuted within each
type, 200 times (seed 20261022), the whole procedure rerun on 200 of the 1,000 null draws per target;
rho = mean estimated / mean realized false calls with a 95% interval over permutations. The decoys'
internal correlation against their targets' is reported per type.

## 3. Application B: a non-immune library, breast cancer prognosis

Library: the 48 sets of the Venet bundle's cancer.signatures object; sets with fewer than 3 genes
measured in every cohort are excluded. Scoring: the first principal component of the gene-wise z-scores
(Venet's own score; its sign is fixed to correlate positively with the mean z and the statistic is
folded), and the mean z as a sensitivity analysis. Only the folded statistic is used: the signatures carry
no declared direction.

Cohorts: NKI and Loi (GSE6532) from the Venet bundle; TCGA-BRCA (primary tumors, as in Application A);
and every one of the following public series that can be obtained with expression and outcome:
GSE2034, GSE7390, GSE11121, GSE2990, GSE1456, GSE3494, GSE4922, GSE20685, GSE21653, GSE42568, GSE17705,
GSE12093, GSE96058 and METABRIC. Inclusion: at least 100 patients and 25 events for the endpoint and at
least 10,000 genes after mapping probes to Entrez genes (the probe with the largest mean per gene); a
series duplicated in another (same patients) is kept once, the larger. Endpoint, one per cohort, the first
available of distant-metastasis-free survival, relapse- or disease-free survival, disease-specific
survival and overall survival. Expression is used as deposited (log2 if not already).

Null: N0, 1,000 draws of the set's own size per cohort. Decoys: 32 co-expressed and 8 random-gene decoys
per signature from genes in none of the 48 sets, measured in every cohort, grown across the cohorts (seed
20261022). Unit, call, FDR and global null (200 permutations) as in Application A. If fewer than 6 cohorts
qualify, the unit of analysis is instead the per-cohort test (C above the 95th percentile of the null) in
each cohort and the FDR is reported per cohort.

## 4. The estimand, and the planted signal read against it

The estimate targets the expected number of signatures called among those whose statistic is
exchangeable with that of structure-matched modules of non-signature genes (H0_s). Under the global null
this holds for every signature, which analysis 3 of round 5 tested. In the planted-signal scenarios of
round 5 the truth was defined by gene content. Here the same replicates (same seeds, same causal sets)
are re-read with a structural truth: for each scenario, a signature's expected directed delta is
estimated as its mean over the 200 replicates (per cohort, then the median over cohorts), and likewise
for each of its pool-A co-expressed decoys; a signature is structurally true when its expected delta
exceeds the 95th percentile of its decoys' expected deltas. The realized false discovery proportion under
this truth is set beside the decoy estimate and the gene-content FDP. This truth definition was chosen
after the round-5 planted results were known and is reported as such. No decision rests on it.

## 5. Scope of the manuscript

The randomized-trial analyses (seven comparisons, the within-arm nulls, I-SPY2 and the three further
trials), the renal contrast, survival under checkpoint blockade, the gastric survival analyses, the
information ceiling and the cross-setting dissociation are removed from this manuscript. They remain in
manuscript v5, Additional file 2 of v5 and the repository. None of them is recomputed for the removal,
and none is cited as a result of this paper.

## Decision rules

* The FDR, the calls and the calibration ratio of Applications A and B are reported whatever they are.
  If rho lies outside 0.8–1.25 in an application, the direction is stated wherever that application's
  FDR is reported and the estimate is not described as calibrated there.
* If V1 fails, the package is corrected until it passes before anything else is run with it.
* The paper's title and abstract describe the estimate as transferring to other libraries and endpoints
  only to the extent that Applications A and B are calibrated under their global null.
