# Analysis plan, round 14: the response calls in checkpoint cohorts outside the collection

Written 2026-11-02, after manuscript v15 and before any sample annotation or expression value of the candidate series
below was read. Its SHA-256 is recorded in PLAN_round14_new_checkpoint_cohorts.sha256 and the file is deposited before
the first sample annotation is read.

## Question

The response analysis calls 21 of 52 signatures (folded AUC; estimated FDR 18%, cluster interval 11-240%, pi0-corrected
4%) in 12 checkpoint cohorts of the ORCESTRA collection. GSE176307 (round 11) is the only checkpoint cohort outside the
collection that has been examined. This plan asks, in further public checkpoint cohorts that played no part in any of
it, whether the called signatures beat their co-expressed decoys, whether they do so more than the uncalled ones, and
what share of the calls fail to replicate beyond what a false call would show.

## What is known when this is written

* Every result of manuscript v15: the calls, decoys and estimates of the 12 cohorts (icb/v12_pkg/response_*), the
  family means of rho for response (inflamed 0.74, other 0.58, confounder 0.62; icb/v15_content/family_by_setting.csv)
  and GSE176307 (inflamed rho 0.91; icb/v17_external/predictions.csv).
* The GEO search below (titles and summaries only; validation/round14_search/): run 2026-11-01 through E-utilities
  (esearch, db = gds), query in round14_search/query.txt, 2,964 series; 303 human series with at least 25 samples
  whose title or summary names a checkpoint agent or checkpoint blockade and a response term (shortlist.csv). No
  sample annotation, supplementary file or expression value of any of them has been read.

## A. Cohorts

A1. Candidates: the 303 shortlisted series, and series the authors know from the literature (GSE274975, GSE177043,
GSE111414; not a systematic addition, stated as such). Read for each, in this order: GEO series record, then sample
annotations (GEO query pages, brief view: no expression values), then, only where response is not in the annotations,
the publication's open full text and tables in PubMed Central.

A2. Eligible when all hold: (1) bulk tumor tissue sampled before checkpoint blockade (one sample per patient, the
earliest before treatment); (2) treatment with anti-PD-1, anti-PD-L1 or anti-CTLA-4 antibodies, alone or combined with
one another, for advanced or metastatic disease (combinations with chemotherapy, targeted agents, radiotherapy,
vaccines, oncolytic viruses or other investigational agents, and neoadjuvant or perioperative settings, are not
eligible); (3) response by RECIST 1.1, iRECIST, irRC, mRECIST or RANO as best overall response, complete or partial
response against stable or progressive disease, or the depositors' responder labels derived from it (labels that count
stable disease as response are used only where nothing else is given, and flagged); (4) expression on a platform
measuring at least 3,000 genes after mapping to symbols, obtainable through the GEO query pages or the GEO download
endpoint (the FTP server is unreachable from this cluster); (5) at least 25 evaluable patients and 5 per class;
(6) no patient in common with the 21 ORCESTRA cohorts or GSE176307. The eligible list, with the reason for each
ineligible candidate that passed (1) and (2), is frozen in an amendment deposited before any expression value is
downloaded. Controlled-access trial data (for example JAVELIN Renal 101, IMmotion150/151, POPLAR/OAK, IMbrave150,
SU2C-MARK, Hartwig) cannot be obtained from this cluster without data-access approval; the plan applies to them
unchanged, by amendment, if approval is granted.

## B. Procedure

As round 11's response application, with the eligible cohorts jointly: the 52 signatures by their own engines (gene
symbols; a set with fewer than three genes measured in some cohort left out), decoys grown by decoyfdr (frozen install
icb/v12_pkg/Rlib6) from the genes measured in every eligible cohort and in no curated signature, ten sets of 32
co-expressed and 8 random-gene decoys (seed 20261121), 1,000-draw exact-size banks per target and cohort (seed
20261122), folded and directed AUC. rho as in round 10 (decoys/v15_content_structure.py). Expression: TPM or
normalized values as deposited, log2(x + 1) where on a linear scale; for raw counts, counts per million then
log2(x + 1). Script applications/v20_new_checkpoint.R; analysis decoys/v20_new_checkpoint.py; outputs icb/v20_new_icb/.

## C. Predictions and tests

Over the eligible cohorts (each cohort weighted equally; statistic folded unless stated):

1. Primary. The 21 called signatures beat their co-expressed decoys: mean rho over the 21 x cohorts above 0.5.
   Reported as holding when the point estimate exceeds 0.5 and as confirmed when its 95% interval over the 13
   correlation clusters (10,000 resamples, seed 916201) excludes 0.5.
2. Called above uncalled: mean rho of the 21 minus that of the other signatures, above 0, with the same interval.
3. Replication. r, the share of the 21 whose median folded Delta AUC (against the bank median) over the cohorts is
   above 0; d, the same share among their co-expressed decoys (ten sets pooled); pi_hat = (1 - r) / (1 - d), bounded
   to [0, 1]. pi_hat estimates the share of false calls from replication alone and is biased upward, because true
   calls can fail to replicate. Prediction: pi_hat below 0.5; reported beside the estimate of 18% (4% corrected),
   with a cluster interval.
4. The inflamed family's mean rho above 0.5 (12 cohorts 0.74), with the same interval.
5. Descriptive: if at least five cohorts are eligible, the calls and decoy estimate of the new cohorts on their own
   (the package's call_units and decoy_fdr, as in round 7), and their overlap with the 21; per-cohort family rho; the
   directed statistic for 1 to 4; the same with GSE176307 added.

If fewer than two cohorts are eligible the tests are not run, and the screen is reported with its reasons. Every result
is reported whatever it shows.

## D. Displays and text

Panels in a main or supplementary figure, a Results paragraph, a Methods paragraph, Additional file 2 sheets.
Numbers already reported do not change.
