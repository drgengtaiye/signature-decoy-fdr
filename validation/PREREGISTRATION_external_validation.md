# External validation of the locked triage predictions — registration before scoring

Written 2026-09-18, before any validation cohort below has been downloaded or scored. Nothing
in this file may change once its SHA-256 has been recorded and the file deposited (private S3
object; the object's LastModified time is the external timestamp).

## 1. What is being validated

`icb/v3_LOCKED_predictions.csv`, created 2026-09-14, SHA-256
`aea506662efef7149cf336004bb1c0c76349206314fc537f950bdac73aa0d9ca` (verified again on
2026-09-18, before this file was written), with model and endpoints in
`icb/v3_LOCKED_predictions.json`. It gives, for 50 signatures, a `predicted_effect`, an
`immune_attributed` flag and a `triage_call`. None of these values is refitted, re-thresholded
or otherwise changed for this validation.

## 2. How the cohorts were found

A systematic GEO search on 2026-09-18 (NCBI E-utilities, db = gds), with the query

    (pembrolizumab OR nivolumab OR atezolizumab OR durvalumab OR "anti-PD-1" OR "anti-PD-L1"
     OR ipilimumab OR "immune checkpoint" OR "checkpoint blockade" OR "checkpoint inhibitor")
    AND "Homo sapiens"[Organism] AND gse[Entry Type]
    AND ("expression profiling by high throughput sequencing" OR "expression profiling by array")

returned 1,006 series (399 with at least 25 samples). Titles and summaries were screened for
patient tumor bulk expression (excluding single-cell, blood, sorted cells, cell lines,
organoids and animal models), and the sample annotations of every remaining candidate, and of
known checkpoint cohorts outside the search terms, were read. No expression value of any
candidate had been downloaded when the lists below were fixed. None of these cohorts was used
at any stage of the development of the signatures' benchmark, the attribution or the model;
the 12 development cohorts all come from the ORCESTRA clinical_icb collection, and the 8
further ORCESTRA cohorts that failed its inclusion rule remain excluded.

## 3. Eligibility (the development rule, unchanged)

Checkpoint blockade (anti-PD-1, anti-PD-L1 or anti-CTLA-4, alone or combined, without
concurrent chemotherapy); a pre-treatment tumor sample with bulk expression; RECIST best
response, responders CR or PR against non-responders SD or PD; at least 25 evaluable patients
with at least 5 in each class; one sample per patient (the pre-treatment one). As in
development, a cohort whose matrix retains fewer than 3,000 genes after filtering is not a
whole-transcriptome cohort and cannot enter the primary analysis.

## 4. The cohorts, fixed now

**Primary set** (whole transcriptome, eligible as above):

| GEO | tumor | therapy | evaluable (CR/PR : SD/PD), from the annotations |
|---|---|---|---|
| GSE161537 | non-small-cell lung cancer | anti-PD-1/PD-L1, second line or later | 76 (20 : 56) |
| GSE159067 | head and neck squamous cell carcinoma | anti-PD-1/PD-L1 | 102 (11 : 91) |
| GSE303300 | cervical squamous cell carcinoma | checkpoint window before chemoradiation; response assessed before chemoradiation | 36 baseline (5 : 31) |
| GSE176307 | metastatic urothelial carcinoma | anti-PD-1/PD-L1 | 89 (16 : 73); enters the primary set only if its matrix passes the 3,000-gene rule, otherwise the panel set |

**Panel set** (secondary; targeted expression panels, below the 3,000-gene rule): GSE215868
(melanoma, NanoString IO 360; 104 evaluable, CR/PR 45 : SD/PD 59) and GSE93157 (lung, head
and neck and melanoma, NanoString PanCancer Immune Profiling; 65, CR/PR 20 : SD/PD 45). Here
the matched null is drawn from the panel's own genes, and signatures keep only the genes the
panel carries (at least 3, as in development).

**Chemo-immunotherapy set** (exploratory only; response is not attributable to the checkpoint
antibody): GSE207422 (neoadjuvant anti-PD-1 plus chemotherapy, lung; major pathological
response against none) and GSE235919 (durvalumab and tremelimumab with chemotherapy,
colorectal; CR/PR against SD).

**Randomized** (secondary): I-SPY2, GSE194040 — paclitaxel plus pembrolizumab (n = 69) against
the HER2-negative paclitaxel controls (n = 179), pathological complete response; and
GSE173839 — durvalumab plus olaparib with paclitaxel (n = 71) against controls (n = 34),
evaluable pCR only.

## 5. Measurement (the manuscript's Methods, unchanged)

Each signature is scored by its own publication's method on the log-scale matrix (GSVA,
ssGSEA or weighted mean of gene z-scores), against a null of 200 random gene sets per rung of
the size ladder (4, 7, 12, 20, 32, 50, 80, 130, 260, 440 genes) per engine, drawn from the
cohort's own genes, seed = sum of the character codes of the cohort name. The measured effect
of a signature in a cohort is its folded AUC minus the median AUC of its matched null.

## 6. Endpoints, fixed now

*Primary (from the locked file).* In each primary cohort, Spearman ρ between
`predicted_effect` and the measured effect over the locked signatures scorable in that cohort;
success in a cohort is ρ > 0 at one-sided P < 0.05. The overall primary test combines the
primary cohorts' one-sided P values by Stouffer's method (unweighted); overall success is a
combined one-sided P < 0.05. Every cohort is reported whatever it shows.

*Sensitivity to signature dependence.* The mean ρ across primary cohorts is referred to a
permutation of `predicted_effect` within the 12 fixed correlation clusters of the 50
signatures (`icb/v6_cluster_restricted.csv`, `cluster_fixed`; 2,000 permutations, one
permutation applied to all cohorts at once).

*Secondary.*
1. (From the locked file.) Difference in measured effect between `immune_attributed` 1 and 0
   signatures, per cohort and pooled over primary cohorts (difference in medians, bootstrap 95%
   interval over signatures).
2. The triage call: among signatures called 'test for treatment response' against 'not worth a
   response cohort', the share whose median measured effect across the primary cohorts exceeds
   0.05 (Fisher's exact test, one-sided).
3. The family contrast: share of signature × cohort tests clearing the null's 95th percentile,
   inflamed/IFN against confounder controls, pooled over primary cohorts.
4. The panel set and the chemo-immunotherapy set analysed as the primary set, reported
   separately and not combined with it.
5. Randomized: for each signature, logistic regression pCR ~ score + arm + score × arm (score
   standardized within the trial); the number of the 17 immune gene sets of the trial library,
   and of the 50 locked signatures, with a checkpoint-favoring interaction at P < 0.05, referred
   to the no-interaction simulation (pCR simulated from the fitted pCR ~ score + arm model for
   each signature jointly through the first K principal components of the score matrix, K = 9,
   and the whole pipeline rerun; 2,000 replicates), as for the survival trials.

## 7. What failure looks like

A combined primary P at or above 0.05, or a negative pooled ρ, means the triage rule has not
been shown to transfer to cohorts it has never seen, and the manuscript will say so in the
abstract and restrict the algorithm's claims accordingly.
