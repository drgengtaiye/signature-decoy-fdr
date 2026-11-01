# Analysis plan, round 11: gene content and structure in independent cohorts

Written 2026-11-01, after manuscript v13 and before any expression value of the survival cohorts below was downloaded.
Its SHA-256 is recorded in PLAN_round11_external_replication.sha256 and the file is deposited before the first download.

## Question

Round 10 (PLAN_round10_content_and_structure.md) found, in TCGA and the 12 response cohorts, that the inflamed
signatures' genes add to what modules of their structure achieve for response to checkpoint blockade but, in most
cancer types, not for survival, and that where they do so for survival (melanoma, sarcoma, cervical cancer) and where
they do not (stomach, lung squamous, colon) varies by cancer type. This plan tests, in cohorts that played no part in
any of it, whether the cancer-type pattern for survival replicates, and whether the inflamed genes add to response in
an independent checkpoint cohort.

## What is known when this is written

* Every result of manuscript v13, including the TCGA per-type values that serve as predictions below
  (icb/v15_content/atlas_by_type.csv, rho_inflamed per type).
* Sample annotations (GEO query pages, brief view: no expression values) of the candidate series listed in
  external/meta/screen_b1.txt and screen_b2.txt, read to apply the rule below.
* GSE84437 (gastric) was analyzed in the separate barriernet work for the prognostic value of a 26-gene immune score
  that shares 10 of its genes with the inflamed signatures; the inflamed signatures themselves and any decoys have never
  been scored in it here. GSE96058 (SCAN-B) was used in round 7 for the breast library only. GSE176307 was scored in
  September 2026 (validation/PREREGISTRATION_external_validation.md) for every signature's AUC against a size- and
  engine-matched null; its decoys have never been grown or scored. No other cohort below has been used in this project.

## A. Cohorts

Survival. Candidates were the large public non-TCGA cohorts known to the authors for each of the 20 TCGA cancer types
(external/meta/screen_*.txt; not a systematic search, and stated as such). Rule: tumors of that type (primary unless
the type, as for melanoma in TCGA, is mostly metastatic), not treated with checkpoint blockade, expression in GEO on a
platform measuring at least 15,000 genes after mapping probes to symbols (the probe with the largest mean kept), a
survival time and event in the GEO sample annotations, at least 80 patients and 25 events; endpoint overall survival if
annotated, otherwise disease-specific, then relapse-, disease- or metastasis-free survival. Among eligible series of a
type, the one with most events; for stomach, the type the manuscript's Discussion turns on, every eligible series.

| type (TCGA) | series | platform | endpoint | patients / events (annotations) |
|---|---|---|---|---|
| STAD | GSE84437 | GPL6947 | OS | 433 / 209 |
| STAD | GSE26253 | GPL8432 | RFS | 432 / 177 |
| SKCM | GSE65904 | GPL10558 | DSS | ~210 / 102 |
| LUAD | GSE72094 | GPL15048 | OS | ~398 / 122 |
| LUSC | GSE157011 | GPL570 | OS | 484 / 264 |
| COAD | GSE39582 | GPL570 | OS | tumors ~566 / 194 |
| CESC | GSE44001 | GPL14951 | DFS | 300 / 38 |
| SARC | GSE21050 | GPL570 | MFS | 309 / 121 |
| HNSC | GSE65858 | GPL10558 | OS | 270 / 94 |
| BLCA | GSE13507 | GPL6102 | OS | primary 165 / 69 |
| OV | GSE32062 | GPL6480 | OS | 260 / 121 |
| ESCA | GSE53625 | GPL18109 | OS | tumors 179 / ~106 |
| BRCA | GSE96058 (SCAN-B) | RNA-seq | OS | 3,273 / 336 (round 7 preparation) |

Not eligible among the candidates: KIRC, KIRP, LGG, GBM, LIHC, PAAD, MESO, UCEC (no series meeting the rule); for the
types above, the other candidates had fewer events, fewer than 15,000 genes (U133A series) or no survival in GEO.
A series that fails the gene, patient or event rule after download is dropped and reported; none is replaced.

Response. GSE176307 (metastatic urothelial carcinoma, checkpoint blockade; 89 evaluable, 16 CR/PR; Salmon TPM, as
prepared in September 2026), alone.

## B. Procedure

As in rounds 7 and 10, unchanged: the 52 signatures by their own engines (gene symbols), decoys grown by decoyfdr
(frozen install icb/v12_pkg/Rlib6) from the genes measured in every survival cohort and in no curated signature, ten
sets of 32 co-expressed and 8 random-gene decoys (seed 20261101), 1,000-draw exact-size banks per target and cohort
(seed 20261102), folded C-index; for GSE176307 its own decoys (seeds 20261103, 20261104), folded AUC. rho, structure,
content and content calls as defined in round 10 (decoys/v15_content_structure.py functions). Script:
applications/v16_external.R; analysis decoys/v17_external_replication.py; outputs icb/v16_external/, icb/v17_external/.

## C. Predictions and tests

1. Primary. Over the 12 cancer types, the Spearman correlation between the inflamed family's mean rho in the external
   cohort(s) of the type (stomach: the mean over its two cohorts) and in TCGA is positive. One-sided P from 100,000
   random permutations of the external values over types (seed 916000). Success: P < 0.05.
2. Stomach: the inflamed family's mean rho is below 0.5 in each gastric cohort (TCGA 0.19). Melanoma: above 0.5
   (TCGA 0.95). Each is reported as holding when the point estimate lies on the predicted side and as confirmed when
   its 95% interval over the 13 correlation clusters (10,000 resamples, seed 916001) excludes 0.5.
3. Response: the inflamed family's mean rho in GSE176307 is above 0.5 (response cohorts 0.74; IMvigor210 0.76), with
   the same interval; and in bladder cancer it is higher for response (GSE176307) than for survival (GSE13507), the
   difference with a cluster interval (TCGA BLCA 0.53 against IMvigor210 0.76).
4. Library level, descriptive: content calls by family for survival over the external cohorts (TCGA: inflamed 1 of
   20, confounder 5 of 8), the decoy estimate of the survival calls with its pi0-corrected form, family means of rho,
   per-cohort single-cohort estimates, and the TCGA-against-external scatter by type.

Every result is reported whatever it shows. The endpoints differ from TCGA's overall survival in five cohorts, which
is a stated limitation, not a reason to drop them.

## D. Displays and text

A new panel set in Figure 8 or a supplementary figure (external against TCGA per type; per-cohort family rho; the
response cohort), a Results paragraph, a Methods paragraph, Additional file 2 sheets, and, if space allows, a clause in
the abstract. Numbers already reported do not change.
