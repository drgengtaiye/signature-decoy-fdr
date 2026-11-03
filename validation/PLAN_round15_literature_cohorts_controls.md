# Analysis plan, round 15: the whole search of round 13, the validation cohorts the literature uses, a law for how often other genes validate, and clinical positive controls

Written 2026-11-03, after manuscript v16 and before any record below was screened, any cohort below was downloaded, and
any model, null model or assay below was scored. Its SHA-256 is recorded in
PLAN_round15_literature_cohorts_controls.sha256 and the file is deposited before the first step. Changes are made only
by dated amendments deposited before the step they change.

## What is known when this is written

* Every result of manuscript v16, including round 13 (PLAN_round13_published_risk_models.md, amendments 1-2): 1,892
  PubMed records in 11 types, 294 screened in the frozen random order, 109 eligible; of 60 independent models 19
  validated, decoys 12%, FDR 39% (29-57%); the pipeline on 15 theme pools FDR 77% (64-94%); 69-77% of decoy pools
  yielded a validated model in LUAD and BRCA, 3-22% in the other nine types.
* The validation cohorts named by the 294 screened records were recorded in round 13 (column validation_cohorts) but
  have not been tabulated. None of the 1,598 unscreened records has been read. No cohort outside round 11 has been
  downloaded for this purpose. No clinical assay of part D has been scored in this project outside the breast library
  of v15 (the van 't Veer and Paik lists of Venet et al.).
* The bulk GEO mirror (ftp.ncbi.nlm.nih.gov) does not answer from the cluster; the GEO query pages
  (www.ncbi.nlm.nih.gov/geo/query/acc.cgi, sample tables) do.

## A. The whole search

A1. Records: the 1,892 records of round 13; the search is not re-run. Each type's records are screened in the frozen
order (round13_search/<TYPE>_order.tsv, SHA256SUMS_search_v2), continuing after the last record screened in round 13,
to the end of the list; no cap. Eligibility, extraction, recording and gene mapping are those of round 13 (A3-A5 and
its amendments) unchanged; round-13 decisions stand. For every record that trains a prognostic model in TCGA of the
type, eligible or not, the validation cohorts it names are recorded (for part B).

A2. A random 10% of the newly eligible models (numpy default_rng(20261201)) is re-extracted in a second, independent
pass; agreement is reported. The table validation/published_models_r15.csv (round-13 models and new ones) is frozen with
its SHA-256 in an amendment before any new model is scored.

A3. Decoys and statistics: round 13 A6-A9 unchanged (32 co-expressed and 8 random-gene decoys per model, grown over TCGA
and the round-11 cohort of the type from genes in no eligible model of the type and no theme pool), on the full table,
seed 20261202 + type index. The round-13 numbers already reported do not change; the round-15 estimates are reported
beside them.

A4. Primary: FDR_L among independent models (paper names neither the round-11 cohort of its type nor a series
containing it) with a 95% interval from 10,000 bootstrap resamples of models with their decoys, stratified by type
(seed 20261203). Secondary: the same for all models; by type; by theme family; by publication year (2019-2021,
2022-2025); by number of genes (2-5, 6-10, 11+); by training-set size (tertiles). Exploratory: citations per year
(PubMed Central cited-by counts, retrieved once on one day) for models with rho >= 0.95 against rho < 0.95.

## B. The validation cohorts the literature uses

B1. Candidates: per type, the cohorts named as validation cohorts by at least 5 records of A1 (round 13 and round 15
together), or among the 5 most named of the type if fewer reach 5, ties included; plus the round-11 cohort. TCGA is
excluded (training data).

B2. Eligibility: the round-11 rule (tumors of the type; no checkpoint blockade; at least 15,000 genes after mapping
probes to symbols, probe with the largest mean kept; survival time and event in the annotations; at least 80 patients
and 25 events; overall survival, else disease-specific, else relapse-, disease- or metastasis-free survival), except
that every eligible cohort is kept, not only the one with most events; and retrievable, that is GEO sample tables
through the query pages, a public HTTPS source (ICGC, ArrayExpress, cBioPortal), or already on disk. A candidate that
fails is listed with the reason; none is replaced.

B3. Null models: per type, for each of the 15 theme pools of round 13 (B1 there), 6 co-expressed and 2 random-gene
decoy pools grown with decoyfdr over TCGA alone from genes in no theme pool and no eligible model of the type (seed
20261204 + type index), 120 pools per type, each run through the round-13 pipeline (univariate Cox P < 0.05, then
LASSO-Cox, lambda.min, fold ids seed 20261205 + pool index) fitted in TCGA alone. Each fitted model is a null model:
a model the pipeline makes from genes of no theme. It is scored in every eligible cohort of its type on its genes
measured there (at least 80%, otherwise not scorable there), as in round 13 A6.

B4. The false-validation rate of a cohort, FVR_c: the share of the type's null models, among those scorable in c,
validated in c under the round-13 rule (L) (median split, log-rank P < 0.05, high score faring worse). Beside it, the
validation rate in c of the published models of the type whose papers did not name c.

B5. H4 (primary for B): within a type, the cohorts' FVRs differ. Cochran's Q over the null models scorable in every
eligible cohort of the type, for each type with at least two cohorts; the Qs and their degrees of freedom are summed
over types for one P. Secondary: Spearman correlation of FVR_c with the cohort's events (one-sided, positive); FVR_c
by type.

B6. H5 (primary for replication): a decoy check in the cohort where a model validated predicts whether it replicates in
another. Units: a published model with an ordered pair (A, B) of eligible cohorts of its type such that it validated
(L) in A and its paper did not name B. rho_A is the share of its co-expressed decoys (A3) whose directed C-index in A
it exceeds. Statistic: replication rate in B for rho_A >= 0.95 minus that for rho_A < 0.95, each model weighted
equally across its pairs; one-sided P from 10,000 bootstrap resamples of models stratified by type (seed 20261206).
Secondary: the replication rate in B of null models and of decoys validated in A.

## C. A law for how often other genes validate

C1. In each eligible cohort c, the genes measured are z-scored and the first K = 10 principal components computed;
survival (the cohort's endpoint of B2) is regressed on the 10 components jointly by Cox regression, giving
gamma_c. For each null model m scorable in c, its standardized score is regressed on the components (loadings a_mc)
and its predicted log hazard ratio per SD is beta_mc = s_m a_mc . gamma_c, with s_m its direction from TCGA. The
predicted probability that it validates under (L) is p_mc = Phi(0.80 sqrt(D_c) beta_mc - 1.96), D_c the cohort's events
(0.80 = the mean difference between the halves of a standard normal score, 1.596, divided by 2), and the predicted
FVR_c is the mean of p_mc.

C2. H6 (primary for C): the predicted FVR_c tracks the observed one across cohorts (Spearman rho > 0, one-sided
permutation P, 10,000 permutations, seed 20261207), and its mean absolute error is below that of an events-only
prediction (leave-one-cohort-out logistic regression of observed FVR on sqrt(D_c) over the other cohorts), one-sided
paired bootstrap over cohorts (seed 20261208). Secondary: K = 5 and K = 20; the area under the curve of p_mc for the
observed validation of each null model, pooled.

## D. Clinical positive controls

D1. Assays: multigene prognostic assays named in ASCO or NCCN guidelines whose genes and algorithm are published:
breast, the 21-gene Recurrence Score, the 70-gene signature, EndoPredict, the PAM50 risk of recurrence and the Breast
Cancer Index; colon, the 12-gene colon Recurrence Score. Genes, weights and direction are taken from the original
publications, with the page or table recorded, before any is scored; where the published algorithm is not a weighted
sum of genes (centroid correlations, ratios), it is approximated by one, stated per assay; reference genes are left
out. An assay whose genes are not printed in a retrievable publication is dropped and reported.

D2. Cohorts: the eligible breast and colon cohorts of part B (including GSE96058 and GSE39582), METABRIC (on disk) and
TCGA-BRCA and TCGA-COAD (none of these assays was trained on TCGA). Decoys: 32 co-expressed and 8 random-gene per assay,
grown over TCGA and the round-11 cohort of the type from genes in no assay, no eligible model and no theme pool (seed
20261209 + assay index), scored and judged as in round 13 A6-A9. All patients; ER-positive patients as a secondary
analysis where annotated.

D3. H7 (primary for D): the assays' decoy FDR of "validates" (round-13 estimator, pooled over assays and cohorts) is
below that of the independent published models of A4; one-sided P from 10,000 bootstrap resamples, assays and
published models resampled separately (seed 20261210). Secondary: each assay's rho per cohort.

## E. Reporting

Every result is reported whatever it shows; a hypothesis whose test cannot be run (too few cohorts or units) is
reported as not run, with the reason. Endpoints that are not overall survival are a stated limitation. Scripts:
applications/v20_*.R, decoys/v20_*.py; outputs icb/v20_round15/. Numbers already reported do not change. Displays: one
main figure for parts B-C, panels added to Fig. 9 for parts A and D, Supplementary Data sheets, Methods.
