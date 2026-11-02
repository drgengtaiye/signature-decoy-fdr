# Analysis plan, round 13: the false discovery rate of published prognostic risk models, and of the pipeline that builds them

Written 2026-11-02, after manuscript v15 and before any paper below was retrieved or any risk model or theme pool was
scored in any cohort. Its SHA-256 is recorded in PLAN_round13_published_risk_models.sha256 and the file is deposited
before the first retrieval.

## Question

A large literature builds a "[theme]-related gene signature" in a TCGA cancer type (genes of a theme such as
ferroptosis or hypoxia; univariate Cox, then LASSO-Cox; a risk score that is a weighted sum of expression) and reports
that it validates in an independent cohort. Part A asks, for a random sample of such published models, what share of
those that validate in an independent cohort would be matched by modules of non-theme genes with the same structure
and the same weights: the decoy false discovery rate of the claim "this model validates". Part B asks the same of the
pipeline itself, run as published on theme gene pools and on matched decoy pools.

## What is known when this is written

* Every result of manuscript v15, including the 52 immune signatures in the 11 independent survival cohorts of round
  11 (external/cohorts/manifest.csv, included rows) and in 20 TCGA types.
* No published risk model and no theme pool below has been scored in any cohort in this project. No paper has been
  retrieved. The search below has not been run.

## Cohorts

The 11 independent cohorts of round 11, one per TCGA type: STAD GSE84437, SKCM GSE65904, LUAD GSE72094, LUSC
GSE157011, COAD GSE39582, CESC GSE44001, SARC GSE21050, HNSC GSE65858, BLCA GSE13507, OV GSE32062, BRCA GSE96058, with
the endpoints prepared in round 11; and the TCGA primary tumors of each type with TCGA-CDR overall survival
(applications/common.R, tcga_cohorts), the training data of the literature.

## A. Published models

A1. Search. PubMed Central through E-utilities (esearch, db = pmc), run once on the day of the deposit, per type:
`(TCGA OR "The Cancer Genome Atlas") AND ("risk score" OR "risk model" OR "prognostic signature" OR "prognostic model")
AND (<type terms>) AND (ferroptosis OR cuproptosis OR pyroptosis OR necroptosis OR anoikis OR disulfidptosis OR
autophagy OR m6A OR "N6-methyladenosine" OR hypoxia OR glycolysis OR senescence OR "immune-related" OR mitophagy OR
"epithelial-mesenchymal") AND "open access"[filter] AND 2019:2025[pdat]`, with type terms STAD "gastric cancer" OR
"stomach adenocarcinoma" OR "gastric adenocarcinoma"; SKCM melanoma; LUAD "lung adenocarcinoma"; LUSC "lung squamous";
COAD "colon cancer" OR "colon adenocarcinoma" OR "colorectal cancer"; CESC "cervical cancer" OR "cervical squamous";
SARC sarcoma; HNSC "head and neck squamous" OR HNSCC; BLCA "bladder cancer" OR "urothelial carcinoma"; OV "ovarian
cancer"; BRCA "breast cancer". The returned PMCIDs and counts are saved (validation/round13_search/).

A2. Order. Per type, the PMCIDs sorted ascending and shuffled with numpy default_rng(20261113 + type index, types in
the order above); screened in that order until 10 eligible models are found or 80 records screened.

A3. Eligibility. (i) Original research that builds a prognostic risk score from tumor gene expression with TCGA tumors
of the type as the training set (alone or merged with other data); (ii) the score is given in the PMC full text (body,
table or figure caption; supplementary files are not read) as a linear formula, each gene with its coefficient;
(iii) 2 to 50 genes, all protein-coding (any lncRNA, miRNA, methylation, copy-number or protein term makes the model
ineligible); (iv) a survival endpoint. One model per paper: the first presented as the main model. Recorded: PMCID,
PMID, year, journal, theme, genes and coefficients as printed, the verbatim formula text, training data, every
validation cohort named, endpoint, and the reason for each ineligible record.

A4. Extraction check. Genes and coefficients are re-parsed by script from the recorded verbatim text; mismatches are
resolved against the full text. A random 20% of eligible models (seed 20261114) is re-extracted in a second,
independent pass; agreement is reported. The curated table validation/published_models.csv is frozen with its SHA-256
in an amendment before any model is scored.

A5. Genes. Symbols are mapped to current symbols through the org.Hs.eg.db alias table (breast/export_orgHs_symbols.R).
A model is scored on its genes measured in both cohorts of its type if they are at least 80% of its genes, and is
otherwise excluded and reported.

A6. Scores. Score = sum over genes of coefficient times the gene's z-score within the cohort (expression on the log
scale as stored), in the external cohort of the type and in TCGA of the type; higher scores mark shorter survival.

A7. Decoys. Per type, decoyfdr (frozen install icb/v12_pkg/Rlib6) grow_decoys over the two cohorts of the type jointly,
from genes measured in both and in no eligible model of that type and no theme pool of part B: 32 co-expressed and 8
random-gene decoys per model, one set, seed 20261115 + type index. A decoy carries its model's coefficients in random
order (the package's rule for weighted-mean signatures). Because the published direction comes from fitting in TCGA,
each decoy's direction is set the same way: the sign of the association of its score with survival in TCGA
(C-index above or below one half).

A8. Statistics, per model or decoy in the external cohort: (L) the literature's: the score split at its median,
log-rank P < 0.05 with the high-score group faring worse ("validated"); (C) the directed C-index (higher score, shorter
survival); rho, the share of the model's co-expressed decoys whose directed C-index it exceeds (ties one half).

A9. Estimates. FDR_L = sum over models of d_m divided by V, with d_m the share of model m's co-expressed decoys
validated under (L) and V the number of models validated; the same with random-gene decoys. Intervals: 10,000
bootstrap resamples of models with their decoys, stratified by type (seed 20261116).

A10. Population and tests. Primary population: eligible models whose paper does not name the external cohort of
its type among its validation cohorts; if fewer than 40, all eligible models (stated). The other population is
reported beside it. H1: models validate more often than their co-expressed decoys (V/M against mean d_m; one-sided
bootstrap P). Primary estimand: FDR_L with its interval. Secondary: mean rho with interval; FDR_L by theme family and
by type; validation rates in TCGA (in-sample) for models and decoys; the same with (C) against the decoys' 95th
percentile in place of (L).

## B. The pipeline

B1. Theme pools (msigdbr 26.1.0 unless listed): ferroptosis WP_FERROPTOSIS; pyroptosis REACTOME_PYROPTOSIS;
necroptosis GOBP_NECROPTOTIC_PROCESS; anoikis GOBP_REGULATION_OF_ANOIKIS; autophagy KEGG_REGULATION_OF_AUTOPHAGY;
mitophagy REACTOME_MITOPHAGY; senescence FRIDMAN_SENESCENCE_UP; hypoxia HALLMARK_HYPOXIA; glycolysis
HALLMARK_GLYCOLYSIS; EMT HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION; inflammatory HALLMARK_INFLAMMATORY_RESPONSE;
lactate GOBP_LACTATE_METABOLIC_PROCESS; m6A (23 regulators: METTL3, METTL14, METTL16, WTAP, VIRMA, RBM15, RBM15B,
ZC3H13, FTO, ALKBH5, YTHDC1, YTHDC2, YTHDF1, YTHDF2, YTHDF3, IGF2BP1, IGF2BP2, IGF2BP3, HNRNPA2B1, HNRNPC, FMR1,
LRPPRC, ELAVL1); cuproptosis (19: NFE2L2, NLRP3, ATP7B, ATP7A, SLC31A1, FDX1, LIAS, LIPT1, LIPT2, DLD, DLAT, PDHA1,
PDHB, MTF1, GLS, CDKN2A, DBT, GCSH, DLST); disulfidptosis (10: GYS1, NDUFS1, OXSM, LRPPRC, NDUFA11, NUBPL, NCKAP1, RPN1,
SLC3A2, SLC7A11). A set name absent from msigdbr 26.1.0 is replaced by the nearest set of the same source and theme,
recorded before any pool is scored.

B2. Pipeline, per (pool, type), on the pool's genes measured in TCGA and the external cohort: univariate Cox on TCGA
overall survival, genes with P < 0.05 kept; if fewer than two, no model. LASSO-Cox (glmnet, alpha 1, cv.glmnet with 10
folds, fold ids drawn with seed 20261117 + unit index, lambda.min); if fewer than two genes have nonzero coefficients,
no model. Score as in A6 with the fitted coefficients; validated as in A8 (L).

B3. Decoy pools, per (pool, type): 20 co-expressed and 5 random-gene decoys grown as in A7 over the two cohorts of the
type (seed 20261118 + type index), each run through B2 unchanged with its own fold ids by the same rule.

B4. Estimates: per unit, model made and validated; FDR_B = sum over units of the decoy validation rate divided by the
number of theme units validated; the share of decoy pools that yield a validated model; bootstrap over units
stratified by type (10,000, seed 20261119). H2: theme pools validate more often than their decoy pools (one-sided
bootstrap P).

## C. Reporting

Every result is reported whatever it shows. Endpoints differ from TCGA overall survival in five external cohorts
(round 11), a stated limitation. Scripts: applications/v19_published_models.R, applications/v19_pipeline.R,
decoys/v19_audit.py; outputs icb/v19_audit/. Displays: a main figure, a Results paragraph, a Methods section,
Additional file 2 sheets, and the abstract if space allows. Numbers already reported do not change.
