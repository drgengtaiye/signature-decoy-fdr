# Venet 2011 signature library: provenance against our 13 breast cancer cohorts

This file goes with `venet_provenance.csv`, which has one row for each of the 45 signatures. The coding comes only from the literature. No expression or outcome data were analyzed, and no analysis results were consulted. The Venet gene-list file (`library/venet_cancer_signatures.csv`, gene symbols only) was used only to check which list in each paper Venet took.

Sources: Europe PMC and PMC full text, publisher pages, and the papers' own supplementary files. All quotes in the CSV are verbatim and give the section and URL. Venet et al. 2011 Text S1 (https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1002240.s003) was used only for its bibliographic entries: title, journal, PMID and gene count for each signature.

## 1. Coding rule

The rule is applied to each cohort as a set of patients, through the established overlaps:

- **D**: the outcome data of the cohort, or of patients it shares with other datasets, were used to pick genes or to fit weights, cut-offs or a model.
- **C**: the gene set was fixed without this cohort's outcomes, but analysis of the cohort's outcome prompted or shaped the signature. In practice this is the paper's own discovery cohort.
- **V**: a finished signature was only tested against outcome in the cohort, as validation or support for a claim in the original paper.
- **blank**: the cohort and its overlapping patients were not used in the original paper. Uses of expression without outcome are blank but are described in the evidence column.

Only the original defining publication counts. Later validations by others, including Venet 2011 for the published signatures, do not.

### How the rule was applied

1. **Overlaps used:**
   - van 't Veer 2002 (78 training plus 19 validation, 117 total) or van de Vijver 2002 (NKI 295) → NKI.
   - The Uppsala cohort (Miller 2005 / GSE3494, Pawitan 2005, Ivshina 2006 / GSE4922, Calza 2006) → both Loi_GSE6532 (about 117 shared patients) and GSE4922.
   - GSE2990 (Sotiriou 2006: KJX64 and KJ125) → Loi_GSE6532 only.
   - The Stockholm cohort (Pawitan 2005) → GSE1456.
   - TRANSBIG / Desmedt 2007 → GSE7390 and Loi_GSE6532 (41 shared patients).
   - Symmans 2010 → GSE17705 and Loi_GSE6532 (Oxford tamoxifen patients).
2. **Dates:** a cohort first published after the signature paper is blank. TCGA, METABRIC, GSE42568 and GSE96058 are blank for all 45 signatures. GSE11121, GSE17705, GSE20685 and GSE21653 were in principle available to the 2009–2011 papers, but none of those papers used them.
3. **What counts as D:** fitting a cut-off on outcome counts, as the rule says. Examples are the log-rank-optimized threshold in GLINSKY and the X-tile cut-points in HUA. Refitting a classifier on a cohort's outcome labels also counts, as with the leave-one-out centroids in KORKOLA. These D codes are listed in Section 5 so that a lenient sensitivity analysis can recode them as V.
4. **What counts as C:** choosing among candidate gene sets, scoring rules or contrasts after looking at the cohort's outcome (WELM, WONG-MITOCHON).
5. **Clustering or filtering inside the cohort:** if the cohort's expression was used without outcome (clustering, variance filtering, module discovery, or non-outcome labels such as grade, p53 status or subtype) before an outcome test in the same patients, the code is V and the case is listed as borderline (Section 4).
6. **Cut-offs fitted on non-breast or non-study cohorts** do not change our codes. These are the prostate data for GLINSKY, Wang/EMC for YU and SHIPITSIN, MSKCC/EMC for TAVAZOIE and the GSE2603 data for MORI. They make `outcome_guided` = yes.
7. **Parent-level harmonization:** the per-paper curation was done in six parallel batches and then harmonized. One code was changed in that step: KORKOLA NKI went from V to D, because the paper's Methods state that the classifier was evaluated by leave-one-out against the outcome means of the remaining tumors. This is a model fitted on the van 't Veer outcome labels.

## 2. Summary matrix (all 45 signatures; blank cells mean not used)

| Signature | PMID | NKI | Loi_GSE6532 | TCGA_BRCA | GSE7390 | GSE11121 | GSE1456 | GSE4922 | GSE20685 | GSE21653 | GSE42568 | GSE17705 | GSE96058 | METABRIC | Outcome-guided | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ABBA | 21082037 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| BEN-PORATH-EXP1 | 18443585 | V | V |  | V |  | V | V |  |  |  |  |  |  | no | high |
| BEN-PORATH-PRC2 | 18443585 | V | V |  | V |  | V | V |  |  |  |  |  |  | no | high |
| BUESS | 17868458 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| CARTER | 16921376 | V | V |  |  |  |  |  |  |  |  |  |  |  | no | medium |
| CHANG | 14737219 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| CHI | 16417408 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| CRAWFORD | 18427120 | V | V |  |  |  | V | V |  |  |  |  |  |  | no | high |
| DAI | 15899795 | D |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| GLINSKY | 15931389 | D |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| HALLSTROM | 18167336 |  | V |  |  |  | V | V |  |  |  |  |  |  | no | high |
| HE | 20028755 | D | D |  | V |  | V | D |  |  |  |  |  |  | yes | high |
| HU | 19111877 | D |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| HUA | 19563758 | D | D |  |  |  |  | D |  |  |  |  |  |  | yes | medium |
| IVSHINA | 17079448 | V | V |  |  |  | V | V |  |  |  |  |  |  | no | high |
| KOK | 19353633 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | low |
| KORKOLA | 17428335 | D |  |  |  |  |  |  |  |  |  |  |  |  | yes | medium |
| LIU | 17229949 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | medium |
| MA | 12714683 |  |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| META-PCNA | 22028643 | V | V |  |  |  |  |  |  |  |  |  |  |  | no | high |
| MILLER | 16141321 | V | V |  |  |  |  | V |  |  |  |  |  |  | no | high |
| MORI | 19483725 |  |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| PAIK | 15591335 |  |  |  |  |  |  |  |  |  |  |  |  |  | yes | medium |
| PAWITAN | 16280042 | V | V |  |  |  | D | V |  |  |  |  |  |  | yes | high |
| RAMASWAMY | 12469122 | V |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| REUTER | 19477427 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| RHODES | 15184677 |  |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| SAAL | 17452630 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| SHIPITSIN | 17349583 | V | D |  |  |  |  |  |  |  |  |  |  |  | yes | medium |
| SORLIE | 12829800 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | medium |
| SOTIRIOU-93 | 12917485 |  |  |  |  |  |  |  |  |  |  |  |  |  | yes | medium |
| SOTIRIOU-GGI | 16478745 | V | V |  |  |  |  |  |  |  |  |  |  |  | no | high |
| TAUBE | 20713713 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| TAVAZOIE | 18185580 |  |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| VALASTYAN | 19524507 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| VANTVEER | 11823860 | D |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| WANG-76 | 15721472 |  |  |  |  |  |  |  |  |  |  |  |  |  | yes | medium |
| WANG-ALK5T204D | 18625725 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| WELM | 17456594 | C |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| WEST | 15869330 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| WHITFIELD | 12058064 |  |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| WONG-ESC | 18397753 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| WONG-MITOCHON | 18199530 | C |  |  |  |  |  |  |  |  |  |  |  |  | yes | high |
| WONG-PROTEAS | 18199530 | V |  |  |  |  |  |  |  |  |  |  |  |  | no | high |
| YU | 18006806 | V | V |  |  |  | V | V |  |  |  |  |  |  | yes | high |

Totals:

| Code | Entries | Signatures | By cohort |
|---|---|---|---|
| D | 13 | 9 | NKI 7 (VANTVEER, DAI, GLINSKY, HU, HE, HUA, KORKOLA); Loi_GSE6532 3 (HE, HUA, SHIPITSIN); GSE4922 2 (HE, HUA); GSE1456 1 (PAWITAN) |
| C | 2 | 2 | NKI (WELM, WONG-MITOCHON) |

No signature is D or C for GSE7390, GSE11121, GSE17705, GSE20685, GSE21653, GSE42568, GSE96058, TCGA_BRCA or METABRIC.

## 3. All D and C entries, with the key quote

| Signature | Cohort | Code | Key quote | Section | Source |
|---|---|---|---|---|---|
| DAI | NKI | D | "Prognostic features were selected based on the training samples by their correlation to outcome and were reselected during each step of leave-one-out cross-validation. ... Within the ER/age high group, we identified a group of 50 prognostic reporter genes that were highly correlated with the outcome" | Materials and Methods; Results | https://aacrjournals.org/cancerres/article/65/10/4059/517951 |
| GLINSKY | NKI | D | "We used the training set to select the prognosis-discrimination cutoff value for a signature based on the highest level of statistical significance in patients’ stratification into poor- and good-prognosis groups as determined by the log-rank test (lowest P value and highest hazard ratio in the training set). ... The cohort of 97 breast cancer patients ... was divided into a training set (78 patients ...) and a test set (19 patients ...) as described in the original publication ( 33 )." | Methods | https://pmc.ncbi.nlm.nih.gov/articles/PMC1136989/ |
| HE | NKI | D | "intercohort multivariate Cox proportional-hazards analysis of each signature gene was performed in three breast cancer datasets. Genes significantly correlated with patient outcomes (p<0.05) were determined for each datasets. ... 6 genes were predictive in all three datasets. This served as the logic and derivation of the new 6-gene model." | Materials and Methods: Selecting the six-gene model from the SpMGS | https://pmc.ncbi.nlm.nih.gov/articles/PMC2866744/ |
| HE | GSE4922 | D | "These included the van de Vijver, GSE4922 and GSE2034 gene sets. [...] Eighty six patients which overlapped with the GSE2990 dataset were excluded, leaving 163 patients in this analysis." | Results; Materials and Methods: GSE4922 dataset | https://pmc.ncbi.nlm.nih.gov/articles/PMC2866744/ |
| HE | Loi_GSE6532 | D | "Eighty six patients which overlapped with the GSE2990 dataset were excluded, leaving 163 patients in this analysis." | Materials and Methods: GSE4922 dataset (D via Uppsala patients in Loi that are not in GSE2990; GSE2990 itself was validation only) | https://pmc.ncbi.nlm.nih.gov/articles/PMC2866744/ |
| HU | NKI | D | "Poor-prognosis genomic gain at 8q22 was detected in all three datasets ( van 't Veer et al., 2002 ; van de Vijver et al., 2002 ; Wang et al., 2005 )." | Figure 1 legend; Experimental Procedures (ACE expression score from poor- vs good-prognosis phenotypes) | https://pmc.ncbi.nlm.nih.gov/articles/PMC2676231/ |
| HUA | NKI | D | "The optimal two cutting-points of RA signature scores were determined by X-tile analysis (Camp et al., 2004)." | Supplemental Experimental Procedures: Analysis of RAR binding and RA expression profiling on breast tumor samples; applied to the 295 NKI patients in Figure 6C-D | https://ars.els-cdn.com/content/image/1-s2.0-S0092867409005042-mmc1.pdf |
| HUA | GSE4922 | D | "The optimal two cutting-points of RA signature scores were determined by X-tile analysis (Camp et al., 2004). [Figure S14: 249 (Ivshina et al., 2006) ... GSE4922]" | Supplemental Experimental Procedures and Figure S14 legend | https://ars.els-cdn.com/content/image/1-s2.0-S0092867409005042-mmc1.pdf |
| HUA | Loi_GSE6532 | D | "Kaplan-Meier curves of relapse-free survival for two breast cancer cohorts which contain 249 (Ivshina et al., 2006) ... GSE4922 (A)" | Figure S14 legend (via the ~117 Uppsala patients shared with Loi) | https://ars.els-cdn.com/content/image/1-s2.0-S0092867409005042-mmc1.pdf |
| KORKOLA | NKI | D | "Classification of tumors for all gene sets was performed using a leave one out cross validation analysis based on correlation of each tumor's gene expression with the means of the remaining "good" or "bad" outcome tumors. (Methods) + We tested the performance of our 21 gene set in the tumor sets of van't Veer et al. and Sotiriou et al. using PAM classification. The gene set classified 65% of the specimens correctly in the van't Veer tumor set (Results)" | Methods, Statistical analysis; Results, Comparison to published predictive gene sets | https://pmc.ncbi.nlm.nih.gov/articles/PMC1855059/ |
| PAWITAN | GSE1456 | D | "In order to maximize the statistical power, we initially used the whole Stockholm cohort (n = 159) as the training set for choosing an optimal gene set and identifying risk groups in the hierarchical cluster analysis." | Materials and methods, Primary analysis and validation data sets | https://pmc.ncbi.nlm.nih.gov/articles/PMC1410752/ |
| SHIPITSIN | Loi_GSE6532 | D | "Univariate Cox proportional hazards regression was carried out on the data for the probes selected from Data Sets 1 and 2 ... to find probes that tend to correlate with distant metastasis-free survival time ... Probes with a Cox regression p-value < 0.2 ... in both data sets were analyzed further. [Data Set 2 = Sotiriou et al., 2006, GSE2990]" | Supplemental Data, Identification of signatures with prognostic relevance | https://ars.els-cdn.com/content/image/1-s2.0-S1535610807000293-mmc1.pdf |
| VANTVEER | NKI | D | "We calculated the correlation between the prognostic category (metastasis vs. no-metastasis) and the logarithmic expression ratio across all 78 samples for each individual gene in the 5,000 significant genes. ... The combined error rate reaches the minimum when we use 70 marker genes from the top of our candidate list." | Supplementary Information, Method of supervised classification | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2F415530a/MediaObjects/41586_2002_BF415530a_MOESM7_ESM.doc |
| WELM | NKI | C | "We found that MSP gene expression alone was not prognostic for outcome in human breast cancer (data not shown); however, ... expression of all three genes ( MSP , MT-SP1 , and MST1R ) would be required for activation of the MSP pathway. ... overexpression of MSP , MT-SP1 , or MST1R alone, or any combination of two of these three genes, did not reveal any significant association with clinical outcome in this data set (data not shown)." | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC1855278/ |
| WONG-MITOCHON | NKI | C | "The mitochondria module was not significantly associated with survival across all 295 samples; however, the 50 tumors that had the highest average expression of the mitochondria module had a statistically significant association with poor survival compared with the 50 tumors that had the lowest average expression ( P < 0.02, Fig. 2 B )." | Results | https://aacrjournals.org/cancerres/article/68/2/369/541675 (doi:10.1158/0008-5472.CAN-07-0382) |

Notes on the D/C entries:

- **VANTVEER:** coded from the defining Nature paper (PMID 11823860). NKI contains 61 of its 78 training patients. The NEJM paper that Venet cites validated the finished profile on the 295 patients, which would be V on its own. D takes precedence.
- **DAI:** the 50 reporters were selected on metastasis outcome within the ER/age-high group of the NKI tumors (311 tumors: the 295 plus 16 from van 't Veer). There was no external validation.
- **HU:** the 8q22 region was retained because it was a poor-prognosis gain, by ACE expression scores computed from poor- versus good-prognosis tumors, in van 't Veer, van de Vijver and Wang. Venet's 13 genes are the known genes in that region.
- **HE:** the 6 genes were selected by per-gene Cox analysis in NKI 295, GSE4922 (163 Uppsala patients after removing 86 shared with GSE2990) and GSE2034. GSE1456, GSE2990 and GSE7390 were validation sets.
  - Loi_GSE6532 = D rests on overlap arithmetic: about 117 Uppsala patients are in Loi, but only 86 were removed, so about 30 Loi patients were probably in He's 163.
  - A sample-level check against He's patient list would settle this. If none are shared, recode Loi as V.
- **HUA:** the gene set comes from MCF-7 cells. Only the two cut-points of the RA score were optimized on outcome with X-tile. The X-tile sentence is generic. Its use for GSE4922 (all 249 patients, hence Loi through the Uppsala overlap) is inferred from the three-group split shown in Figure S14.
- **SHIPITSIN:** probes had to have Cox P < 0.2 for distant metastasis-free survival in both Wang (GSE2034) and Sotiriou 2006 (GSE2990, 125 node-negative untreated patients, all of them in Loi). NKI (141 node-negative untreated) was used for an expression-variability filter (no outcome) and then for validation.
  - Matching Venet's 56 genes to signature A rests on gene content; main-text Table 1 was not reachable.
- **PAWITAN:** developed by leave-one-out cross-validation and DLDA on the Stockholm 159 (5-year distant metastasis or death). It was validated on Uppsala (76 node-positive treated, 135 node-negative untreated) and on van 't Veer 78.
- **KORKOLA:** the 21 genes were fixed on the UCSF cohort. Testing in the van 't Veer data used PAM classification with leave-one-out against the remaining good and poor outcome tumors (see Section 5).
- **GLINSKY:** the 11 genes were selected on prostate outcome. The breast cut-off was chosen by lowest log-rank P in the van 't Veer 78 training patients (see Section 5).
- **WELM (C):** the three genes (MSP/MST1, MT-SP1/ST14, MST1R/RON) were chosen from pathway biology.
  - In NKI, MSP alone, each gene alone and all two-gene pairs were tested against outcome and were not significant, and the all-three-high rule was adopted afterward.
  - The MSP probe was picked from four variants that had been tested against outcome.
- **WONG-MITOCHON (C):** the module was found in NKI 295 expression by its co-activation with the wound signature, without outcome. After it failed to predict survival across all 295 patients, the paper reported a top-50 versus bottom-50 contrast.

## 4. Borderline cases (the "conservative" sensitivity set)

Each line is a case where a stricter reading would raise the code: blank → C/D, V → C/D, or C → D. The cases are grouped by how strong the argument for raising them is.

**Suggested sensitivity sets:**
- **Main:** use the CSV codes.
- **Conservative:** also recode the Tier 1 and Inherited lines to their stricter code.
- **Most conservative:** also recode Tier 2 and Tier 3.

### Tier 1: outcome analysis in the cohort, or the cohort's own patients, helped make the signature

The paper reports an outcome analysis of the cohort before the final scoring or definition was fixed, or the cohort's patients were the development set (with non-outcome labels) and the main survival claim was tested in-sample.

| Signature | Cohort | Current | Stricter | Reason |
|---|---|---|---|---|
| BEN-PORATH-EXP1 | NKI | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-EXP1 | Loi_GSE6532 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-EXP1 | GSE4922 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-EXP1 | GSE7390 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-PRC2 | NKI | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-PRC2 | Loi_GSE6532 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-PRC2 | GSE4922 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| BEN-PORATH-PRC2 | GSE7390 | V | C | Cohort was in the 6-study compendium in which grade/ER-based (non-outcome) enrichment patterns preceded the choice of ES exp1 (up) plus PRC2 targets (down) as the combined 'ES signature' used for the survival test; the outcome-stratified enrichment (Fig. 3c) showed 'only weak associations' before the Kaplan-Meier definition was applied. The gene sets themselves were fixed from ES-cell literature. |
| CARTER | NKI | V | C | The van 't Veer breast expression data (NKI patients; 61 of the 78 training patients are in NKI 295) were one of the three datasets used to rank genes for CIN70 (non-outcome tFA correlation), and outcome was then tested in the same patients and in NKI 295. |
| CHI | NKI | V | C | High/low groups were defined by unsupervised clustering within NKI 295 before the outcome test, and the score threshold of 0 was examined against NKI survival with spline Cox models and endorsed (not changed). A strict reading treats the threshold check as outcome-shaped (C; the most conservative reading could call it a cut-off evaluated on outcome). |
| IVSHINA | GSE4922 | V | C | Uppsala patients were the development set (grade-supervised gene selection) and the main survival claim was tested in-sample. |
| IVSHINA | Loi_GSE6532 | V | C | Shared Uppsala patients: development set with grade labels plus in-sample survival test. |
| META-PCNA | NKI | V | C | NKI was Venet's discovery cohort: the observation that random and published signatures are prognostic in NKI overall survival led the authors to build meta-PCNA, and the R2 = 0.9 relation was established in NKI. The 131 genes and the 1% cut-off were fixed from normal tissues without NKI outcome. |
| MILLER | Loi_GSE6532 | V | C | Shared Uppsala patients were the development set (expression with p53-status labels) and the headline survival claim was tested in-sample in the same patients. |
| MILLER | GSE4922 | V | C | Same as Loi: Uppsala development patients, in-sample disease-specific survival test. |
| SOTIRIOU-GGI | Loi_GSE6532 | V | C | Loi contains the KJX64 training patients whose expression and grade labels selected the GGI genes (no outcome), alongside the KJ125 outcome-validation patients; a reading that counts any development use of a cohort's patients would code Loi as C (or D). |
| WELM | NKI | C | D | The three-gene AND rule was adopted after single genes and all pairs were tested against NKI outcome and failed, and the MSP oligonucleotide used for the main analysis was picked among four outcome-tested variants: effectively a model/variant choice made on NKI outcome. The gene set itself was chosen from pathway biology. |
| WONG-MITOCHON | NKI | C | D | The extreme-group contrast (top 50 vs bottom 50) used for the outcome claim was adopted after the full-cohort survival test in NKI was null, i.e., a cut-off/scoring choice made with NKI outcome in view; the gene list itself was fixed without outcome (expression + wound-signature labels in NKI). |
| WONG-PROTEAS | NKI | V | C | NKI is the discovery cohort: module membership was defined from NKI 295 expression and selected for coactivation with the wound signature (a classifier whose prognostic value in NKI motivated the search), then outcome was tested in-sample. |
| YU | GSE1456 | V | C | Oncomine molecular-concept enrichment of the PRC2 gene set in the Pawitan (Stockholm) 'downregulated in dead' survival concept (Table S1, Fig. S2) is an outcome analysis of this cohort that, together with the Wang concept, prompted the breast-cancer evaluation and the breast-based refinement of the signature. |
| YU | NKI | V | C | Oncomine concept 'Under-expressed in ER+ Breast Carcinoma with Relapse (vandeVijver)' was linked to the 87-gene signature in the concept map (Fig. S2) that motivated the breast analysis; the main-text P value cited is for the Wang concept. |

### Inherited provenance (ABBA only)

ABBA is a meta-signature of genes that recur in 42 published lists. Several of those lists were selected on outcome in our cohorts. The ABBA paper itself analyzed only NKI, so under the original-paper rule these are blank or V. A gene-level contamination analysis may still want to treat them as D or C.

| Signature | Cohort | Current | Stricter | Reason |
|---|---|---|---|---|
| ABBA | NKI | V | D | One of the 42 constituent gene lists is the van 't Veer 2002 Nature prognosis gene table (SuppTable4), selected by outcome correlation in the 78 NKI training patients (61 in NKI 295); the Teschendorff 2006/2007 lists (combined analyses of major public datasets) may also have used NKI outcome. Gene votes therefore partly reflect NKI outcome. Also the cluster 1 vs clusters 2/3 grouping used in the Cox model was chosen after seeing NKI survival (stricter C at minimum). |
| ABBA | GSE1456 | blank | D | Constituent list Pawitan 2005 (64-gene, PMID 16280042) was trained on Stockholm outcome; the ABBA paper itself never analyzed GSE1456. |
| ABBA | Loi_GSE6532 | blank | D | Constituent list Loi 2008 (PMID 18498629) is an outcome-based predictor trained on 255 tamoxifen-treated ER+ samples (Loi's cohorts); Sotiriou 2006 GGI list was grade-trained on GSE2990 samples. ABBA itself never analyzed Loi. |
| ABBA | GSE4922 | blank | C | Constituent list Miller 2005 p53 signature was defined on Uppsala expression with mutation labels (outcome then tested there); Pawitan 2005 was validated on Uppsala. ABBA itself never analyzed GSE4922. |

The constituent lists also include Sotiriou 2006 GGI (grade-trained on GSE2990 patients → Loi), Wang 2005 (GSE2034, not ours) and Paik 2004 (NSABP, not ours). Whether the Teschendorff 2006/2007 constituent lists used NKI outcome was not checked.

### Tier 2: clustering or filtering inside the cohort before the outcome test

Here the cohort's expression only was used to form groups or filter genes, and then outcome was tested in the same patients. No outcome is shown to have been used.

| Signature | Cohort | Current | Stricter | Reason |
|---|---|---|---|---|
| BUESS | NKI | V | C | High/low interferon groups were defined by unsupervised clustering within NKI 295 (expression only) immediately before the outcome test. |
| CHANG | NKI | V | C | Activated/quiescent classes in the van 't Veer 78 patients were defined by unsupervised clustering of CSR genes within the same patients (expression-only, first dendrogram bifurcation) immediately before the outcome test; no outcome was used. |
| CRAWFORD | GSE1456 | V | C | Per-cohort variance filter and unsupervised clustering (expression only) defined the groups before the survival test; Dataset S5 also reports per-cohort survival associations of individual signature probes (post hoc, not used to define the list). |
| CRAWFORD | GSE4922 | V | C | Same in-cohort filtering and clustering (GSE4922 and GSE3494). |
| CRAWFORD | Loi_GSE6532 | V | C | Same, via shared Uppsala patients. |
| CRAWFORD | NKI | V | C | Same in-cohort filtering and clustering on the Rosetta data. |
| HALLSTROM | GSE1456 | V | C | Risk groups were defined by in-cohort unsupervised clustering (expression only) before the survival comparison; no outcome used to define groups. |
| HALLSTROM | GSE4922 | V | C | Same in-cohort clustering before the survival comparison. |
| HALLSTROM | Loi_GSE6532 | V | C | Same, via the Uppsala patients shared with GSE4922. |
| IVSHINA | NKI | V | C | G2a/G2b groups in van 't Veer and van de Vijver data were defined by unsupervised clustering within the cohort before the outcome test. |
| MILLER | NKI | V | C | Risk groups in the van 't Veer 97 were formed by unsupervised clustering within the cohort (expression-only) before the outcome test, and both a 2-cluster and a 4-cluster split were reported. |
| SHIPITSIN | NKI | V | C | An expression-variability filter on the NKI node-negative untreated subset (Data Set 3) was applied during gene selection (genes had to pass in both Data Sets 1 and 3) before the NKI outcome test; no NKI outcome was used. |
| SORLIE | NKI | V | C | Weak: the van 't Veer tumors were hierarchically clustered and subtype-labelled within the cohort (expression only) before the survival test, and the basal gene cluster was also displayed from the van 't Veer data (Figs. 2F, 3B). Subtype labels came from fixed Norway/Stanford centroids, so the standard reading is V. |
| VALASTYAN | NKI | V | C | The positivity rule (multiple of six genes in the top/bottom 15%) was set within this cohort with no stated rationale; no evidence it was tuned on outcome, but it was defined in the only test cohort. |
| WANG-ALK5T204D | NKI | V | C | In-cohort expression filter (90 of 271 genes) and in-cohort hierarchical clustering (expression only) defined the groups tested for outcome. |
| WEST | NKI | V | C | Tumor groups A and B were chosen post hoc from an unsupervised dendrogram of the same NKI patients and each compared with all remaining tumors; the paper says they were picked for their DTF/SFT expression contrast, but branch selection could have been informed by outcome. |
| WONG-ESC | NKI | V | C | The two-class split used for the outcome test was defined by hierarchical clustering of ESC-module expression within the same 295 NKI tumors (expression-only, cohort-specific stratification) before the in-sample survival test. |

### Tier 3: speculative, or depends on something not documented

| Signature | Cohort | Current | Stricter | Reason |
|---|---|---|---|---|
| PAIK | NKI | blank | C | Weak: the 250-gene candidate pool was seeded partly from published microarray prognosis studies; if van 't Veer 2002 (NKI patients) was one of them, its outcome analysis shaped the candidate pool. Not verifiable without the NEJM full text (403); SCUBE2 is shared with the 70-gene list. |
| SOTIRIOU-93 | NKI | blank | D | Only if Venet's SOTIRIOU-93 list were the 93-probe-element overlap with the van 't Veer 231 genes (genes originally selected on van 't Veer 78 NKI outcomes). The 343-symbol size argues strongly against this; listed only for completeness. |
| SOTIRIOU-93 | Loi_GSE6532 | blank | D | Possible undocumented patient overlap: the 99 development patients were John Radcliffe Hospital (Oxford) cases from January 1993 to December 1994, and Loi_GSE6532 contains Oxford untreated (OXFU) and Oxford tamoxifen (OXFT) JRH samples, including GSE2990. No paper states an overlap and it is not among the established overlaps; if a sample-level check finds shared patients, Loi would be D (the 485-element list was selected on relapse-free survival in these 99 patients). |

Note on Tier 3: Sotiriou 2003 (99 John Radcliffe Hospital, Oxford, patients operated January 1993 to December 1994) was also a validation set for KORKOLA and SAAL. If those patients overlap the Oxford samples in Loi (OXFU, OXFT, GSE2990), SOTIRIOU-93 would be D for Loi, and KORKOLA and SAAL would be V for Loi. No paper documents such an overlap, and Sotiriou 2006 treats the 2003 "NCI" dataset as independent.

## 5. The opposite direction: D codes resting only on cut-offs or refitting inside the cohort

These D codes come from fitting a cut-off or a classifier on the cohort's outcome. The gene set itself was fixed elsewhere, so a gene-set-only reading would make them V. Venet-style scoring (the median of PC1, with no published cut-off) does not carry these cut-offs into the score.

| Signature | Cohort | Current | Lenient | Reason |
|---|---|---|---|---|
| GLINSKY | NKI | D | V | The 11 genes were selected on prostate cancer outcome (not our cohorts); only the breast cut-off was chosen by lowest log-rank P in the van 't Veer 78 training patients. |
| HUA | NKI | D | V | Gene set fixed in MCF-7 cells; only the two cut-points of the RA score were optimized on outcome with X-tile (application to GSE4922 inferred from Figure S14). |
| HUA | GSE4922 | D | V | Gene set fixed in MCF-7 cells; only the two cut-points of the RA score were optimized on outcome with X-tile (application to GSE4922 inferred from Figure S14). |
| HUA | Loi_GSE6532 | D | V | Gene set fixed in MCF-7 cells; only the two cut-points of the RA score were optimized on outcome with X-tile (application to GSE4922 inferred from Figure S14). |
| KORKOLA | NKI | D | V | The 21 genes were fixed on the UCSF cohort; only the class centroids (good/poor means, leave-one-out) were refit on van 't Veer outcome labels. A reading that counts only gene selection would make this V. |

Also lenient-sensitive: HE Loi_GSE6532 = D depends on the overlap arithmetic described in Section 3.

## 6. PMID corrections and list identity

| Signature | Venet PMID | PMID used | Finding |
|---|---|---|---|
| VANTVEER | 12490681 | **11823860** | 12490681 is van de Vijver et al. 2002 NEJM, which validated the already-fixed 70-gene profile in NKI 295. The 70 genes were selected in van 't Veer et al. 2002 Nature (PMID 11823860), which is the defining paper. Venet's 60 symbols match the 70-gene profile. |
| PAWITAN | 16846532 | **16280042** | 16846532 is Calza et al. 2006 Breast Cancer Res (412 Stockholm + Uppsala patients). That paper validates the ~500-gene Norway/Stanford intrinsic set and contains no 46-gene list. Venet's 46 symbols (CDKN1C, HOXA5, IGF1, PDGFD, SLIT2, PTN, PRC1, DTL, BUB1B, CDC2, CDC20, ...) match the 64-gene list in Additional file 1 of Pawitan et al. 2005 Breast Cancer Res 7:R953 (PMID 16280042), which is the defining paper. |
| META-PCNA | none | **22028643** | Defined in Venet et al. 2011 PLoS Comput Biol itself (PMC3197658). Development used the normal-tissue compendium of Ge et al. 2005 (GDS1096). Venet then tested it against outcome in NKI 295 and Loi GSE6532, so both are V. |
| RHODES | 15184677 | 15184677 | PMID correct. Venet's 67 genes are the "meta-signature of neoplastic transformation" (cancer versus normal), not the 69-gene undifferentiated-cancer meta-signature from the same paper. |
| SOTIRIOU-93 | 12917485 | 12917485 | PMID correct, but the label "93" is ambiguous. The paper's "93 probe elements" are the overlap with van 't Veer's 231 genes (56 unique genes) and cannot produce Venet's 343 symbols. The only list of matching size and content is the 485 probe elements associated with relapse-free survival at P < 0.05 in the 99 Oxford patients (Supplementary Table 9, which could not be downloaded). |
| CHANG | 14737219 | 14737219 | PMID correct. Chang 2004 PLoS Biol defines the CSR gene set. Chang 2005 PNAS (PMID 15701700) later applied it to NKI 295 and fitted a cut-off there, but it is a later paper and does not count. |
| WANG-ALK5T204D | 18625725 | 18625725 | PMID correct. The expression signature was derived in BT474 cells, not MCF10A/HER2. |

The other 38 PMIDs were checked against the papers and are correct.

## 7. Uncertainty and access limits

Confidence is high for 35 rows, medium for 9 and low for 1 (KOK).

- **KOK:** full text not reachable (Wiley 403, no open-access copy). The gene set comes from mammospheres of ER-positive cell lines (abstract). NKI 295 use (122 ER-positive untreated and 69 ER-negative patients) and a second NKI-AVL tamoxifen series of 129 patients (not one of ours) come from search-engine excerpts. D or C is unlikely because the abstract says the set was defined from cell lines only.
- **Paywalled main text:** NEJM (PAIK, LIU, and the van de Vijver companion of VANTVEER), Lancet (WANG-76), Nature Genetics (CARTER, RAMASWAMY) and Cancer Cell (SHIPITSIN) main texts returned HTTP 403. Coding uses abstracts plus the papers' own supplementary methods, which were reachable for CARTER, RAMASWAMY, SHIPITSIN and VANTVEER. For PAIK and WANG-76 the design is unambiguous: NSABP and Rotterdam GSE2034, none of our cohorts.

Per-signature notes:

- **ABBA** (high): Under the user's rule only the original paper's own use counts, so constituent-signature training cohorts are listed as borderline rather than coded. Whether the Teschendorff 2006/2007 constituent lists used NKI outcome was not checked in their full texts. Venet list has 111 symbols vs the paper's 117 genes (alias mapping).
- **BEN-PORATH-EXP1** (high): The survival test used the combined ES signature (ES exp1 over AND PRC2 targets under), not ES exp1 alone; coded V for each component set. Pawitan n=159 identifies the Stockholm (GSE1456) series.
- **BEN-PORATH-PRC2** (high): The survival test used the combined ES signature (ES exp1 over AND PRC2 targets under), not PRC2 targets alone; coded V for each component set.
- **BUESS** (high): Venet's BUESS list has 30 symbols versus 42 genes in the paper (symbol mapping); does not affect coding.
- **CARTER** (medium): Main text paywalled; coding uses abstract and the paper's Supplementary Methods, Table 1 and Figures 2-4. Identity of the three ranking datasets is read from asterisk positions in Suppl. Fig. 3. Whether any breast analysis (e.g., grade 1/2 stratification in the main text) used a different dataset could not be checked, but all named breast datasets are listed above.
- **CHANG** (high): Venet's 355-symbol CHANG list cannot be tied by gene count alone to the 459-UniGene 2004 list versus the 380-gene NKI-matched subset of Chang 2005 PNAS; both are the same 2004 CSR gene set, so the coding is unaffected. Chang 2005 PNAS also fit a centroid-correlation threshold on a training half of NKI 295, but that later paper does not count and Venet uses only the gene list.
- **CRAWFORD** (high): Whether the Dutch Rosetta set is the van 't Veer 2002 (78/97/117) or the NKI 295 data is ambiguous (the URL says vantveer, but overall survival, node-positive and ER data suggest NKI 295); either maps to NKI. Venet's 377 genes match the 379-gene human list.
- **GLINSKY** (high): D rests on the breast cut-off only (fitted on outcome in the van 't Veer 78, 61 of whom are in NKI 295); the 11 genes were not selected on breast outcome. A gene-set-only (lenient) reading would make NKI V. Gene selection used prostate outcome set 1, which is not one of our cohorts.
- **HALLSTROM** (high): Supplementary Figure 2 (recurrence) could not be downloaded (PMC binary files behind a challenge page); it likely shows GSE2034 and/or GSE4922 recurrence. Venet's 78-gene HALLSTROM list mixes genes from Tables 1 and 2 (PI3K-repressed and non-repressed E2F1 targets).
- **HE** (high): Loi_GSE6532 = D depends on the overlap arithmetic: He removed only the 86 GSE4922 patients that overlap GSE2990, while Loi contains ~117 Uppsala patients, so ~30 Loi patients were likely in the 163-patient gene-selection set. If the parent's sample-level overlap shows no Loi patient among He's 163, recode Loi_GSE6532 as V (GSE2990 and GSE7390-shared patients were validation only).
- **HU** (high): Supplemental Experimental Procedures (exact phenotype definitions, e.g., metastasis within 5 years, and which van 't Veer patients) could not be downloaded; the coding does not depend on them.
- **HUA** (medium): D rests only on outcome-optimized cut-points (X-tile) for the three-group P/N/U split; the 1,345-gene set itself was fixed in MCF-7 cells without patient data, so a reading that ignores cut-offs not carried into the gene list would recode NKI, GSE4922 and Loi_GSE6532 as V. The X-tile sentence is generic ('the patient sample data'); its application to the GSE4922 cohort (same P/N/U grouping in Figure S14) is inferred, not stated per cohort. Main-text Figure S14 call-out says only 'two independent breast tumor cohorts'.
- **IVSHINA** (high): Table 1 gene symbols were not machine-readable on the publisher page; Venet's 17 symbols are consistent with the 18-probe PAM classifier (plus or including the 6-probe SWS set). Survival figures shown are for the SWS classifier; the PAM classifier's survival results are described only as concordant. GSE1456 code relies on the Stockholm G2 subset (58 of 159 patients).
- **KOK** (low): Full text not read. Use of NKI-295 (122 ER+ untreated, 69 ER-) and the identity of the second ER+ cohort (NKI-AVL tamoxifen series, n=129) come from search-engine excerpts, not verbatim reading. Any overlap of the NKI-AVL tamoxifen series with the NKI 295 patients is unknown (would still be V). Abstract states the gene set was defined from cell lines only, so D/C is unlikely.
- **KORKOLA** (medium): Which van 't Veer subset (78, 97 or 117) was used is not stated. Sotiriou 2003 Oxford patients may overlap Loi OXFU/GSE2990 Oxford samples (same hospital and group), but this is not an established overlap.
- **LIU** (medium): NEJM full text and supplement returned HTTP 403 (curl and WebFetch). NKI 295 use and the fixed cut-off come from the abstract plus a secondary summary. Any additional breast cohorts used in the supplement (e.g., Wang/GSE2034) could not be checked; none of our other cohorts is expected.
- **MA** (high): The match of Venet's 30 symbols to Table 1 (39 clones) is by content (IDH2, KNSL1/KIF11, ANKT/NUSAP1, UBE2C, ANLN, CENPA, PSMD12, TOPK/PBK, RRM2, TOP2A, PLK, NEK2, CDKN3, RACGAP1, BIRC5, CKS2 and others in Table 1 are in Venet's list).
- **MILLER** (high): Venet's 18 symbols are the annotated subset of the 32-gene DLDA classifier (Fig. 1A); supplementary tables (PNAS site) were not read.
- **MORI** (high): The ROC-chosen 0.5 cut-off used outcome in GSE2603/GSE5327, which are not among our 13 cohorts.
- **PAIK** (medium): NEJM full text and supplementary appendix returned HTTP 403; coding rests on the abstract and the first author's 2007 review. No breast microarray cohort among our 13 was used for development or validation.
- **PAWITAN** (high): If the Venet PMID (Calza 2006) were taken at face value, that paper used Stockholm 159 + Uppsala 253 outcome to evaluate intrinsic subtypes (V-type use of GSE1456, GSE4922, Loi) but it does not define Venet's 46 genes; coding follows Pawitan 2005.
- **RAMASWAMY** (high): Main article text (Nature Genetics) is paywalled; coding uses the paper's own Supplementary Information, which describes the breast analysis in full. outcome_guided=yes refers only to the lung reduction step (outcome-associated clusters), not to any of our cohorts.
- **SAAL** (high): Sotiriou 2003 Oxford patients may overlap Loi OXFU/GSE2990 Oxford samples; not an established overlap, not coded.
- **SHIPITSIN** (medium): Main text (Cancer Cell) blocked (HTTP 403); coding uses the paper's Supplemental Data PDF (verbatim methods). Main-text Table 1 (gene lists of signatures A/B) could not be read, so the match of Venet's 56-gene list to signature A is inferred from gene content; if Venet instead used a SAGE-only CD44+ list without the outcome filter, Loi would become V.
- **SORLIE** (medium): Which figure Venet took the 15 basal genes from is inferred from labelled genes (CHI3L2 appears only in the Norway/Stanford Fig. 1E labels; GABRP, SLPI, ZDHHC5, C1orf198 are unlabelled in the figures), so the gene list most likely comes from the Norway/Stanford basal cluster (full list in supporting Fig. 6, not read). The outcome test in van 't Veer was of subtype membership, not of the 15-gene cluster score itself.
- **SOTIRIOU-93** (medium): Supplementary Table 9 (485 probe list) could not be downloaded (PNAS 403, PMC download challenge), so the match of Venet's 343 symbols to the 485 list rests on list size and content. Overlap of JRH Oxford 1993-1994 patients with Oxford samples in GSE2990/Loi is not documented in the literature read.
- **SOTIRIOU-GGI** (high): Venet's 90 symbols correspond to the 97-gene/128-probe GGI list. Exact Oxford/Uppsala split of KJX64 and KJ125 rows in Table 1 was only partly legible on the publisher page.
- **TAUBE** (high): The NKI dataset is cited as reference 51 (van de Vijver M, Oncologist 2005), not the 2002 NEJM paper; it is presumably the NKI 295 set but the patient count is not stated. The outcome test was negative but still counts as V (outcome tested in the original paper).
- **TAVAZOIE** (high): None for our cohorts. Supplementary Fig. 15 (per-cohort curves) not read, but the main text states the cohorts are MSK and EMC only.
- **VANTVEER** (high): None for the coding. The NEJM text describing the 61 overlapping patients could not be retrieved (HTTP 403); the overlap is documented by the user's established overlap.
- **WANG-76** (medium): Lancet full text returned HTTP 403; abstract only. The design (single Rotterdam series, GSE2034) is unambiguous, so no cohort among our 13 is involved.
- **WANG-ALK5T204D** (high): Venet's list has 239 symbols vs 271 genes in the paper's Table S1 (symbol mapping); not checked gene by gene.
- **WHITFIELD** (high): Venet's 587 symbols vs 874 UniGene clusters reflects symbol mapping; does not affect coding.
- **WONG-ESC** (high): The Segal et al. 2004 cancer compendium used for the tumor-vs-normal module map (Figure 2A) may include van 't Veer 2002 arrays with annotations; the paper does not list them. This would not change the NKI code (V).
- **WONG-MITOCHON** (high): The paper does not state explicitly why 50/50 was chosen (the module-map figure also highlights the samples with significant coregulation), so whether the cut-off was outcome-driven is inferred from the sentence order. Venet's list has 217 genes vs 218 in the paper (symbol mapping).
- **WONG-PROTEAS** (high): Venet's list has 46 genes vs 50 in the paper (symbol mapping).
- **YU** (high): The Sotiriou et al. breast data set (n=94) cited as Sotiriou 2006 JNCI is assumed to be a subset of GSE2990 patients (all in Loi_GSE6532); exact patient subset not stated. The main-text P = 9.7e-8 matches the Wang ER+ relapse concept, not the Pawitan concept.
