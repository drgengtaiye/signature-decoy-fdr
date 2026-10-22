# Signature x cohort provenance: notes

Companion to `signature_provenance.csv` (52 signatures x 12 ORCESTRA cohorts). Each original publication's Methods, Results and figure legends were read, and the supplement where it was reachable. Sources: PMC JATS full text, publisher HTML, or author manuscripts. Every row's `evidence` cell gives the quote or section and the URL.

## Coding rule
- **D**: the cohort's ICB outcome data were used to pick genes or to fit weights, cutoffs or a model.
- **C**: the gene set was fixed without this cohort's outcomes, but analysis of the cohort itself prompted or shaped the signature. In practice this is the paper's own discovery cohort.
- **V**: a finished signature was only tested in the cohort.
- **Blank**: the cohort was not used.

Following the coding rule, a signature fixed on non-ICB data (TCGA, mouse, in vitro, chemotherapy or other non-ICB cohorts) is coded V in every ICB cohort where it was only tested. This holds even when those cohorts were the first, or the only, evidence for the ICB claim. Those cases are listed under "Borderline" below.

## All D and C entries (leave these out in the rerun)

| Signature | Cohort | Code | Evidence |
|---|---|---|---|
| CRMA_Shukla (29656892) | Van_Allen | D | Results: "A differential gene expression analysis comparing NB (n=22) and CB (n=13) groups identified 457 genes upregulated in NB"; the Van Allen ipilimumab RNA-seq is called the "discovery cohort". Validation was CheckMate 064 (not among the 12). Hugo = V. |
| MPS_Guijarro (32284588) | Van_Allen | D | Methods "MPS development": genes were added stepwise "to classify patient samples in the Van Allen dataset… We chose the largest gene list that had a P value <0.05 as derived from Fisher's exact test". The 33rd-percentile cutoff was also set in Van Allen. Hugo and Riaz = V. |
| ImmuneCells_Xiong (33033253) | Hugo | D | "we developed an ICT responsiveness signature based on the scRNA-seq dataset and a bulk gene expression dataset - GSE78220 … in the discovery dataset". The peer-review file says GSE78220 "was used to create the signature". Riaz, Gide and Liu = V ("Neither of these datasets were used to develop the ImmuneCells.Sig"); caveat below. |
| IMS_Cui (33542239) | Riaz, Hugo, Gide | D | "three external datasets … which we used as our discovery cohort (n = 120)" (Riaz n=51, Hugo n=28, Gide n=41, anti-PD-1 monotherapy only). Methods "Identification of the IMS genes": one-sided t-test of PD vs PR/CR after normalising to the IFN-γ score, run per dataset and combined. Puch (the paper's own PUCH cohort), Liu, Van_Allen and Kim = V. |
| TGFB_Mariathasan (29443960) | Mariathasan | C | The 19 Pan-F-TBRS genes come from TGF-β-stimulated fibroblasts in vitro. The TGF-β hypothesis came from IMvigor210 outcome analyses ("Lack of response was associated with a signature of TGF-β signalling in fibroblasts"; TGFB1 and TGFBR2 were the top non-response genes). The Pan-F-TBRS was then shown to be "associated with non-response in excluded tumours" in IMvigor210. This sits close to D. |
| MHC-II_Liu (31792460) | Liu | C | A fixed list of 13 MHC-II HLA genes. Results: "expression of all 13 MHC-II-associated HLA genes was higher in responders (P = 0.0002)". It was also the top feature in Liu's outcome-selected logistic model, which is not the SignatureSets score. Riaz = V (the paper's validation cohort; not significant there). |
| MHC-I_Liu (31792460) | Liu | C | A fixed list (HLA-A/B/C/E/F, TAP1, TAP2, B2M). "MHC-I antigen-presentation genes all trended toward having higher expression in responders (… P = 0.02)". Medium confidence, because the evidence is only a trend. |

No D or C entries exist for Braun, Jung, Kim, Limagne1, Miao1 or Puch.

## Borderline entries (coded V or blank; a stricter reading could make them C or D)
- **C-ECM_Chakravarthy, Hugo and Riaz (V; first reading coded C).** The 58 genes come from TCGA tumour vs normal. Hugo, Riaz and Snyder 2017 were the only ICB cohorts and were pooled. The paper's separate 19-gene "practical signature" was selected on pooled Hugo+Riaz+Snyder outcomes, which would be D if that subset were used. SignatureSets uses all 58 genes.
- **STAT1_Care, Mariathasan (V; first reading C).** The genes come from Care 2015 (DLBCL, no ICB). Zemek 2019 tested them in IMvigor210 first. Its PFS classifier (CABS) built responder/progressor centroids from IMvigor210's own labels. Riaz = V.
- **ImmunoScore_Hao, Van_Allen, Hugo and Riaz (V; first reading C).** The weights were fitted on ovarian-cancer survival without ICB, then the score was meta-analysed across four ICB datasets.
- **PTEN-MITF_Cabrita, Gide and Riaz (V; first reading C).** The centroid and cutoff were fixed in the Lund cohort (no ICB), then Gide and Riaz were pooled for one survival analysis.
- **KDM5A_Wang, Liu and Riaz (V; first reading C).** The 8 genes come from TCGA plus the Leeds cohort (no ICB). Cross-validated models were fitted within each cohort, and the supplement could not be accessed.
- **TLS_Cabrita, Van_Allen, Gide and Riaz (V).** The 9 genes come from the Lund cohort (no ICB), selected by an IHC grouping. These three cohorts and the authors' Danish cohort carry the whole ICB claim.
- **APM_Wang, Hugo and Van_Allen (V).** The genes come from a literature review. These are the paper's only patient-level ICB tests, alongside Snyder 2017.
- **CYT_Davoli, Van_Allen (blank; first reading C, low).** Van Allen was used only for aneuploidy (SCNA) level and a separate 23-gene set. The 7-gene CYT score is never tested against ICB outcome.
- **ImmuneCells_Xiong, Riaz, Gide and Liu (V).** Only the 108-gene list was fixed in advance. The authors' code refits the nearest-centroid classifier inside each "validation" cohort (in-sample AUCs), and regimen-specific z-score cutoffs (0.19 and 0.1) were set in the Gide subsets. These would be D if the published classifier or cutoffs were used rather than an unweighted GSVA score.
- **Cutoffs only, which do not matter for continuous scores:**
- IMPRES: the suggested cutoff of 8 came from pooled ICB data that include Hugo, Riaz and Van Allen.
- IMS_Cui: Youden cutoffs were chosen within Puch, Liu, Van Allen and Kim.
- TIDE/CTL: cutoff choices were compared on Hugo and Van Allen as a robustness check.
- **Expression-only uses (V, no outcome test):**
- CD8_Sade-Feldman: Van_Allen and Riaz, correlation with CD3 only.
- ImmuneScore_Roh: Van_Allen, agreement with ESTIMATE only.

**Conservative sensitivity set:** add the first eight bullets to the D/C list above.

## Cohort identities and overlaps
- **ICB_Puch** is Cui C et al., *NPJ Genom Med* 2021;6:7 (PMID 33542239). It is the "newly generated RNA-seq dataset from 55 … melanoma patients treated with anti-PD-1 monotherapy at Peking University Cancer Hospital (PUCH)", with raw data at GSA-Human HRA000524.
- Confirmed from the ORCESTRA API (https://orcestra.ca/api/clinical_icb/available): the ICB_Puch publication is PMID 33542239, Zenodo 10.5281/zenodo.7332226.
- Confirmed from the local `ICB_Puch.rds`: 55 samples, `data_set = PUCH`, melanoma, PD-1.
- It is the same paper that defines IMS_Cui. There, Puch was the internal validation cohort (V), not a discovery cohort.
- **ICB_Limagne1** is Limagne E et al., *Cancer Cell* 2022 (PMID 35051357). The local rds holds 70 anti-PD-1 NSCLC samples, GSM5718805 onward, from GSE190266 "[France4]" (SuperSeries GSE190270), submitted by INSERM UMR 1231 in Dijon in December 2021.
- No signature uses it.
- VIGex used "Fumet.1/Fumet.2", which are Dijon NSCLC data from Fumet 2018 *Br J Cancer*. Whether those patients overlap with Limagne1 was not checked. If they do, VIGex x Limagne1 would be V; that changes nothing for the leave-out.
- **Braun overlaps Miao1.** Braun 2020 Methods: "For the CheckMate 009 cohort, the previously published TPM matrix was used" (ref 19 = Miao 2018). No signature is D or C on either cohort.
- **Liu and Van Allen** share the dbGaP study accession phs000452 but are different cohorts (Liu: 144 anti-PD-1 patients; Van Allen: ipilimumab). The papers using them were disambiguated by citation, sample size and drug.
- **Riaz** is CA209-038. Some papers used only the ipilimumab-naive subset (MPS_Guijarro, for example).

## What the curation and PredictIO document
- Neither documents development datasets.
- **SignatureSets** `data-raw/signature_information_io.csv` has these columns only: signature, association, RNA.type, method, cancer.type, score.function, PMID, curated, information, immunotherapy, gitHub.status.
- **PredictIO** (`bhklab/PredictIO`, `data/signatures/signature_INFO.txt`) has only Signature, method, association and PMID. The Bareche 2022 *Ann Oncol* full text returned 403. Its abstract says only that signatures were "computed with respect to the original publication".
- **Attribution differences.** PredictIO cites Macrophage as PMID 33520406 (Budczies), PDL1 as 30127393 (TIDE) and T_cell_inflamed as 30804515 (Ott 2019). The gene lists match SignatureSets, and the matrix follows the SignatureSets PMIDs.

## Signatures not fully verified (medium or low confidence)
- **CIN25/CIN70_Carter (2006).** Abstract only, because the paper is paywalled. It predates every cohort and uses only non-ICB microarrays, so the blanks are safe.
- **CYT_Davoli.** Main text read; supplementary tables S5 and S6 were not reachable.
- **STAT1_Care.** Read from the Zemek 2019 author manuscript; science.org returned 403.
- **TLS_Cabrita.** Paywalled, and no author manuscript was found. The coding rests on the abstract, the Fig. 3 and Extended Data legends, Extended Data Tables 4–5 and the Reporting Summary.
- **KDM5A_Wang.** Main text read from an archived copy of the lab's PDF; the supplementary methods were not reachable.
- **WNT-beta-catenin_Takeuchi.** Abstract, the NCC press release and the Crossref reference list only (none of the 12 cohort papers is cited). How the 17 genes were derived is not verified.
- **GeparSixto_Karn.** Abstract only. It is a "predefined" GeparSixto chemotherapy-trial immune profile tested in GeparNuevo.
- **Teff-IFNG_Fehrenbacher.** Abstract only. POPLAR was published on 10 Mar 2016, before every cohort except Van Allen, which it does not cite.
- **Not checked at all:**
- Thompson 2020 Supplementary Table 1 (Inflammatory_Thompson genes)
- Dangaj Key Resources table
- Xiong's exact gene-selection steps within Hugo (the code starts from the final list)

## SignatureSets problems found along the way (outside the matrix)
- **IRG_Yang.** The ±1 weights follow tumour-vs-normal direction, not the Cox coefficients. Five genes have conflicting signs, so the "sensitive" label may be inverted.
- **Direction labels.** Myeloid-DC_Helmink is labelled "resistance", but its median is higher in responders. CYT_Davoli is labelled "resistance", which looks like the wrong direction.
- **Tcell-inflamed_Ayers.** SignatureSets weights every gene 1; the paper used regression weights.
- **PDL1_Nishino.** The 3-gene set (CD274, PDCD1, PDCD1LG2) is not defined in either cited paper.
- **Inaccurate descriptions.**
- Blood_Friedlander is whole-blood qPCR, not RNA-seq.
- NonResponse/Response_Chen was never tested in a bulk human ICB cohort.
- proliferation_Pabla's original claim is U-shaped, not monotone.

