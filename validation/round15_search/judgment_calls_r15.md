# Round 15 screening: judgment calls reported by the readers, to resolve in the central check

## c10_OV (83 screened, 32 eligible, trains_tcga_model TRUE 79)
- PMC12167356: GRIP2 coefficient typeset "1" + "1.80464867091524"; recorded 11.80464867091524 as printed -> check full text.
- PMC7720540: 5-gene model, printed formula has 4 terms and "DHRS95" (Discussion: DHRS9) -> ruled "coefficients not given in text".
- PMC8267297: cohort of the multivariate Cox fit not stated (TCGA + unnamed GEO in screening) -> counted TCGA-trained, eligible (round-13 precedent SKCM PMC9715341).
- PMC10795461: Methods sentence CLDN4 0.0839, formula (twice) -0.0839 -> formula used.
- PMC11711411: "IL1B × (−0.3618734) × IL21", "+" meant -> two terms.
- PMC7749778: print splits joined (CH25H, HSPB7, LOC158830 = CXorf65, PPM2C alias of PDP1).
- PMC9679378: formula ends with a constant (−2.427119), left out (does not change ranks).
- PMC8595115: first score is a PC (ICI) score -> ineligible; later 18-gene Cox model coefficients only in figure/supplement.
- Conventions: ICGC cohorts as "OV-AU" only when named, else "ICGC"; TCGA internal test splits left out; IMvigor210 only where survival by risk group is shown.

## c11 backfill of the 294 round-13 records (r13_backfill.csv; trains_tcga_model TRUE 272, FALSE 22; 142 TRUE with validation_used)
- 294 rows, 292 distinct PMCIDs: PMC8938586 (SKCM rank 11 and CESC rank 18) and PMC11276689 (STAD 16, BLCA 14; an ovarian paper) sit under two types -> every merge keys on (type, pmcid). Two rows appended by the reader's own script (tmp_c11_backfill/append_dups.py) because add_rows_r15.py refuses a repeated PMCID.
- Convention differs from c10_OV: checkpoint-treated cohorts always left out of validation_used here (OV reader kept IMvigor210 where survival by risk group was shown). Immaterial for part B: B2 excludes checkpoint-treated cohorts.
- TRUE judgment calls: PMC7563367 (STAD), PMC8938586 (per-type models in a pan-cancer paper); PMC10637804 (BRCA; one pooled pan-TCGA model).
- FALSE judgment calls: PMC9878554 (LUAD, subtypes only), PMC11635229 (LUAD, pathology image features), PMC9515658 (COAD, trained on GSE40967).
- validation_used unclear: PMC9941880 (SKCM, GSE65904 batch-corrected with training), PMC7880322 (BLCA, "unnamed GEO cohort", likely GSE13507), PMC8743383 (BLCA, GSE31864 as printed, probably GSE31684), PMC9798726 (BRCA, GSE158309 printed, n matches GSE20685), PMC11164503 (OV, GSE70769 and GSE102073 both named), PMC8923768 (STAD, GSE62254 probably not with TCGA coefficients), PMC9691887 (HNSC, GSE102349), PMC11911325 (LUSC, protein-staining cohort excluded); empty for PMC11574080, PMC12148105, PMC12666366, PMC10576492.
- Misprinted accessions recorded in corrected form with the misprint in notes (except PMC8743383).

## c09_SKCM_CESC (164 screened: SKCM 104 / 30 eligible / 77 TCGA model; CESC 60 / 28 / 59)
- PMC12219082: text formula lists CLEC2D twice; 9th term SOD2 from Table 1 and abstract.
- PMC8595480: formula prints "GPX" -> GPX2 from the Discussion's gene list; GSE78220 and MSKCC counted in validation_used (OS tested, not only response).
- PMC8253645: text cites Table 1 but gives GAPDH 0.238 vs Table 1 0.327653 -> Table 1 for all eight genes (round-13 OV precedent).
- PMC11319947 (eligible): printed weights 1.29, 1.12, 0.77 look like hazard ratios (DUOX1 described as protective) -> taken as printed. TO RESOLVE: HR-as-coefficient rule.
- PMC8446628 (ineligible): HR-like weights applied to median-split 1/2 coding; formula quoted in notes.
- PMC8927698, PMC7186619, PMC8517401: fitting cohort / use of GEO in gene selection ambiguous -> TCGA taken as training; GSE59455, GSE65904 recorded as validation though used in gene selection.
- PMC12481447: training assumed TCGA merged with GSE30759 (IL1B not in formula_text, explained in notes).
- PMC9091987: first model is diagnostic -> first survival model (5-gene hypoxia) used.
- PMC9288220: text coefficients are Table 3's truncated -> text used.
- PMC10667679: Table 1 HRs all > 1 though six coefficients negative -> coefficients as printed.
- trains_tcga_model TRUE for non-standard models: PMC12605100 (logistic immune-subtype classifier), PMC7298832 (slide-image model), PMC7563367 (pan-cancer with a per-type SKCM score).
- PMC9852722: validation "ICGC dataset" with no project named; may overlap TCGA.

## c07_BLCA (169 screened, 69 eligible, trains_tcga_model TRUE 160)
- PMC12008088 (eligible): final 20-gene score from a 10-algorithm search (RSF + elastic net) but printed as a linear formula; last gene printed "HFG" -> HGF.
- PMC9493074: 5-gene hexosamine Cox score fitted in GSE13507 + GSE32894 only -> "not trained on TCGA".
- PMC10853254: training TCGA-BLCA merged with GSE13507 (655) -> no external validation recorded.
- Training TCGA merged with GEO (eligible): PMC9343957, PMC9575963, PMC12221931, PMC10694890, PMC11593262, PMC9577248.
- PMC7967990: "prepronociceptin" -> PNOC. PMC9210750: "PD-1" kept as printed (paper's table: PDCD1) -> map at A5.
- PMC10412633: text "fourteen" genes, formula 15 terms -> 15 recorded. PMC11582591: text seven genes, formula 4 -> 4 recorded.
- PMC6526186: external validation is of the composite index (score + age + stage), not the score alone.
- PMC7350589: snoRNA model -> "not gene expression: snoRNA model".
- PMC7563367: pan-cancer, ineligible; trains_tcga_model TRUE (BLCA-specific score).
- PMC9814967: stomach paper whose abstract says "BLCA" in error -> type mismatch.

## c04_COAD_b (150 screened, 47 eligible, trains_tcga_model TRUE 139)
- PMC9513614: DLAT coefficient printed "08,399" (not a number) -> recorded -0.8399, eligible. TO RESOLVE: keep or exclude (a printed value that is not a number).
- PMC8034457: "‐VIM × ‐0.70744" (minus operator and negative coefficient) -> recorded negative for VIM and FLNA. TO RESOLVE: arithmetic reading would be positive.
- PMC12604944: "HHMR" -> HMMR from the abstract (round-13 ACSF -> ACSF2 precedent).
- PMC9531163: comma decimals "-0.397,447" -> -0.397447, matches the cited Table 4.
- PMC10234596: missing operator before HOXC6 -> "+".
- PMC11830669: PSMD12 "(-) 0.5928…" -> negative; TCGA + GEO merged for gene screening, Methods name TCGA as training.
- PMC6926093: 2-gene LASSO in one consensus-cluster subgroup of TCGA-COAD -> subset of the type; "risk signal values" read as coefficients.
- PMC9731216: JMJD7-PLA2G4B readthrough kept as protein-coding.
- PMC7993742: IGHG1 -> ineligible by analogy with IGKC (amendment 2 rule 1).
- PMC8888821: genes picked by RSF, final model linear multivariate Cox -> eligible.
- PMC8847155 and PMC7969881: near-duplicate papers (same five glycolysis genes, near-identical coefficients) -> both eligible. NOTE for estimates: duplicates.
- validation_used: TCGA-READ listed for PMC10488217; GSE17538 for PMC8714135 (survival curve in drug section); AC-ICAM short name for PMC12267444; PMC12695801 (TCGA-READ only) -> trains_tcga_model FALSE.

## c01_LUAD_a (162 screened, 71 eligible, trains_tcga_model TRUE 156)
- TCGA LUAD+LUSC trained together counted for LUAD (mirror of amendment 2 rule 2): PMC6953332 (DFS), PMC12321723.
- Several models: PMC7333141 first of four Table-4 models (m1A, 3 genes); PMC8567015 3-gene main model (title/abstract) though a 20-gene formula is printed first; PMC9243308 first of three (11-gene mRNA), +23.20828 centering constant dropped. TO RESOLVE: "first presented as the main model" when title/abstract disagree with print order.
- Not a Cox fit but used for survival: PMC8182055 (LASSO-logistic for nodal metastasis), PMC12686163 ("diagnostic" LASSO-logistic). TO RESOLVE: eligibility of non-survival-fitted scores.
- ML pipelines ending in a printed linear formula (eligible): PMC12402859 (LASSO + survival-SVM; 13 of 22 genes printed by full name, mapped to symbols; least certain), PMC11550309 (elastic net), PMC12066106 (ridge).
- PMC9730711 comma decimals "0.105,367", "0.134,121" recorded as printed -> parse as 0.105367, 0.134121.
- PMC10767103: three tiny all-negative coefficients although higher score = worse OS.
- PMC10726488: "IL-1A" -> IL1A.
- PMC12603507: Methods call GSE68465 training, Results build in TCGA -> TCGA-trained, eligible. PMC9574025: the reverse -> not trained on TCGA.
- TCGA merged with GEO then split, validation_used empty: PMC11840156, PMC10031379.

## c03_COAD_a (150 screened, 58 eligible, trains_tcga_model TRUE 138)
- PMC8687688: formula in TeX markup ({\rm{0}}.{\rm{24}}); coefficients checked by eye (flagged by the script check).
- PMC8610139: Table 1 column headed "expression"; signed 7-decimal values read as LASSO coefficients (text cites the table) -> eligible.
- PMC10052662: three formulas; combined 3-gene model (HSPA1A, CDKN2A, UCN3) taken, not the first-printed 2-gene. TO RESOLVE with the "main model" rule.
- PMC9792841: 6-gene Cox formula printed first, 4-gene LASSO signature recorded (abstract and later results use it). TO RESOLVE with the "main model" rule.
- PMC8456021: "Status of PARPBP…" read as expression; only validation an in-house immunofluorescence cohort.
- PMC10583017: JMJD7-PLA2G4B readthrough accepted.
- PMC10976740: main result is clustering; 5-gene risk score recorded.
- PMC11781745: printed signs contradict the stated high-risk direction -> as printed.
- PMC9421078: coefficients labelled "lambda = …" read as LASSO coefficients.
- PMC11191489, PMC8812245: READ-only -> type mismatch (PMC8812245 trains_tcga_model TRUE via a comparison COAD model).
- PMC11839610: pooled pan-GI training -> type mismatch; trains_tcga_model TRUE. (Also listed under STAD.)
- PMC8670907: fitting cohort not stated, GEO inferred -> not trained on TCGA.
- PMC8635899, PMC12278666, PMC12872855, PMC9440785: validation cohorts differ between sections; what the text ties to validation recorded.
- IMvigor210 in validation_used only where survival assessed: PMC9271306, PMC8907673, PMC10829884.

## c05_BRCA (214 screened, 84 eligible, trains_tcga_model TRUE 200)
- PMC10097725 corrected by the reader: "no prognostic model" -> "lncRNA model" (ALMS1-IT1), trains_tcga_model TRUE.
- PMC12480286: TMEM65 6.355 × 10^-5 recorded 6.355e-05 (script flag explained).
- Eligible, to look at: PMC9035888 (Table 2 coefficients look unsigned; HRs imply −0.565 IFNG, −0.197 TP63; as printed) TO RESOLVE; PMC8888901 ("status of" read as expression per Methods); PMC6997532 (text vs cited Table 2 for RAC2 −0.180 vs −0.081 -> Table 2; "LCAM" -> L1CAM); PMC9019977 (IGLL5 treated as protein-coding); PMC11953232 (101-combination ML search, chosen StepCox+Enet printed linear -> eligible; relevant to round 16); PMC9437322, PMC9479131 (TCGA merged with GEO/METABRIC); PMC8397558 ("IFNA" -> IFNA13, "NKX2 − 3" -> NKX2-3 from the same paragraph); PMC12640789 (all positive though some genes described protective).
- Ineligible, to look at: PMC9357945 (6-gene model, coefficients for 4), PMC10718145 (scrambled formula), PMC9352481 (age term in the model), PMC10905160 (GBP1P1 pseudogene), PMC9243266 (0/1 status).
- Training unclear: PMC8793540 (LASSO in GSE86374 -> trains_tcga_model FALSE), PMC10374203 (GSE42568 training per Methods vs abstract -> not trained on TCGA), PMC9465332 (not stated, TCGA-TNBC assumed).

## c06_STAD (183 screened, 77 eligible, trains_tcga_model TRUE 169)
- PMC7982000 (eligible): printed weights are exactly Table 3's multivariate HRs (0.49, 1.75, 0.55, 1.7) -> as printed. TO RESOLVE: HR-as-coefficient rule.
- PMC7593539: text CTLA4 −0.1370 vs cited Table 1 −0.370 (ln HR 0.691) -> Table 1 (OV PMC9050317 precedent).
- PMC10148812 (eligible): 24-gene StepCox + Enet from 98 ML combinations, printed linear -> eligible (relevant to round 16).
- Symbol corrections from the same paper: PMC9465170 ANXA5 printed twice, sixth term AXL; PMC10008502 "INGBE" -> INHBE; PMC9452725 "NFC1" -> NCF1; PMC12280864 "(HSP)A9" -> HSPA9.
- Main model: PMC8642032 13-gene RiskScore (not the PCA FRGScore printed first); PMC7880323 first model (MSI-L/MSS); PMC10990040 "multiplied by the exponential of IL6" read as expression.
- PMC9008303 ineligible (age and gender terms). PMC12515644 ineligible (SEER clinical nomogram), trains_tcga_model TRUE (TCGA-derived ARSB/PDCD1 ratio).
- Training merged with GSE84433/GSE84437: PMC9719959, PMC11655940, PMC12835628, PMC8728957, PMC9452725, PMC10669360, PMC10283039.
- Other-cancer TCGA validation listed: TCGA-PAAD (PMC9485603), TCGA-COAD+READ (PMC11760997).
- Accession typos as named with notes: PMC10522364 "GSE88437" (probably GSE84437); PMC9465170 "GSE15495" -> GSE15459 (Methods); PMC12776241 GSE87211; PMC9795222 ACRG labelled GSE64951 -> "ACRG". PMC11103507: third database not named.

## c02_LUAD_b (161 screened, 69 eligible, trains_tcga_model TRUE 155; validation_used filled 103)
- PMC10292897: IGHD (table names it "Idiopathic growth hormone deficiency") -> ineligible as Ig segment.
- PMC9302251: TCGA 3-gene model taken (presented first); title/conclusion favor a GEO-trained 7-gene model. TO RESOLVE with the "main model" rule.
- PMC12043920: first "prognostic model" is diagnostic -> the following Cox survival model taken.
- Symbol corrections from the paper's own text: PMC10457049 "SLC3A11" -> SLC3A2; PMC10292873 "SERINE2" -> SERPINE2; PMC8753132 ID2 split around a figure caption.
- Linear formula with RF-chosen genes or score: PMC10363522, PMC11735815 -> eligible.
- Odd but as printed: PMC11095029 (all ±0.001), PMC10051631 (two HR directions disagree with signs), PMC10085621 (all 12 positive).
- TCGA merged with GEO: PMC10051631 (GSE72094), PMC12170285 (GSE68571).
- PMC8059184: training "the LUAD dataset", TCGA implied.
- Pseudogenes ruled ineligible (PMC11496652 CYP4Z2P; PMC8643173 PSMD10P1); snoRNA not linear (PMC8649667). TO RESOLVE: pseudogene rule (BRCA reader did the same for GBP1P1).
- Checkpoint cohorts in validation_used only as OS validation: PMC8017122 (IMvigor210, GSE78220, GSE91061).
- Accessions as printed though suspect: PMC8497135 GSE13507 (bladder), PMC12133857 GSE5081 (probably GSE50081), PMC9367639 "TCGA-NSCLC".

## c08_HNSC_LUSC (162 screened: HNSC 122, LUSC 40; 64 eligible = HNSC 46 + LUSC 18; trains_tcga_model TRUE 157)
- Reader's own first pass missed formulas printed as TeX or with spaced letters; re-scan of non-lncRNA ineligible rows made PMC8851317, PMC8969218, PMC8249947 eligible and PMC12746546 "lncRNA model" (MIR9-3HG). -> SYSTEMIC RECALL CHECK needed in all chunks.
- PMC12575726: abstract 15 genes (incl. PLAU), formula 14 -> 14; "=" before CEACAM1 read "+"; PIK3R2 5.2505 as printed.
- PMC9523032: XML drops operators before KIF13B, FCGBP -> minus inferred (-0.446, -0.238).
- PMC7493716: text NRTN +0.1896 vs cited Table 6 -0.1896 -> table. PMC11220558: text multivariate signs opposite to univariate Table 3 -> text.
- Main model: PMC9948032 IRGPI (not composite); PMC10324233 cell death index (not the autophagy signature listed first); PMC8576700 Model-2 (not Model-1 printed first). TO RESOLVE with the "main model" rule.
- Symbols: PMC7485220 "GCCR" -> GCGR (GCCR is an NR3C1 alias); PMC8817543 "PLG1" -> PLCG1 (GSE41643 probably GSE41613); PMC10344514 protein names -> CLDN18, ADAMTS8; PMC9715950 old symbols kept (C9orf144B -> SPATA31F1, coding status unverified).
- PMC8109104: 0/1 median coding -> ineligible. PMC6969412: no explicit formula, LASSO coefficients from Table 1; GSE65858 in training unclear. PMC11382407: nasopharyngeal paper with a TCGA-HNSC model -> counted.
- Pooled training: PMC12321723 TCGA LUAD+LUSC counted for LUSC (ALSO counted for LUAD by c01) -> TO RESOLVE: one model under two types. PMC12168404 TCGA-LUSC + GSE50081.
- validation_used leaves out survival checks in checkpoint cohorts of other cancers.
