# Amendment 4 to PLAN_round15_literature_cohorts_controls.md: part D, before any assay is scored (2026-11-06)

Written after the assays' genes and algorithms were extracted and before any assay or decoy was scored. Only gene identifiers of the cohort files were read, to check coverage.

## The table
round15_assays/assays_r15.csv holds every printed row with its source, page or table and a verbatim quote; round15_assays/assays_notes_r15.md holds each algorithm verbatim. All six assays are kept. Frozen for scoring:
round15_assays/assays_r15_frozen.csv. SHA-256 prefixes: assays_r15.csv ade0c3919746309e; assays_r15_frozen.csv ea05d456cc91c01e; assays_alias_r15.csv aa440e7defc842dc; er_status_r15.csv f3f0022d3937d496; reference_genes_r15.csv c25b899f95415e4c; assays_notes_r15.md e2c323d4f21c8713.

The plan asks for the original publications. Three could not be retrieved or do not print everything, and the substitute is stated rather than hidden:
- 21-gene Recurrence Score: Paik 2004 was not retrievable; the algorithm is taken from the reprint in Habel 2006 (Fig. 1), a study with the developer's co-authors.
- PAM50 ROR-S: coefficients from Parker 2009 as printed; centroids from the developers' file for that paper (Internet Archive capture), because the journal's supplement could not be retrieved.
- Breast Cancer Index: the polynomial from Habel 2013, a reprint of Jerevall 2011; the current linear coefficients are unprinted, so the control is the earlier version of the index. MGI definition from Ma 2008.
Original: 70-gene (van 't Veer 2002 supplement), EndoPredict (Filipits 2011), 12-gene colon score (O'Connell 2010).

## Approximations, as plan D1 requires
Each score is the sum of weight x within-cohort z-score, higher meaning worse, as rule A6 of round 13. Reference genes, intercepts, rescaling, clipping and risk categories are dropped.
- 21-gene: group coefficients expanded to genes; the floors on the HER2 and proliferation groups are dropped.
- 12-gene colon: the same; the floors on FAP, MKI67, MYBL2 and GADD45B are dropped.
- EndoPredict: equation (B) as printed; EPclin (tumor size, nodes) is dropped.
- 70-gene: weight = −r, the published correlation with outcome, in place of the unprinted template and cosine correlation. 58 of 70 reporters have a gene symbol, giving 54 genes (repeats averaged); 12 reporters without one are left out. r's sign is anchored on the 18 and 52 of Marchionni 2013.
- PAM50: ROR-S by replacing each Spearman correlation with a dot product: weight = 0.05 C_basal + 0.12 C_her2 − 0.34 C_luma + 0.23 C_lumb, the printed coefficients (the developers' code has 0.142 for LumB). Tumor size, proliferation and rank correlation are dropped. Basal markers get negative weights: that is what the linearization gives.
- Breast Cancer Index: first-order terms of the cubic, f ≈ 0.4431 mgi + 0.4972 hi, hi = z(HOXB13) − z(IL17RB); the cubic term is dropped. MGI is the first principal component (Ma 2008) with equal loadings assumed, unit norm, so each of the five genes has 0.4431/√5 = 0.198. This scale was chosen before scoring, on that definition; loadings are not printed.

These are linear stand-ins for commercial assays. A result here is about these sums, not about the assays as sold.

## Decisions on coverage and cohorts
1. A gene that a cohort file carries under an older symbol counts as measured (org.Hs.eg.db aliases mapping to that gene alone; assays_alias_r15.csv); the current symbol wins when both exist.
2. Primary for H7: TCGA and the cohorts of the independent set of amendment 3 (BRCA: GSE20685, GSE42568, GSE58812, GSE21653, GSE20711, GSE96058, METABRIC; COAD: GSE39582, GSE17538, GSE38832, GSE103479, GSE72970). Plan D2 says "the eligible cohorts", and these drop overlapping patients, which pooling would count twice. The overlapping built cohorts are scored and reported as a sensitivity. GSE58812 is triple-negative, and it stays, as registered.
3. ER-positive patients, secondary: annotated cohorts with at least 25 events among ER-positive patients. These are GSE21653, GSE96058 and METABRIC. TCGA, GSE20685 and GSE58812 carry no usable annotation or no ER-positive patients in the files on disk; GSE20711 and GSE42568 have fewer than 25 events.
4. Assay genes, and decoys, are read per cohort at the 80% rule of round 13 A5.

H7, its bootstrap and its seed are as registered. With at most six assays the bootstrap over assays is coarse; this is a stated limitation of H7.
