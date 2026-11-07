# Clinical positive controls, round 15, part D: assay extraction notes

Status: extraction only. Nothing was scored. All approximations below are proposals for acceptance or change.

All six assays are kept, none is dropped.

| assay (csv name) | source | kind | rows | weights |
|---|---|---|---|---|
| 21-gene RS (RS21_breast) | Habel 2006, PMID 16737553, Fig 1 | reprint of Paik 2004, PMID 15591335 | 16 | group coefficients printed, per-gene values derived by arithmetic |
| 70-gene (MammaPrint70) | van 't Veer 2002, PMID 11823860, Supp Table 2 and Methods | original | 70 reporters | approximated, -r from Table S2 |
| EndoPredict | Filipits 2011, PMID 21807638, Eq. B | original | 8 | printed per gene |
| PAM50 ROR (PAM50_ROR) | Parker 2009, PMID 19204204, Eq. 1; centroids from the developers' file | coefficients original, centroids reprint | 50 | approximated, derived from printed coefficients and centroids |
| Breast Cancer Index (BCI) | Habel 2013, PMID 23497539 (polynomial); Ma 2008, PMID 18451222 (MGI) | reprint of Jerevall 2011, PMID 21559019, for the polynomial | 7 | H/I approximated; 5 MGI weights left blank |
| 12-gene colon RS (RS12_colon) | O'Connell 2010, PMID 20679606, Fig 3 legend | original | 7 | group coefficients printed, per-gene values derived by arithmetic |

## Conventions

- `weight` is per unit of gene level in a sum where higher means worse. `direction` is the sign of that, or for blank weights the printed direction.
- Scales differ between assays and I did not harmonise them. RS21, RS12 and EndoPredict are per reference-normalised log2 unit. BCI H/I is per standardised unit (mean 0, SD 1). MammaPrint -r is unitless. PAM50 weights multiply gene-centred expression in the centroid scale. A sum across genes is fine within an assay; weights are not comparable between assays.
- `gene_symbol` is the current NCBI Gene symbol as retrieved on 2026-11-06 (raw results in `external/raw_r15/assays/NCBI_gene_alias/` and `NCBI_lookup_mp70/`). Cohort files may carry older names. Printed to current, outside MammaPrint: Cathepsin L2 / CTSL2 to CTSV; STK15 to AURKA; IL17BR to IL17RB; CDCA1 to NUF2; KNTC2 to NDC80; ORC6L to ORC6; reference ATP5E to ATP5F1E. Table of the MammaPrint names is at the end of its section.
- `quote`: every text quote in the two csv files was machine-checked as a substring of the retrieved source after whitespace normalisation, and every table row (MammaPrint, PAM50) was copied by script from the retrieved file with cells joined by one space. The RS21 quotes come from Figure 1 of Habel 2006, which is an image: I read it, and its coefficients agree with the text of PMC7437067 and PMC3978487. Equations of Filipits 2011 and Parker 2009 are also images; I read them and they agree with the text of the Filipits Supplementary Appendix and of Cheang 2012 (PMC3743660).
- Reference genes are in `reference_genes_r15.csv`, not in `assays_r15.csv`. MammaPrint has none: it uses log ratios to a reference pool on a two-colour array.

## 21-gene Recurrence Score (Oncotype DX)

Source: original Paik 2004 NEJM (PMID 15591335) could not be retrieved: `nejm.org` returned 403 to curl and to `WebFetch`, and it is not in PMC. Used the reprint in Habel 2006, Breast Cancer Res (PMC1557737), a study with Genomic Health co-authors that cites Paik 2004. Algorithm, Figure 1, as printed in the image:

- HER2 Group Score: "(0.9 × GRB7) + (0.1 × HER2). If HER2 Group Score is less than 8 then the HER2 Group Score is considered equal to 8"
- ER Group Score: "([0.8 × ER] + [1.2 × PgR] + Bcl2 + SCUBE2])/4"
- Proliferation Group Score: "(SURV + KI-67 + MYBL2 + Cyclin B1 + STK15)/5. If the Proliferation Group Score is under 6.5 then the Proliferation Group Score is considered equal to 6.5"
- Invasion Group Score: "(Cathepsin L2 + Stromelysin 3)/2"
- "RSU = + 0.47 × HER2 Group Thresholded Score – 0.34 × ER Group Score + 1.04 × Proliferation Group Thresholded Score + 0.10 × Invasion Group Score + 0.05 × CD68 – 0.08 × GSTM1 – 0.07 × BAG1"

Transformations: "normalized relative to the expression of the five reference genes (b-actin, GAPDH, GUS, RPLPO, and TFRC)"; "A plus sign indicates increased expression is associated with increased recurrence risk.". Categories: "low risk (Recurrence Score <18), intermediate risk (Recurrence Score 18–30) and high risk (Recurrence Score ≥31)". The rescaling is not in Figure 1. Bernhardt 2020 (PMC7437067) prints "RS = 0 if RSu < 0 RS = 20 × RSu - 6.7 if 0 ≤ RSu ≤ 100 RS = 100 if RSu > 100" (its LaTeX form has parentheses, 20 × (RSu − 6.7)).

Proposed sum: expand the group scores, RS = sum of coefficient × within-group weight × gene. Example: GRB7 0.47 × 0.9 = 0.423; PGR -0.34 × 1.2 / 4 = -0.102; each proliferation gene 1.04 / 5 = 0.208.

Drops: the floor at 8 on the HER2 group and at 6.5 on the proliferation group; the rescale to 0-100 and its clipping; the risk categories.

Not settled: the rescale constant comes only from a later paper, and does not affect weights. Bernhardt 2020 (four proliferation genes, /4) and Balko 2012 (PMC3693569, three genes, states two were unavailable) print shortened groups; I used the five-gene form of Habel, which agrees with Madden 2013 Table 3 (PMC3978487).

## 70-gene signature (MammaPrint)

Source: van 't Veer 2002, Nature (PMID 11823860). The article body is behind a paywall. Its Supplementary Methods and Tables were retrieved from `static-content.springer.com`. The gene-defining statement is original: "The 70 optimal marker genes are the 70 genes with the highest absolute correlation coefficients." I took the top 70 by absolute r of Table S2 (`415530a-s9.xls`), the smallest being 0.356001 (NM_018354), the next 0.355105 (NM_000286).

Algorithm, Supplementary Methods: "We calculated the correlation between the prognostic category (metastasis vs. no-metastasis) and the logarithmic expression ratio across all 78 samples for each individual gene". Classification, as restated by the developers in Glas 2006 (PMC1636049): "the 'good prognosis template' was (re-)constructed using the average expression for each of the 70 genes in tumors from the 44 'good outcome' patients"; "A patient with a cosine correlation to the good prognosis template higher than 0.4 (the previously determined threshold [ 6 ]) was assigned to the good-profile group." For a van de Vijver subgroup the threshold used is "the previously determined threshold of 0.55".

Direction: the sign convention of r is not stated in the retrieved supplement. It is anchored by Marchionni 2013 (PMC3662649): "18 genes had expression values positively associated with good prognosis, while 52 were associated with metastatic recurrence." The top 70 of Table S2 have exactly 18 r > 0 and 52 r < 0, and the genes it names as good (GSTM3, MS4A7) have r > 0 and as bad (OXCT1, HRASLS, MELK, Contig40831_RC) r < 0. So direction = -sign(r).

Proposed sum: weight = -r from Table S2. The real algorithm correlates a tumour with a template of average log ratios in 44 good-outcome tumours; those template values are not printed anywhere I could retrieve, so the published outcome correlation stands in for them. Alternative for you to consider: unit weights +1/-1 from the same signs.

Drops: the template, the cosine correlation (normalisation by the tumour's own vector), the 0.4 threshold, and the 78-patient training data.

Not settled:
- 12 rows have a weight but a blank `gene_symbol`. Eleven are reporters I could not link to a gene: Contig63649_RC, Contig38288_RC, Contig32125_RC, Contig55377_RC, Contig35251_RC, Contig48328_RC, Contig40831_RC, Contig51464_RC, Contig24252_RC, Contig56457_RC, Contig20217_RC. The developers' gene-level table (Tian 2010, PMC2999994, Table 1) leaves exactly eleven names unassigned: LIN9, RASSF7, RECQL5, RUNDC1, JHDM1D (KDM7A), ZNF533 (ZNF385B), QSCN6L1 (QSOX2), C9orf30 (MSANTD3), LOC100288906, LOC730018, LOC100131053. So the set is known and the pairing is not. I did not pair them by guess. Table S2 prints Contig51464_RC as FLJ22477 and Contig56457_RC as TMEFF1, but neither name is in Tian's list as such. The twelfth, AA555029_RC, is listed by Tian as "No significant similarity found". Sequence mapping failed: NCBI BLAST reset the TLS connection from this cluster and the UCSC search API timed out. ESTs have no NCBI Gene link.
- Rows are reporters, not genes, because the assay is defined on 70 reporters. Three genes repeat: DIAPH3 (3 reporters), IGFBP5 (2), ECI2 (2). Proposed: average within symbol at scoring. 58 rows carry 54 distinct symbols, all in Tian's list; 54 + 12 = 66 = Tian's count of distinct names.
- Contig28552_RC to DIAPH3 rests on its Table S2 description (same clone as AL137718) plus "DIAPH3, which was unidentified in the original paper, appears three times in the 70 gene set."
- Hedjazi 2013 (PMC3728725) flags Contig11065_RC and Contig49388_RC as in the 70-gene list; they are not in the 231 of Table S2. I followed Table S2.
- The Agendia array design file (A-MEXP-318) was checked and not used: it gives identical UniGene and EMBL entries to unrelated reporters.

| rank by abs(r) | reporter | Table S2 name | r | symbol | evidence for symbol |
|---|---|---|---|---|---|
| 1 | AL080059 |  | -0.527145 | TSPYL5 | accession to NCBI Gene (nuccore elink) |
| 2 | Contig63649_RC |  | -0.468129 | (blank) | not settled |
| 3 | Contig46218_RC |  | -0.432539 | DIAPH3 | Carrivick 2006 Table 4 |
| 4 | NM_016359 | LOC51203 | -0.424927 | NUSAP1 | accession to NCBI Gene (nuccore elink) |
| 5 | AA555029_RC |  | -0.424122 | (blank) | not settled |
| 6 | NM_003748 | ALDH4 | 0.420671 | ALDH4A1 | accession to NCBI Gene (nuccore elink) |
| 7 | Contig38288_RC |  | -0.414971 | (blank) | not settled |
| 8 | NM_003862 | FGF18 | 0.410964 | FGF18 | accession to NCBI Gene (nuccore elink) |
| 9 | Contig28552_RC |  | -0.409259 | DIAPH3 | Table S2 description clone DKFZp434C0931, same as AL137718; Carrivick 2006: DIAPH3 appears three times in the 70 |
| 10 | Contig32125_RC |  | 0.409054 | (blank) | not settled |
| 11 | U82987 | BBC3 | 0.407002 | BBC3 | accession to NCBI Gene (nuccore elink) |
| 12 | AL137718 |  | -0.404979 | DIAPH3 | accession to NCBI Gene (nuccore elink) |
| 13 | AB037863 | KIAA1442 | 0.402335 | EBF4 | accession to NCBI Gene (nuccore elink) |
| 14 | NM_020188 | DC13 | -0.400072 | CMC2 | accession to NCBI Gene (nuccore elink) |
| 15 | NM_020974 | CEGP1 | 0.399987 | SCUBE2 | accession to NCBI Gene (nuccore elink) |
| 16 | NM_000127 | EXT1 | -0.399515 | EXT1 | accession to NCBI Gene (nuccore elink) |
| 17 | NM_002019 | FLT1 | -0.398066 | FLT1 | accession to NCBI Gene (nuccore elink) |
| 18 | NM_002073 | GNAZ | -0.395461 | GNAZ | accession to NCBI Gene (nuccore elink) |
| 19 | NM_000436 | OXCT | -0.392119 | OXCT1 | accession to NCBI Gene (nuccore elink) |
| 20 | NM_004994 | MMP9 | -0.391694 | MMP9 | accession to NCBI Gene (nuccore elink) |
| 21 | Contig55377_RC |  | 0.3906 | (blank) | not settled |
| 22 | Contig35251_RC |  | -0.390407 | (blank) | not settled |
| 23 | Contig25991 | ECT2 | -0.390367 | ECT2 | Table S2 gene name ECT2 |
| 24 | NM_003875 | GMPS | -0.386515 | GMPS | accession to NCBI Gene (nuccore elink) |
| 25 | NM_006101 | HEC | -0.385893 | NDC80 | accession to NCBI Gene (nuccore elink) |
| 26 | NM_003882 | WISP1 | 0.384479 | CCN4 | accession to NCBI Gene (nuccore elink) |
| 27 | NM_003607 | PK428 | -0.384392 | CDC42BPA | accession to NCBI Gene (nuccore elink) |
| 28 | AF073519 | SERF1A | -0.38334 | SERF1A | accession to NCBI Gene (nuccore elink) |
| 29 | AF052162 | FLJ12443 | -0.380834 | LPCAT1 | accession to NCBI Gene (nuccore elink) |
| 30 | NM_000849 | GSTM3 | 0.380831 | GSTM3 | accession to NCBI Gene (nuccore elink) |
| 31 | Contig32185_RC |  | -0.379169 | GPR180 | Marchionni 2013 Additional file 1 Table 10 |
| 32 | NM_016577 | RAB6B | -0.376227 | RAB6B | accession to NCBI Gene (nuccore elink) |
| 33 | Contig48328_RC |  | 0.375252 | (blank) | not settled |
| 34 | Contig46223_RC |  | 0.374289 | RTN4RL1 | Marchionni 2013 Additional file 1 Table 10 |
| 35 | NM_015984 | UCH37 | -0.373876 | UCHL5 | accession to NCBI Gene (nuccore elink) |
| 36 | NM_006117 | PECI | 0.37329 | ECI2 | accession to NCBI Gene (nuccore elink) |
| 37 | AK000745 |  | -0.373061 | MTDH | accession to NCBI Gene (nuccore elink) |
| 38 | Contig40831_RC |  | -0.372931 | (blank) | not settled |
| 39 | NM_003239 | TGFB3 | 0.371524 | TGFB3 | accession to NCBI Gene (nuccore elink) |
| 40 | NM_014791 | KIAA0175 | -0.370858 | MELK | accession to NCBI Gene (nuccore elink) |
| 41 | X05610 | COL4A2 | -0.370857 | COL4A2 | accession to NCBI Gene (nuccore elink) |
| 42 | NM_016448 | L2DTL | -0.369416 | DTL | accession to NCBI Gene (nuccore elink) |
| 43 | NM_018401 | HSA250839 | 0.368349 | STK32B | accession to NCBI Gene (nuccore elink) |
| 44 | NM_000788 | DCK | -0.367698 | DCK | accession to NCBI Gene (nuccore elink) |
| 45 | Contig51464_RC | FLJ22477 | -0.367446 | (blank) | not settled |
| 46 | AL080079 | DKFZP564D0462 | -0.367388 | ADGRG6 | accession to NCBI Gene (nuccore elink) |
| 47 | NM_006931 | SLC2A3 | -0.366486 | SLC2A3 | accession to NCBI Gene (nuccore elink) |
| 48 | AF257175 |  | 0.3659 | ECI2 | accession to NCBI Gene (nuccore elink) |
| 49 | NM_014321 | ORC6L | -0.365812 | ORC6 | accession to NCBI Gene (nuccore elink) |
| 50 | NM_002916 | RFC4 | -0.365585 | RFC4 | accession to NCBI Gene (nuccore elink) |
| 51 | Contig55725_RC |  | -0.365347 | CDCA7 | Marchionni 2013 Additional file 1 Table 10 |
| 52 | Contig24252_RC |  | -0.364986 | (blank) | not settled |
| 53 | AF201951 | CFFM4 | 0.363953 | MS4A7 | accession to NCBI Gene (nuccore elink) |
| 54 | NM_005915 | MCM6 | -0.363852 | MCM6 | accession to NCBI Gene (nuccore elink) |
| 55 | NM_001282 | AP2B1 | 0.363326 | AP2B1 | accession to NCBI Gene (nuccore elink) |
| 56 | Contig56457_RC | TMEFF1 | -0.361645 | (blank) | not settled |
| 57 | NM_000599 | IGFBP5 | -0.361285 | IGFBP5 | accession to NCBI Gene (nuccore elink) |
| 58 | NM_020386 | LOC57110 | -0.360778 | PLAAT1 | accession to NCBI Gene (nuccore elink) |
| 59 | NM_014889 | MP1 | -0.360044 | PITRM1 | accession to NCBI Gene (nuccore elink) |
| 60 | AF055033 | IGFBP5 | -0.359944 | IGFBP5 | accession to NCBI Gene (nuccore elink) |
| 61 | NM_006681 | NMU | -0.359698 | NMU | accession to NCBI Gene (nuccore elink) |
| 62 | NM_007203 | AKAP2 | -0.359573 | PALM2AKAP2 | accession to NCBI Gene (nuccore elink) |
| 63 | Contig63102_RC | FLJ11354 | 0.359255 | DHX58 | Marchionni 2013 Additional file 1 Table 10 (LGP2); NCBI alias of DHX58 |
| 64 | NM_003981 | PRC1 | -0.358259 | PRC1 | accession to NCBI Gene (nuccore elink) |
| 65 | Contig20217_RC |  | -0.357877 | (blank) | not settled |
| 66 | NM_001809 | CENPA | -0.357716 | CENPA | accession to NCBI Gene (nuccore elink) |
| 67 | Contig2399_RC | SM-20 | -0.356602 | EGLN1 | Table S2 name SM-20; NCBI alias of EGLN1; Tian 2010 Table 1 EGLN1 |
| 68 | NM_004702 | CCNE2 | -0.356595 | CCNE2 | Table S2 gene name CCNE2; NCBI record for the 2002 accession returned no title |
| 69 | NM_007036 | ESM1 | -0.356542 | ESM1 | accession to NCBI Gene (nuccore elink) |
| 70 | NM_018354 | FLJ11190 | -0.356001 | TMEM74B | accession to NCBI Gene (nuccore elink) |

## EndoPredict

Source: Filipits 2011, Clin Cancer Res 17:6012 (PMID 21807638), original, full text and Supplementary Appendix retrieved from `aacrjournals.org`. p.6013, equations read from the images:

- (A) ΔCt(GOI) = 20 − Ct(GOI) + [Ct(CALM2) + Ct(OAZ1) + Ct(RPL37A)]/3
- (B) su = 0.41 · ΔCt(BIRC5) − 0.35 · ΔCt(RBBP8) + 0.39 · ΔCt(UBE2C) − 0.31 · ΔCt(IL6ST) − 0.26 · ΔCt(AZGP1) + 0.39 · ΔCt(DHCR7) − 0.18 · ΔCt(MGP) − 0.15 · ΔCt(STC2) − 2.63
- (C) s = 0, if 1.5 · su + 18.95 < 0; s = 15, if 1.5 · su + 18.95 > 15; s = 1.5 · su + 18.95, otherwise
- (D) sclin = 0.35 · t + 0.64 · n + 0.28 · s

"The EP risk score ranges from 0 to 15; higher values indicate a higher risk of recurrence." "Thresholds for EP and EPclin to discriminate patients into low and high risk of distant recurrence were developed in the training set at 5 and 3.3, respectively". ΔCt rises with expression because Ct(GOI) is subtracted, so the sign in (B) is direction directly. Check: the four intercepts of Supplementary Appendix p.8 (-1.20 + 1.19 - 1.63 - 0.99) sum to -2.63.

Proposed sum: equation (B) as printed, no approximation of gene weights.

Drops: the intercept -2.63, the scaling and clipping (C), the cut-offs, and EPclin (tumour size and nodes).

Not settled: nothing on weights.

## PAM50 risk of recurrence

Source: Parker 2009, J Clin Oncol 27:1160 (PMID 19204204, PMC2667820). The full text is read from the PMC page because the publisher does not allow XML download. Equations, read from the images:

- (1) ROR-S = 0.05 · basal + 0.12 · HER2 + −0.34 · LumA + 0.23 · LumB
- (2) ROR-C = 0.05 · basal + 0.11 · HER2 + −0.23 · LumA + 0.09 · LumB + 0.17 · T

Text: "A ROR score was assigned to each test case using correlation to the subtype alone (1) (ROR-S) or using subtype correlation along with tumor size (2) (ROR-C):" "the final algorithm consisted of centroids constructed as described for the PAM algorithm24 and distances calculated using Spearman's rank correlation." Equation (1) in text in Cheang 2012: "ROR-S = (0.05) × basal like + (0.12) ×HER2-E + (−0.34) × LumA + (0.23) × LumB", with "low risk if the ROR-S score was less than 23, moderate risk if ROR-S score was between 23 and 53, and high risk if the ROR-S score was 53 or more". ROR-C training cut-offs in Parker Appendix Fig A5: "low-risk scores less than 29 (black), moderate-risk scores between 29 and 53 (green), and high-risk scores ≥ 53 (red)".

Genes and centroids: gene names are legible in Parker Appendix Fig A3 (image) and match Okimoto 2024 Table 1 (PMC10908178) as a set. Centroid values are in the Data Supplement `JCO.2008.18.1370_supplemental_tables1-6.xls`, which I could not get: PMC answers curl and `WebFetch` with a browser proof-of-work page, which I did not try to defeat, and `ascopubs.org` returned 403. Used instead the developers' file `pam50_centroids.txt` from their site for this paper (`genome.unc.edu/pubsup/breastGEO`, host no longer resolves), by the Internet Archive capture of 2010-06-13. The copy inside `PAM50.zip`, captured 2013-02-24, is token-identical. Hence `source_kind` = reprint for these rows. By eye, the heat-map colours in Fig A3 agree in pattern (for example LumA high for ESR1, NAT1, SLC39A6). The values were not compared with the supplement, which I do not have.

Proposed sum (ROR-S chosen because it is gene-only): replace each Spearman correlation with a dot product, so weight(g) = 0.05 C_basal,g + 0.12 C_her2,g - 0.34 C_luma,g + 0.23 C_lumb,g. Applies to gene-median-centred log expression.

Drops: the rank correlation (and so the outlier robustness), the Normal centroid, the tumour size term of ROR-C and the proliferation term of ROR-P, the rescale to 0-100, the cut-offs, and subtype assignment itself.

Not settled:
- **The coefficients differ between the paper and the developers' own code.** `subtypePrediction_distributed.R` in the same archive has `genomic <- 0.04210193*basal + 0.12466938*her2 + -0.35235561*luma + 0.14213283*lumb`, with cut-offs -0.15 and 0.1 and rescale 100*(x + 0.35)/0.85. The LumB term is 0.142 there against 0.23 printed. I used the printed values, per the rule. Using the code values flips the sign of CDH3 (-0.101 to +0.027) and ERBB2 (+0.0007 to -0.0062) and changes magnitudes moderately. Your call.
- ERBB2 gets a near-zero weight (+0.0007), and basal markers get negative weights (KRT5 -0.309, KRT14 -0.321, KRT17 -0.358, FOXC1 -0.340), because their LumA, LumB and Normal contrasts dominate under this formula. That is what the linearisation gives, not a claim about those genes in the paper. FGFR4 (0.030) and TMEM45B (0.039) are also small. If you prefer not to rely on these signs, a sign-free alternative is to score ROR by actual centroid correlation, which the plan excludes unless you amend it.
- Reference genes (ACTB, PSMC4, RPLP0, MRPL19, SF3A1) come from US 9,631,239 B2 by the developers, because the retrievable part of Parker 2009 does not name them.
- Rejected source: Hamaneh 2023 (PMC10052604) prints 50 x 5 centroid values in text, but they are retrained, "the centroids are computed differently", so they are not this algorithm.
- A search result describes a centroid object in the Bioconductor package `genefu`; it was not downloaded because it is stored as `.rda`.

## Breast Cancer Index

Source. Genes: Zhang 2013, Clin Cancer Res 19:4196 (PMID 23757354), p.4197: "The genes analyzed were HOXB13, IL17BR (H/I), BUB1B, CENPA, NEK2, RACGAP1, RRM2 (MGI), ACTB, HMBS, SDHA, and UBC (reference genes)." Original cubic model: Jerevall 2011 (PMID 21559019, PMC3111159): "We used a polynomial function of H : I to approximate the restricted cubic spline", without printing its coefficients. Those are printed in Habel 2013 (PMC3672697), a Kaiser Permanente study funded by bioTheranostics that cites Jerevall 2011: "f h i , m g i = 0 . 4431 * m g i + 0 . 4972 * h i - 0 . 09 * h i 3" (that is, f = 0.4431 mgi + 0.4972 hi - 0.09 hi^3), with BCI = 0 if f < -2.5, 2f + 5 if -2.5 <= f <= 2.5, 10 if f > 2.5, and groups "low (BCI <5), intermediate (5 ≤ BCI <6.4), and high (BCI ≥ 6.4)".

Current linear version, Zhang 2013 p.4198: "The raw BCI risk index was produced by a linear combination of H/I and MGI weighed by the corresponding coefficients derived from a Cox regression model including both H/I and MGI as independent variables." with "low risk if BCI < 5.0825, intermediate risk if 5.0825 ≤ BCI < 6.5025, and high risk if BCI ≥ 6.5025". The coefficients are not printed there or in any paper I retrieved.

MGI, Ma 2008 (PMID 18451222) p.2602: "standardized to a mean of 0 and a SD of 1 across samples within each data set and then combined into a single index per sample as the first principal component". Cut-point "using the mean (0) as cutpoint", and "The cutpoint of 0.06 for HOXB13:IL17BR". Supplementary Tables 4 and 5 list the genes without loadings. Zhang 2013 p.4202 describes it as "MGI is the additive expression values of 5 cell-cycle genes".

H/I: original Ma 2006 (PMID 17008703) returned 403 from `ascopubs.org`, not in PMC; only its abstract was retrieved. The definition comes from the developers' patent US 11,021,754 B2 (assignee Biotheranostics; inventors include Erlander, Ma, Sgroi), as a patent and not peer reviewed: "The HOXB13:IL17BR ratio was calculated as the difference in standardized expression levels between HOXB13 and IL17BR as described previously". Goetz 2006 agrees on the first step: "which were z-transformed before taking the ratio".

Proposed sum: first-order terms of the cubic, f approx 0.4431 mgi + 0.4972 hi, with hi = z(HOXB13) - z(IL17RB). So HOXB13 +0.4972, IL17RB -0.4972, per standardised unit.

Drops: the cubic term -0.09 hi^3, whose gradient 0.4972 - 0.27 hi^2 changes sign at |hi| = 1.36, so the linear form misdescribes tumours with extreme H/I; the rescale and clipping; the cut-offs; the reference-cohort means and SDs for standardisation (not printed; each cohort's own would stand in).

Not settled:
- **The 5 MGI weights are blank.** Direction +1 is printed. The loadings are not printed anywhere retrieved, and the patent says "generation of the formula for the index is dependent upon the dataset, reference gene, and genes of the MGI". Assuming equal loadings, the scale relative to H/I still depends on what mgi is: mean gives 0.4431/5 = 0.0886 per gene, a unit-norm first component gives 0.4431/sqrt(5) = 0.198, a plain sum gives 0.4431. These differ five-fold against H/I, so I did not pick one. You must choose before freezing, or BCI scores on H/I alone.
- The coefficients are from the cubic model, not the current linear BCI, whose coefficients are unprinted. Accepting these means the control is the earlier version of the index.
- H/I weights are only meaningful on gene-wise z-scores, unlike the log2 scale of RS21, RS12 and EndoPredict.

## 12-gene colon Recurrence Score

Source: O'Connell 2010, J Clin Oncol 28:3937 (PMID 20679606, PMC2940392), original, Figure 3 legend on the PMC page, verbatim:

"Expression of individual genes is given a threshold as follows: if FAP, Ki-67 or MYBL2 measurement is less than 6.5 CT, it is considered to be 6.5; if GADD45B measurement is less than 5 CT, then it is considered to be 5."

"stromal group score = (BGN + FAP + INHBA) ÷ 3; cell cycle group score = (Ki-67 + C-MYC + MYBL2) ÷ 3."

"RSu = + 0.15 × stromal group score − 0.30 × cell cycle group score + 0.15 × GADD45B. A plus sign indicates that increased expression is associated with increased risk of recurrence, and a minus sign indicates that increased expression is associated with decreased risk of recurrence."

"RS = 44 × (RSu + 0.82); if RS is less than 0, then RS = 0, and if RS is greater than 100, then RS = 100."

Proposed sum: expand the groups. BGN, FAP, INHBA 0.15/3 = 0.05 each; MKI67, MYC, MYBL2 -0.30/3 = -0.10 each; GADD45B +0.15.

Drops: the floors on FAP, Ki-67, MYBL2 (6.5) and GADD45B (5); the rescale and clipping; risk groups; the separate treatment score.

Not settled: Yothers 2013 (PMC3871512) prints the rescale with the opposite sign, "RS = 44 × (RSu − 0.82)", against + 0.82 in O'Connell. It does not change weights, so I left it. Clark-Langone 2010 (PMC3016296) spells Ki-67 as MK167 in places and MKI67 in others.

## Retrieval record

Downloads are under `external/raw_r15/assays/<source>/`; scripts stayed outside them, in the session scratchpad. No e-mail or personal identifier was sent. What could not be retrieved, and what was done instead:

- PMC binary supplements (`/articles/instance/.../bin/...`): a browser proof-of-work page ("Preparing to download") to curl and `WebFetch`. Not bypassed. Affects Parker 2009 tables.
- `nejm.org`, `ascopubs.org`, `cell.com`: 403 challenge pages to curl; `nejm.org` also 403 to `WebFetch`. No browser impersonation was attempted. Affects Paik 2004, Ma 2006, Ma 2004. Nature's article body for van 't Veer 2002 is paywalled; only its open supplements were used.
- NCBI BLAST: TLS reset from this cluster. UCSC API search: timed out.
- `.rda` and `.rds` files: none loaded. The Bioconductor `mammaPrintData` and `seventyGeneData` source archives were downloaded and only text and `.xls` members extracted (`ArrayNomenclature_contig_accession.xls`, `ProbeSeq.xls`, vignettes). From the NKI 295 table I read only the substance and gene-name columns of the 70 reporters, not any patient value. Those archives were used for reporter annotation only, not as a source of weights.
- A guessed DOI and a guessed ScienceDirect PII returned wrong or missing pages on the way; those downloads were deleted and nothing from them was used.
- USPTO PDFs of the patents are scanned images with no text layer. The patent text quoted here is from Google Patents (`GooglePatents_*`).
