# Result of PLAN_round14_new_checkpoint_cohorts.md: no eligible cohort

Recorded 2026-11-02. The screen of section A ran on the 303 shortlisted series and the three nominated by the authors
(GSE177043 was already in the shortlist; 305 candidates). Its table, with the criterion reached and the reason for every
candidate, is round14_search/screen.csv. No expression value of any candidate was downloaded or read.

No candidate met all six criteria of A2. Under section C ("if fewer than two cohorts are eligible the tests are not
run, and the screen is reported with its reasons") the predictions were not tested.

Stage reached: 254 candidates were excluded on title and summary (cell lines, mice, blood, single-cell series,
neoadjuvant settings, checkpoint blockade combined with other agents), 27 on the series record, 11 on the sample
annotations and 13 on the publication. Twelve passed criteria (1) and (2):

| Series | Cohort | Fails | Reason |
|---|---|---|---|
| GSE212549 | NIVACTOR, recurrent/metastatic HNSCC, nivolumab; 80 tumors, Clariom D (GSM tables) | (3) | 12 PR against 67 SD/PD in the paper, but no per-patient label in GEO, the PMC text or the supplement (eTable 4 is aggregate; read through Europe PMC's supplementary-file service, because PMC's download page asks for a browser check) |
| GSE160638 | melanoma, PD-1 blockade, 36 biopsies | (3) | no response labels in GEO and no PMC full text |
| GSE110390 | durvalumab (Study 1108), NSCLC and urothelial | (3) | no labels; processed data cover 21 genes |
| GSE235092 | Merkel cell carcinoma, 63 patients (31 CR/PR in the paper) | (3), (4) | labels only in a publisher supplementary table; expression only in supplementary files on the GEO FTP server, which is unreachable from this cluster |
| GSE276284 | ovarian clear-cell, durvalumab + tremelimumab, 29 | (3), (4) | as GSE235092 |
| GSE215868 | melanoma, 104 (45 CR/PR) | (4) | NanoString IO360, 770 genes |
| GSE161537 | NSCLC, 76 (20 CR/PR) | (4) | HTG EdgeSeq panel, 2,549 genes |
| GSE224450 | nasopharyngeal, nivolumab + ipilimumab | (4), (5) | NanoString; 22 patients |
| GSE181815 | thymic carcinoma | (5) | 10 patients |
| GSE190265, GSE190266, GSE78220 | NSCLC, melanoma | (6) | already in the collection (Limagne, Hugo) |

Near misses with (1) or (2) unclear: GSE283829 (NSCLC; regimen not stated), GSE318164 (upper-tract urothelial; labels and
expression supplementary only), GSE274975 (NSCLC; mixed regimens, labels in the supplement), GSE282471 (melanoma;
"clinical benefit" counting stable disease and long survival), GSE273829 (renal; arms include VEGFR inhibitors).

Obtainability is a limit of this cluster as well as of the data: series whose expression exists only as supplementary
files (GSE235092, GSE276284, GSE318164, GSE283829, GSE274975) fail (4) here because the GEO FTP server cannot be
reached; GSE235092 would also need per-patient labels from its publisher's supplement. Controlled-access trial data
(JAVELIN Renal 101, IMmotion150/151, POPLAR/OAK, IMbrave150, SU2C-MARK, Hartwig) need data-access approval. Either
would enter by amendment under the plan's rules unchanged.

During screening the reader saw, incidentally, the GSE274975 paper's own biomarker table, the gene counts of a
differential-expression paragraph of the GSE235092 paper and a depositor-computed "HOT" score in the GSE161537 and
GSE159067 annotations; none concerns the 52 signatures and none was used.
