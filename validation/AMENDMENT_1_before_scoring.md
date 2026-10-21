# Amendment 1 to the external-validation registration — written before any scoring

Registration: `PREREGISTRATION_external_validation.md`, SHA-256
b9918beadbad3c99aca2f5d4341a74ca58e439b0cba8df3493ea7ff011c19059, deposited 2026-09-18T23:26:04Z.
This amendment changes no endpoint, no analysis and no threshold. It records what was found on
downloading the expression files, before any signature was scored and before any expression
value was related to any outcome.

1. **Section 3's whole-transcriptome rule applied.** The GEO files of GSE161537 (2,559 rows),
   GSE159067 (2,567 rows) and GSE303300 (`htg_log2cpm`) are HTG EdgeSeq targeted panels of about
   2,550 genes, although GEO lists them as RNA-seq (the panel is read out on a sequencer). Under
   Section 3 — fewer than 3,000 genes is not a whole-transcriptome cohort and cannot enter the
   primary analysis — all three move to the panel set, which is analysed as registered (null
   drawn from the panel's own genes; reported separately, not combined with the primary set).
2. **The primary set is therefore GSE176307 alone** (metastatic urothelial carcinoma; Salmon
   gene-level TPM over 58,387 genes; 89 evaluable, CR/PR 16 : SD/PD 73), provided it passes the
   3,000-gene rule after filtering. The overall primary test (Stouffer over primary cohorts)
   reduces to that cohort's one-sided test.
3. **No other cohort qualifies.** The GEO supplementary files of the remaining
   whole-transcriptome candidates without response in their sample annotations were checked:
   GSE110390 (21 genes only), GSE212549, GSE282471, GSE285963, GSE302495, GSE202069, GSE245372,
   GSE276284, GSE313101, GSE281729, GSE206422, GSE235092, GSE213065, GSE318308 and GSE298296
   carry no RECIST best response in GEO; GSE183924 is adjuvant (no measurable disease);
   GSE202361 has 2 responders (fewer than 5). None is added.
4. **Panel set, as now constituted** (secondary): GSE161537, GSE159067, GSE303300 (HTG EdgeSeq
   precision immuno-oncology panel), GSE215868 (NanoString IO 360) and GSE93157 (NanoString
   PanCancer Immune Profiling). The panel-set summary combines its cohorts by Stouffer's method
   exactly as the primary set would have been, and is labelled secondary wherever it appears.
