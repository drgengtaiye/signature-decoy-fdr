# Amendment 3 to PLAN_round15_literature_cohorts_controls.md: the validation cohorts (B1-B2)

Written 2026-11-03, after the candidates of B1 were listed and their GEO sample annotations read, and before any
expression value of a new cohort was downloaded or any null model, published model or assay was scored in one. The null
models of B3 had been fitted in TCGA alone (applications/v20_null_models.R; 1,060 null models from 1,320 decoy pools).
Its SHA-256 is recorded in AMENDMENT_3_round15_literature_cohorts_controls.sha256 and it is deposited before the first
expression download.

## Candidates (B1)

validation/round15_cohorts/candidates_r15.py lists, per type, the cohorts named in `validation_used` by the records
that train a TCGA model (1,598 records of round 15 and the 294 of round 13): 92 candidates, 79 (type, GEO series) pairs
over 76 series and 13 others (candidates_r15.csv).

## Eligibility (B2)

The rule was applied to the GEO sample annotations (round15_cohorts/eligibility_r15.csv, made by
round15_cohorts/tmp_b2/screen_b2.py): 41 of 79 pairs are eligible, from 40 series, and include the 11 round-11 cohorts.
The others fail with fewer than 80 patients or 25 events (21), no survival time or event in GEO (10), an Affymetrix
U133A platform below 15,000 genes (6), or the wrong tumor or data (4); a pair can fail for more than one reason. Of the 13 other candidates, METABRIC (BRCA),
ICGC OV-AU (OV; also the "ICGC" ovarian cohort of 11 records) and CGCI-HTMCP-CC (CESC, GDC open data) are retrievable
and are judged when built; IMvigor210 is checkpoint-treated; 8 are in-house or not identifiable
(eligibility_nongeo_r15.csv).

Decisions fixed in doing so:

1. A time of zero or less is not a survival time (the round-11 loader); a time derived as the difference of two
   annotated dates is (GSE41271: last follow-up minus surgery).
2. An event coding the series does not state is read as the one under which the events have the shorter median time
   (GSE14333 DFS_Cens 0 = recurrence; GSE72970 "os censored" 1 = death).
3. Tumors of the type follow the TCGA type's histology and site: serous carcinoma for OV (GSE26193 then has 79 patients
   and fails), no Ewing sarcoma for SARC (GSE17674, GSE63157), colon sites for COAD (rectal-only GSE87211 fails;
   rectosigmoid tumors of GSE72970 are left out).
4. The authors' own exclusion flags are not applied (GSE31210 keeps 226 patients): the literature uses series whole.
5. Platforms are checked for 15,000 genes after building; a cohort that falls below is dropped and reported.
6. Cohorts that share patients. Every eligible named cohort gets its own false-validation rate (B4), flagged with its
   overlap, because that is what a reader of a paper meets; two series with the same patients and the same data are
   built once (GSE40967's expression samples are GSE39582's; GSE103091's patients with survival are GSE58812's). In the
   tests that need independent cohorts (H4 Cochran's Q, the law of part C, the replication pairs of H5), each set of
   overlapping cohorts contributes only its cohort with most patients: GSE84437 (not GSE84433), GSE17538 (not GSE17536
   or GSE14333), GSE39582, GSE58812, GSE71118 (not GSE21050), GSE41271 (not GSE42127), GSE13213 (not GSE11969). GSE38832
   may share patients with GSE17537 but this cannot be checked and it is kept.
7. Non-GEO cohorts are built from their public files: METABRIC from cBioPortal (overall survival); OV-AU from the ICGC
   25K open release (exp_seq normalized read counts of primary solid tumour specimens, one per donor; overall survival
   from donor_survival_time and donor_vital_status); CGCI-HTMCP-CC from GDC open STAR counts (overall survival from
   vital status, days to death or to last follow-up). Each is held to the B2 rule.
## Frozen files (SHA-256)

    770c6d2b396fc2837a66fcae80440e28431f8a793986849bf4772659d43cc148  round15_cohorts/candidates_r15.csv
    d8ec79d7b56212612318821eb6511c57895bb37ab161002aaaf1912437865aac  round15_cohorts/eligibility_r15.csv
    f75a357063be86a42194ca8dece5d0ddcc9f8519dffb6b74283d93c7d227e65a  round15_cohorts/eligibility_nongeo_r15.csv
    c4479132b1d8646139e29b37c4ba43350a3db1c57b2f23c6235ce1016b4fabf8  round15_cohorts/cohorts_r15.json
