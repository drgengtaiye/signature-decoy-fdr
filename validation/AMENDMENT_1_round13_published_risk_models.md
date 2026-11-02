# Amendment 1 to PLAN_round13_published_risk_models.md: the search

Written 2026-11-02, during the first screening pass of part A and before any published model, theme pool or decoy was
scored in any survival cohort. Its SHA-256 is recorded in AMENDMENT_1_round13_published_risk_models.sha256 and it is
deposited before the amended search is run.

## What was seen

The search of A1 was run in PubMed Central (db = pmc), where the type and theme terms match anywhere in the full text.
It returned 2,697 (SARC) to 15,889 (BRCA) records per type, most of them about other cancers that mention the type in
passing: of the first 30 COAD records screened, 19 were ineligible for type mismatch and 1 was eligible; of the first
16 STAD records, 1 was eligible; about 13 of the first 80 COAD records were about colorectal cancer. At that rate the
cap of 80 screened records would yield one to three models per type. PMC also returns at most 9,999 identifiers per
query. Screening read only titles, abstracts and the papers' own text; no cohort was opened and nothing was scored.

## Change

A1 is replaced by: PubMed (db = pubmed) through E-utilities, run once, per type:
`(TCGA OR "The Cancer Genome Atlas") AND ("risk score" OR "risk model" OR "prognostic signature" OR "prognostic model")
AND (<type terms, each tagged [tiab]>) AND (<theme terms, each tagged [tiab]>) AND 2019:2025[dp] AND
"pubmed pmc open access"[filter]`, with the type and theme terms of the plan unchanged except for the [tiab] tag on
each. The TCGA and risk-score terms still match all fields. PMIDs are converted to PMCIDs through the PubMed summary
(article ids); the full text is fetched from PMC as before.

A2 then orders PMIDs (not PMCIDs) sorted ascending and shuffled with numpy default_rng(20261113 + type index), screened
until 10 eligible models or 80 records. A3 to A10 and part B are unchanged.

The superseded search outputs and the 46 records screened under it (STAD 16, COAD 30) are kept in
validation/round13_search/superseded_query1/ and are not used; a paper found by both searches is screened afresh under
the amended order.
