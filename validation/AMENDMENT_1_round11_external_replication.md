# Amendment 1 to PLAN_round11_external_replication.md — written after the downloads, before any scoring

Plan SHA-256 14f0567f6ffcbb5d1a9288e6e61a2cedbaa3f107a2de3a4908bc55eeac7e5acb, deposited 2026-11-01T02:54:05Z. This
amendment changes no prediction, test or threshold. It records what was found on downloading and preparing the
cohorts, before any signature, decoy or null set was scored and before any expression value was related to outcome.

1. Download route. The GEO FTP server is not reachable from this cluster; each series' sample data tables and each
   platform table were taken from the GEO query pages (acc.cgi, view=data), the sample annotations from the brief
   view saved at screening (external/meta/). Scripts: external/scripts/download.sh, parse_soft.py, prep_scanb.py.
2. The rule applied after mapping (Section A). STAD GSE26253 maps to 13,509 genes (fewer than 15,000) and ESCA
   GSE53625 is on a platform (GPL18109) whose GEO table carries no gene symbols: both are dropped, not replaced.
   Stomach is therefore represented by GSE84437 alone, and the primary correlation runs over 11 types (STAD, SKCM,
   LUAD, LUSC, COAD, CESC, SARC, HNSC, BLCA, OV, BRCA). Patients and events after the annotations' filters
   (external/cohorts/manifest.csv): GSE84437 431/207, GSE65904 210/102, GSE72094 398/113, GSE157011 483/264,
   GSE39582 556/187, GSE44001 300/38, GSE21050 289/101 (20 with time 0 excluded), GSE65858 253/85 (primary tumors),
   GSE13507 165/69 (primary bladder cancer), GSE32062 260/121 (GPL6480 samples), GSE96058 3,273/336.
3. Signatures measured in too few genes. Over the 11 cohorts (11,373 genes common to all), CRMA (Shukla) keeps 1 of
   its 5 genes in some cohort and PDL1 (Nishino) 2 of 3. As in the breast application (sets with fewer than three
   genes measured in some cohort excluded), both are left out of the survival application; the other 50 are scored,
   and every family mean is over the signatures that remain (inflamed 20; confounder 8; other 22). The response
   application (GSE176307) keeps every signature with at least three measured genes.
