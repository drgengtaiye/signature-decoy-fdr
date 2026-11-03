# Amendment 1 to PLAN_round15_literature_cohorts_controls.md: the screen, the rules fixed during it, and the frozen table

Written 2026-11-03, after the screening of part A and before any new model, decoy, cohort or assay was scored or
downloaded for parts A-D (only gene identifiers of the TCGA and round-11 cohorts were read, to map symbols, as in round
13). Its SHA-256 is recorded in AMENDMENT_1_round15_literature_cohorts_controls.sha256 and it is deposited before any new
model is scored.

## The screen (A1)

The 1,598 records of the round-13 search not screened in round 13 (STAD 183, SKCM 104, LUAD 323, LUSC 40, COAD 300,
CESC 60, HNSC 122, BLCA 169, OV 83, BRCA 214; SARC had none left) were fetched from PMC (none failed) and reduced by
script to excerpts (round15_search/extract_snippets_r15.py: title, abstract, accessions, windows around printed
products of a coefficient and a gene, coefficient tables). On the 109 round-13 eligible models every gene of 103 lay in
the excerpts; readers searched the full text whenever the excerpts did not settle a record. Ten readers screened every
record in the frozen order under round15_search/CURATION_BRIEF_r15.md, one row each (coverage checked by script). One
reader found on a re-scan that formulas typeset in TeX had been missed and corrected four rows; a script then searched
the full text, MathML and TeX included, of all 351 records rejected for want of a usable formula, and none of the 13 it
flagged was a missed formula (round15_search/recall_scan_r15.csv). Readers ruled 629 models eligible. A further reader
recorded, for the 294 round-13 records, whether each trains a TCGA model and the cohorts it validates in (for part B;
round15_search/r13_backfill.csv); round-13 decisions are unchanged.

## Second extraction (A2)

A random 63 of the 629 (numpy default_rng(20261201) over the models sorted by type in the plan's order and then rank;
round15_search/second_pass_sample_r15.csv) were re-extracted by a reader who did not see the first extraction. All 63
agreed in genes and coefficients; the only difference was the printed spelling of one gene (IL-1A recorded as IL1A).

## Rules fixed during screening (applied to every new record; round 15 rule numbers)

1. Rules 1-3 of amendment 2 of round 13 stand.
2. Protein-coding means NCBI gene type "protein-coding" in org.Hs.eg.db (EGSOURCEDATE 2025-Sep24), besides the
   immunoglobulin and T-cell-receptor segments of round 13. Pseudogenes (CYP4Z2P, PSMD10P1, GBP1P1) and genes now typed
   ncRNA (SFTA3; DIRC1, SPATA8) make a model ineligible. Readthrough genes typed protein-coding (JMJD7-PLA2G4B) and
   IGLL5 count as protein-coding. UGT1A, a complex locus of protein-coding genes in a round-13 model, stands.
3. Printed numbers: a comma between digits is a digit-group separator ("0.105,367" is 0.105367); "× 10^-n", "× 10−n"
   and "e−n" are powers of ten; a minus operator printed before a negative coefficient is not applied twice (the signed
   coefficient as printed); when the text cites a table and contradicts it, the table is used (round 13); a constant
   added to the score is left out (ranks do not change).
4. The main model is the one the title or abstract presents as the paper's signature; when they present none or several,
   the first printed survival model. Re-read for every row whose notes record a choice between models: 11 of 13 stood;
   LUAD PMC9302251 (main model the 7-gene GEO-cohort signature: not trained on TCGA) and HNSC PMC9948032 (main model a
   composite of two within-cohort z-scored scores: not a fixed linear formula) became ineligible.
5. Coefficients must come from a fit to a survival endpoint: a score fitted to a binary endpoint (diagnosis,
   metastasis) is "no prognostic model" even if later used for survival (LUAD PMC8182055, PMC12686163 became ineligible).
6. A printed coefficient with more than one reading makes the model ineligible (COAD PMC9513614, DLAT "08,399"). One
   with a single literal reading is used as printed (OV PMC12167356, GRIP2 11.80464867091524).
7. A paper enters once. When the same model is eligible under two types (a training set pooling two TCGA types), it
   enters under the type where it was first found eligible: round 13 before round 15, then the plan's order of types
   (LUAD PMC6953332 stays round-13 LUSC; PMC12321723 enters under LUAD, not LUSC).
8. Weights that are the paper's printed hazard ratios are scored as printed (the claim audited is the printed formula),
   and the two such models (STAD PMC7982000, CESC PMC11319947; found by checking every model whose weights are all
   positive and span 1) are left out in a sensitivity analysis. Other unusual signs are used as printed.
9. Printed names read as the paper itself describes them where the alias is ambiguous or absent: PD-1 PDCD1, NOV CCN3,
   HSP27 HSPB1, LOC158830 CXorf65; misprinted symbols corrected by a reader from the same paper's text are listed in each
   row's notes and in round15_search/judgment_calls_r15.md.
10. Near-duplicate publications (COAD PMC8847155, PMC7969881) both enter: the audit samples publications.
11. Scores with non-expression terms (age, sex, stage) or genes coded 0/1 stay ineligible, as the readers ruled.
12. For part B, `validation_used` leaves out internal TCGA splits and series merged into training; readers differed on
    checkpoint-treated cohorts, which B2 excludes in any case.

These rules changed 11 rows (round15_search/apply_rules_r15.py; screen_r15_final.csv keeps each reader's decision
beside the one that stands): 619 new models are eligible.

## The frozen table (A2-A3)

validation/published_models_r15.csv (applications/v20_freeze_models.R): the 109 round-13 rows unchanged and 619 new
ones parsed and mapped as in round 13 (amendment 2, items 4-6): 728 models, 707 scored (21 below the 80% rule), 678 with
three or more measured genes. The primary population of A4 (scored, three or more genes, paper naming neither the
round-11 cohort of its type nor a related series anywhere in its text) holds 427 models (60 from round 13, 367 new):
STAD 38, SKCM 11, LUAD 96, LUSC 20, COAD 58, CESC 20, SARC 4, HNSC 27, BLCA 33, OV 37, BRCA 83. Column hr_weights marks
the two models of rule 8.

## Frozen files (SHA-256)

    925469fa86a68975bc22fc09d434247b5e5b6f7da1145d63f33fd43e10675b5c  published_models_r15.csv
    24401d48bd7856115e35c69b7586794b86b8933e954f1cb753834525a5052853  round15_search/screen_r15_final.csv
    b4359a8a6123bfc6e9c40821b921f18ac03c2854e579eac2bf578dbfd23bb736  round15_search/screen_r15_all.csv
    2d3507775abfafc5fd202afea34b5af8b4bbb390b9fa4b442925669178f40c76  round15_search/r13_backfill.csv
    370c9fff2bb9dd9ff565d659222a1a1980cef5d3ff54c462f79263e73615dc1f  round15_search/second_pass_sample_r15.csv
    cc04f6b61ab145e082c80252fbcf0d530e3d366075e992089fd041eb9f544b03  round15_search/second_pass_r15.csv
    5f746110067315f49420c48fa5e54b462f23bc5ec3fa1ab0b36397c70c1765e3  round15_search/second_pass_agreement_r15.csv

The per-chunk screening files, the recall scan, the judgment-call log, the scripts and the brief are listed with their
SHA-256 in round15_search/SHA256SUMS_freeze_r15.
