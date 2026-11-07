# Result of the acceptance test of `check_models()` (decoyfdr 0.3.0)

Run 7 November 2026 against the criteria of TOOL_check_models_acceptance.md (SHA-256 d56647c2..., hashed before the
run). Outputs in icb/v24_tool (repro_units/scores/decoys_<TYPE>.csv.gz, repro_report.csv, logs/). Package installed in
icb/v24_tool/Rlib; the frozen install icb/v12_pkg/Rlib6 was not touched. Run time of `check_models()` per cancer type
(all its frozen models, decoy growth, every built cohort of the type), 6 cores: 31 s (SARC) to 1.7 min (LUAD).

| Criterion | Checks | Passed |
|---|---|---|
| 1. Decoys: names and order, genes, coefficients (8 significant digits), directions | 44 | 44 |
| 2. Part A in TCGA and the external cohort: C-index (1e-10), log-rank P (1e-8 relative), high group worse, validated | 99 | 99 |
| 3. Per model: d_coex, d_rand, n_coex, rho, validated in TCGA | 77 | 77 |
| 4. Estimates: validated count, decoy share, FDR, rho for the 427 independent models and all 678 with decoys | 22 | 22 |
| 5. Part B5 (142,570 model and decoy rows in 41 cohorts): rows, scorable, validated, C-index (1e-6) | 45 | 44 |

Criterion 4 in numbers: independent models, 180 of 427 validated, FDR 0.3951389 with co-expressed decoys (round 15:
0.395139); all models with decoys, 306 of 678, FDR 0.3184232 (0.318423). `validation_fdr()` of the package gives the
same point estimates; its bootstrap interval within cancer type (R's generator, 2,000 draws) is 35.4-44.2% against the
paper's 35.7-44.2% (numpy, 10,000 draws), as the criteria anticipated.

## The one failed check

Criterion 5 failed on one row of 142,570: decoy dcy_coex_BRCA_PMC11537572_11 in BRCA_GSE96058, directed C-index
0.4716990062 from the function and 0.4716978230 in part B, a difference of 1.18e-6 against the tolerance of 1e-6.
Its validation (not validated, log-rank P 0.3187) is the same in both. Cause, checked: part B read the decoys'
coefficients rounded to 8 significant digits (the criteria name this as the reason for the tolerance); 14 of this
decoy's 18 coefficients change by up to 4.4e-8 when rounded, and rescoring it with the rounded coefficients gives
0.4716978230, part B's value to ten digits. In GSE96058 (3,273 patients) one pair of patients whose scores the rounding
reorders moves the C-index by about 1e-6, so the tolerance set before the run was too tight for the largest cohort.
Six rows in all differ by more than 1e-9 (all decoys, all below 1.2e-6, none with a different validation).

The known exception stated before the run behaved as stated: BRCA_PMC11332736 in BRCA_GSE96058 is scored on 7 genes
(part A) by the function and was scored on 8 in part B (C-index 0.5620 against 0.5461; validated in both).

`validated` agrees in every row of parts A and B5. The function is released with these two differences documented;
neither changes a number of the paper.
