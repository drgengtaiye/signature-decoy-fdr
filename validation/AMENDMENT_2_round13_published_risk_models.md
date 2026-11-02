# Amendment 2 to PLAN_round13_published_risk_models.md: the frozen models and the rules of implementation

Written 2026-11-02, after the screening of part A and before any published model, theme pool or decoy was scored in any
cohort (one test of v19_pipeline.R ran the theme pools themselves, without decoys, through the pipeline in CESC to
check the code: 15 pools, one validated; nothing else). Its SHA-256 is recorded in
AMENDMENT_2_round13_published_risk_models.sha256 and it is deposited before any model is scored.

## The screen

Under Amendment 1 the search found 1,892 records (STAD 220, SKCM 136, LUAD 353, LUSC 68, COAD 329, CESC 79, SARC 17,
HNSC 140, BLCA 201, OV 103, BRCA 246). Screening in the A2 order read 294 records and found 109 eligible models: 10 per
type except SARC, 9 of whose 17 records were eligible. Every record, eligible or not, with its reason, is in
round13_search/screen_<TYPE>.csv.

## Rules fixed during screening (applied to every type)

1. Immunoglobulin and T-cell-receptor gene segments (IGHV, IGKV, IGLV, IGKC, TRAV, TRDC and the like) are not
   protein-coding (their biotypes are IG_ and TR_ genes) for A3(iii). Four models are ineligible for this reason
   (COAD ranks 3 and 9, LUSC ranks 7 and 25; BRCA rank 31 was so ruled from the start) and screening continued.
2. A training set that is a subset of the TCGA type (muscle-invasive bladder, triple-negative breast), or the type
   merged with other data (TCGA-OV with two GEO series; TCGA LUAD with LUSC for lung squamous), counts as TCGA tumors
   of the type. Gene-pair models, hazard ratios without coefficients, scores of principal components, clusters,
   splicing events or machine-learning models are not linear formulas of gene expression. A score printed as the
   exponential of a linear predictor is taken as the linear predictor (same ranks). When text and a table both give
   the coefficients, the text is used unless it cites the table and contradicts it.
3. Coefficients are parsed with every minus-like character (U+2212, U+2013, U+2014, U+02D7 and the like) read as a
   minus and brackets removed; obvious print splits inside symbols ("KLF 9", "ATG1 0") are joined.

## Second extraction (A4)

A random 22 of the 109 eligible models (numpy default_rng(20261114); round13_downloads/second_pass_sample.csv) were
re-extracted from the full text by a second reader who did not see the first extraction. Twenty agreed in genes and
coefficients; two did not, and both were resolved against the full text in favor of the second reading: OV PMC9050317
prints CXCL11 as -0.161 in the formula but cites Table 1, which gives -0.01606 (the other coefficients are Table 1's,
rounded), so Table 1's value is used; BRCA PMC7917807 prints its last gene as "ACSF", which is not a symbol, where the
same paragraph lists ACSF2. Agreement before resolution: 20 of 22 (91%) in genes and coefficients, 21 of 22 in genes.

## Implementation of A5 to A8 and B2 to B3

4. Genes (A5): a printed symbol measured in both cohorts of the type is kept; otherwise it is mapped through
   org.Hs.eg.db (official symbols and aliases to Entrez Gene) to the one symbol of that gene measured in both. A gene
   measured under different symbols in TCGA and the external cohort counts as unmeasured. With the 80% rule, 106 of the
   109 models are scored (16 genes of 14 models unmeasured); three fall below 80% (STAD PMC11276922, HNSC PMC9304892,
   BLCA PMC11625674).
5. Population (A10): a model's paper uses the external cohort when it names, anywhere in its text, the external
   cohort's accession or a GEO series that contains it or that it contains (STAD GSE84426, GSE84433, GSE147163; LUSC
   GSE157009, GSE157010; COAD GSE40967; BRCA GSE81540): broader than "among its validation cohorts", so that models
   trained or selected with the cohort are also left out. 63 scored models are independent in this sense, 60 of them
   with three or more genes, so the primary population is the independent models.
6. decoyfdr grows decoys only for sets of three or more genes; the five scored two-gene models (SKCM PMC9715341, LUSC
   PMC10258956, OV PMC10837903, BRCA PMC9855618, BRCA PMC7877483) are scored and reported but have no decoys and do
   not enter the estimates.
7. A decoy's coefficients are its model's in random order, drawn with Mersenne-Twister seed 20261115 + type index,
   decoy by decoy in the order grow_decoys returns them; its direction is that of its score's C-index in TCGA (C at or
   below one half: higher score, shorter survival).
8. The median split puts scores equal to the median in the lower group. In B3 the fold ids of decoy pool j (1-20
   co-expressed, 21-25 random) of theme pool p in type t use seed 20261117 + 100000 + 10000 t + 100 p + j.

## Frozen files (SHA-256)

    3dea385e5b652e513554f27b1680b32b05199548a464182b2beef760bac3355b  published_models.csv
    f70ac1e592710c7a109ba83c6458e839975febc05d042b2e293bbe236f4f13a7  round13_search/screen_BLCA.csv
    b84e0a9d2344e7904b8b0e6171c014deb7353eaac9f35e82cf9ad003e89da228  round13_search/screen_BRCA.csv
    7a5235e189ba84ae2021f64eef231199874ce24b84ae495fe5009f14c353ac4d  round13_search/screen_CESC.csv
    fdc9d4c935af0eae847745a578b7dc5acf2745b03601b203a23cbdf6e3557a2f  round13_search/screen_COAD.csv
    9883a4cc04ba0edfffdca166e34bddbb61b953333afec50044f6c9b04217e1aa  round13_search/screen_HNSC.csv
    56505c4a2a066ddbbdafed272cedc7cf18b139a6928a3255649e8023bf33926b  round13_search/screen_LUAD.csv
    874935734ad9bdbc917f033cc09ad2d4bb56757770b59df77478fbb89ad55b82  round13_search/screen_LUSC.csv
    f391e8f4d2c0de6605d658357abc7e181fbfc68088e97655861c73f80e6c6379  round13_search/screen_OV.csv
    756b6125e84f6cfdbec41aab66663fab89ca214ce7e5601aefa318ee46f211e7  round13_search/screen_SARC.csv
    7100f5e5e50b9f0cf84ed00ddf50bf769fd3741c0198dfd629b5602c43144ad1  round13_search/screen_SKCM.csv
    933fcc28503f2911b6d3c4964eea1137dfcea5f3ec0871c54fb8830024a236dd  round13_search/screen_STAD.csv
    fee3017ae3a7884a1851fde62d1e0f4a00c52337083ada0d524103f0994698ea  round13_search/second_pass_agreement.csv
