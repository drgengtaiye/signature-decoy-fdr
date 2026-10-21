Analysis code for: Matched decoy modules give immune gene-expression signatures a false discovery rate

Scripts were run on one HPC system and carry its absolute paths; set the data root at the top of
each script (and D in scripts/v4kit.py) before running elsewhere. The same files are at
https://github.com/drgengtaiye/signature-decoy-fdr.

scripts/     benchmarks, dependence, randomized-trial analyses and figures
decoys/      the target-decoy calibration
trials/      treatment-interaction models and the no-interaction simulation for IMmotion150 and PRINCE
validation/  deposited analysis plans with their SHA-256 hashes; the GEO search, the I-SPY2 analysis,
             and the seven registered cohorts used by pre-specified tests A and B
manuscript/  build the manuscript and Additional files 2-4 from the previous version of the paper

Run order (each step reads the outputs of the steps before it)
1  Benchmarks. scripts/null_matched_v2.R: response AUC in each checkpoint cohort with its size- and
   engine-matched null (native_response_matched_v2.csv). scripts/null_matched_prog.R: C-index in 20 TCGA
   cancer types (native_prognosis_matched.csv). scripts/build_dissociation_v3.py: per-signature effects
   (v3_dissociation.csv). scripts/standardised_effects.py: effects in null standard deviations
   (v4_standardised_effects.csv). scripts/perm_inference.py: label-permutation inference and effective
   number of tests (v4_permutation_inference.csv). scripts/v5_clusters.py: correlation clusters
   (v5_signature_clusters.csv). scripts/cf1_cohort_stats.R, cf2_fit.py: the closed-form null
   (PLAN_closed_form_null.md). scripts/positive2.R: consensus, supervised model and single signatures.
2  Target-decoy calibration. decoys/run_null_v2.sh runs null_matched_v3_coexpr.R (1,000-draw N0/N1/N2
   banks per cohort), then agg_null_v3.R (null_v3_tidy.csv). calib_decoys_v3.R grows the 2,912 decoys;
   run_calib_v3.sh runs calib_score_v3.R per cohort; calib_v3_fdr.py estimates the FDR
   (calib_v3_units.csv, calib_v3_fdr.csv, calib_v3_per_signature.csv); calib_v3_robust.py gives the paired
   and matching-sensitivity estimates; calib_v3_rcurve.py the call rate by decoy co-expression.
   run_null_natpool_v4.sh builds the natural-pool N1/N2 banks (null_natpool_v4.R) and reads the decoys
   against them (calib_read_natbank.R). run_calib_v3_own.sh / run_calib_v3_ownnat.sh and
   calib_v3_ownbank_compare.py check the shared-bank shortcut; pool_shift_check.R tests the draw pool;
   export_decoy_genes_v3.R writes the decoy gene lists for Additional file 4.
3  Randomized trials. trials/interaction2.py (interaction_all_trials.csv), arm_specific.py,
   v4_arm_split.py; trials/v7_*.py the no-interaction simulation and its calibration (v7_lib.py is shared);
   scripts/v8_within_arm_null.py, v8b_arm_power_check.py, wa1_perm.py (within each arm), pi1_pooled.py
   (pooled interaction), arm_split.py (the two atezolizumab arms), tc1_gene_arms.py, tc2_perm.py (gene
   level). I-SPY2: validation/ispy2_prep.py, score_ispy2_native.R, ispy2_interaction.py,
   ispy2_within_arm_null.py (ispy2_interaction.py also fits a second library that belongs to a
   separate paper; its tables here keep the trial library only).
4  Survival. scripts/null_matched_icbsurv.R, analyse_icbsurv.py (under checkpoint blockade,
   PLAN_icb_survival.md); scripts/v2_random_logrank.R (random axes in three gastric cohorts).
5  Tests A and B (PLAN_three_tests.md). validation/prep_cohorts.py and score_validation.R prepare and
   score the seven registered cohorts (validation_scored.csv); validation_by_signature.csv lists them
   with each signature's effect (written by an analysis that belongs to a separate paper and is not
   included). scripts/ta1_spread.R, ta2_test.py; tb1_axis.R, tb2_test.py.
6  Figures. Figure 1 fig1_v6.py; Figure 2 fig2_gb.py; Figure 3 fig3_v6.py; Figure 4 fig_targetdecoy.py;
   Figure 5 fig4_v6.py; Figure 6 fig5_v6.py; Figure 7 fig6_v6.py (file names keep the previous
   numbering). Additional file 1: supp01_scoring.py, supp_gastric.py, supp_gastric_cox.py,
   supp_trials.py, supp_ceiling.py, figS_ispy2.py.
7  Manuscript. manuscript/build_manuscript_GB.py rebuilds the manuscript from the previous version,
   taking its new sections from build_background.py, build_targetdecoy.py and build_front_v2.py;
   build_manuscript_GB_v2.py applies the journal format (abstract, section order, tables) and the
   text fixes of the split, build_refs_GB.py then renumbers the references, and
   build_manuscript_GB_v3.py sets the code-availability statement. build_addfile1.py,
   build_addfile2.py, build_addfile3.py and build_addfile4.py build Additional files 1-4, and
   build_package.py assembles the submission.
