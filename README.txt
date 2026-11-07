Analysis code for: Matched decoy modules estimate the false discovery rate of claims about gene-expression
signatures

The R package decoyfdr (decoyfdr/) implements the procedure; install it with
remotes::install_github("drgengtaiye/signature-decoy-fdr", subdir = "decoyfdr"). Its tests reproduce the
response-benchmark calls and estimates from the per-cohort statistics shipped in decoyfdr/inst/extdata.

The analysis scripts were run on one HPC system and carry its absolute paths. INPUTS.csv lists every input file
with its source; lay them out under a local root and run  python scripts/set_root.py /path/to/root  to rewrite
the paths in every script. ENVIRONMENT.txt gives the R and Python versions used. Seeds in the scripts are arbitrary
labels; several have the form of calendar dates.

decoyfdr/      the R package (R code, C++ for the C-index, tests, vignette)
applications/  the package run on the response cohorts, on 20 TCGA cancer types, on breast cancer cohorts and on
               the MSigDB Hallmark library
breast/        download and preparation of the breast cancer cohorts and of the Venet et al. signature library
external/      screening, download and preparation of the independent survival cohorts and of GSE176307
scripts/       benchmarks, dependence between signatures, figures
decoys/        the target-decoy estimate on the response benchmark and its calibration
validation/    analysis plans with their SHA-256 hashes; provenance of the immune and breast signatures

Run order (each step reads the outputs of the steps before it)
1  Benchmarks. scripts/null_matched_v2.R: response AUC in each checkpoint cohort with its size- and engine-matched
   null (native_response_matched_v2.csv). scripts/null_matched_prog.R: C-index in 20 TCGA cancer types
   (native_prognosis_matched.csv). scripts/build_dissociation_v3.py: per-signature effects (v3_dissociation.csv).
   scripts/standardised_effects.py: effects in null standard deviations (v4_standardised_effects.csv).
   scripts/perm_inference.py: label-permutation inference and effective number of tests
   (v4_permutation_inference.csv). scripts/v5_clusters.py: correlation clusters (v5_signature_clusters.csv).
   scripts/cf1_cohort_stats.R, cf2_fit.py: the closed-form null (validation/PLAN_closed_form_null.md).
   scripts/positive2.R: consensus, supervised model and single signatures. Signed null
   (validation/PLAN_calibration_round5.md): scripts/v10_signed_null.R draws signed nulls for the four weighted-mean
   signatures with weights of both signs and v10_apply_signed.py puts them into the benchmark tables (the unsigned
   rung kept as ladder_*); scripts/score_matrices.R writes the score matrices that perm_inference.py and
   v5_perm_family.py read, and v10_signed_ranks.py adds the signed sets to them; then v10_dependence_aware.py,
   v10_label_vs_matched.py and v10_prog_cluster_boot.py.
2  Target-decoy estimate on the response benchmark (first decoy set). decoys/run_null_v2.sh runs
   null_matched_v3_coexpr.R (1,000-draw N0/N1/N2 banks per cohort), then agg_null_v3.R (null_v3_tidy.csv).
   calib_decoys_v3.R grows the 2,912 decoys; run_calib_v3.sh runs calib_score_v3.R per cohort; calib_v3_fdr.py
   estimates the FDR (calib_v3_units.csv, calib_v3_fdr.csv, calib_v3_per_signature.csv); calib_v3_robust.py gives
   the paired and matching-sensitivity estimates; calib_v3_rcurve.py the call rate by decoy co-expression.
   run_null_natpool_v4.sh builds the natural-pool N1/N2 banks (null_natpool_v4.R) and reads the decoys against them
   (calib_read_natbank.R). run_calib_v3_own.sh / run_calib_v3_ownnat.sh and calib_v3_ownbank_compare.py check the
   shared-bank shortcut; pool_shift_check.R tests the draw pool; export_decoy_genes_v3.R writes the decoy gene lists.
3  Calibration (validation/PLAN_calibration_round5.md). decoys/v10_export_n0.R and v10_null_decomp.py compare the
   nulls of Figs 3 and 4; run_v10_store.sh runs v10_store.R per cohort (ranks of every target, decoy and N0 set,
   reproducing the published banks), then v10_store_json.R; v10_calib_sim.py repro | null | signal reruns the whole
   procedure on the observed labels, on 1,000 label permutations and on planted signal (v10_calib_lib.py holds the
   engine), v10_calib_signal_totals.py summarizes the latter; v10_fdr_intervals.py gives the dependence-aware
   intervals and the FDR with each signature's origin cohorts left out (validation/signature_provenance.csv);
   v10_cohort_boot_diag.py checks the cohort bootstrap.
4  Package and applications (validation/PLAN_round6_applications.md, AMENDMENT_2_round6_structural_truth.md).
   R CMD INSTALL decoyfdr. applications/common.R holds the loaders. applications/v11_response_rerun.R reruns the
   response benchmark with the package; v11_tcga_prognosis.R and v11_tcga_resume.R the prognostic application;
   breast/fetch_breast_data.sh, export_orgHs_symbols.R, export_venet_signatures.py, prep_breast_cohorts.py and
   verify_breast_cohorts.py prepare the breast cohorts and library, and applications/v11_breast.R,
   v11_breast_pc1_strength.R and v11_breast_absolute.R are the breast application; v11_post.R and
   v11_decoy_r_by_cohort.R give secondary readings. decoys/v11_structural_truth.py reads the planted-signal
   replicates against the expected-delta truth.
5  The estimand, the bias and ten decoy sets (validation/PLAN_round7_estimand_and_bias.md; decoyfdr 0.2.0).
   applications/v12_rerun.R response | tcga | breast | breast_abs reruns every application with ten independent
   decoy sets (set 1 reproduces step 4 exactly; the check is printed); applications/v12_breast_origin.R leaves out
   each breast signature's origin cohorts (validation/venet_provenance.csv, frozen with its SHA-256 before it was
   applied); applications/v12_hallmark.R measures the decoy pool at library scale and runs the MSigDB Hallmark
   library. decoys/v12_call_truth.py signal <0-7> | null reads the response planted signal against the
   call-probability truth and adds the pi0-corrected estimate and the coverage of the intervals;
   decoys/v12_tcga_planted.R <0-5> plants survival signal with new causal gene sets in TCGA; v12_truth_summary.py
   summarizes both. decoys/v14_standard_fdr.py sets the decoy estimate beside the nominal (Benjamini-Hochberg,
   Storey) and permutation (SAM-type) estimates of the same calls, from the stored outputs of this step
   (validation/PLAN_round9_standard_fdr_comparison.md; icb/v14_standard/). decoys/v15_content_structure.py splits each
   signature's gain into what its co-expressed decoys achieve (structure) and what its genes add (content), by endpoint,
   family and cancer type, from the same stored outputs (validation/PLAN_round10_content_and_structure.md;
   icb/v15_content/).
6  Independent cohorts (validation/PLAN_round11_external_replication.md and AMENDMENT_1, deposited before the
   downloads and before the scoring). external/screen_meta2.py reads the candidates' sample annotations
   (list_a.txt, list_b.txt); external/download.sh fetches each series' data tables and platform tables through the
   GEO query page; external/parse_soft.py maps probes to symbols and writes each cohort's expression and survival
   from the fields named in cohorts.json; prep_scanb.py takes GSE96058 from the breast preparation;
   prep_validation_cohorts.py is the September 2026 preparation of GSE176307. applications/v16_external.R survival |
   response (run_v16.sbatch) runs the package on the 11 survival cohorts and on GSE176307 (icb/v16_external/);
   decoys/v17_external_replication.py computes the plan's tests from them (icb/v17_external/).
   Split halves (validation/PLAN_round12_reproducibility_by_type.md, deposited before it was computed, and three
   amendments: a tie rule, one round-10 cell accepted, and a post hoc bootstrap; AMENDMENT_2's first deposited version
   is kept as _superseded). applications/v18_split_half.R tcga | external | tcga_boot | external_boot (run_v18.sbatch)
   splits every TCGA type and independent survival cohort 500 times into halves stratified by event status and
   bootstraps its patients 200 times, from the per-patient scores of the step-5 TCGA fit and of the independent cohorts
   (recomputed by the steps of fit_units and checked against step 6); decoys/v18_reproducibility.py computes the
   reliabilities, the discordance test, the power simulation and the single-signature rates (icb/v18_split/).
7  Published prognostic risk models (validation/PLAN_round13_published_risk_models.md with two amendments, deposited
   before any paper was retrieved and before any model was scored). validation/round13_search/round13_search_amended.py
   runs the PubMed search and orders the records (search_summary.json, <TYPE>_order.tsv); the screening tables
   screen_<TYPE>.csv record every record read, with the formula extracted, and second_pass_*.csv the independent
   re-extraction of a random 22 (compare_second_pass.py). applications/v19_freeze_models.R maps the eligible models'
   genes to the cohorts' symbols and writes validation/published_models.csv (frozen; SHA-256 in AMENDMENT_2).
   applications/v19_published_models.R <TYPE> scores each model and its decoys in TCGA and the independent cohort of
   the type; applications/v19_pipeline.R <TYPE> runs the published pipeline (univariate Cox, LASSO-Cox) on 15 theme
   gene pools and their decoy pools (run_v19.sbatch; v19_common.R holds the shared pieces); decoys/v19_audit.py gives
   the estimates (icb/v19_audit/).
8  Checkpoint cohorts outside the collection (validation/PLAN_round14_new_checkpoint_cohorts.md, deposited before any
   sample annotation was read): validation/round14_search/ holds the GEO search (query.txt, gds_ids.txt,
   geo_filter.py, shortlist.csv) and the screen of every candidate (screen.csv); no cohort was eligible and the tests
   were not run (validation/RESULT_round14_new_checkpoint_cohorts.md).
9  Published models, validation cohorts and clinical assays (validation/PLAN_round15_literature_cohorts_controls.md and
   amendments 1-4, each deposited before the step it changes). Part A, the whole search: validation/round15_search/
   (screen, second extraction, rules), applications/v20_freeze_models.R (validation/published_models_r15.csv),
   applications/v20_published_models.R <TYPE> and decoys/v20_audit.py (estimates), decoys/v20_citations.py
   (exploratory citations). Parts B-C, validation cohorts: validation/round15_cohorts/ (candidates, eligibility, GEO
   metadata in meta_summary/, built cohorts), applications/v20_null_models.R, v20_cohort_scores.R, v20_law.R and
   decoys/v20_cohorts_audit.py (false-validation rates, H4, H5, the law). Part D, clinical assays:
   validation/round15_assays/ (genes, weights and stated approximations; RUNLOG_v20_assays.md), applications/
   v20_assays.R and decoys/v20_assays_audit.py (H7 and sensitivities). Outputs: icb/v20_round15/. A plan for the
   machine-learning pipelines of the same literature (validation/PLAN_round16_ml_pipeline.md) is deposited and was not run.
10 Decoy size (validation/PLAN_round17_decoy_size.md and AMENDMENT_1_round17_decoy_size.md): applications/
   v22_size_check.R (icb/v22_size/).
11 The law out of sample (validation/PLAN_round18_law_out_of_sample.md): applications/v23_law_split.R and
   decoys/v23_law_split_audit.py (icb/v23_law_split/).
12 Figures (validation/PLAN_round8_figure_displays.md and its four amendments; the panels added in the redraw are
   descriptive under validation/NOTE_figures_v17_descriptive_only.md, except those that show rounds 15 and 18).
   scripts/pixl_nc.py holds the style and pixl_nc_check.py the checks. Figure n is scripts/fig<n>.py: Figures 2-7
   are the content and replication analyses (rounds 5-12), Figure 8 the published models, clinical assays and
   pipeline (rounds 13 and 15), Figure 9 the validation cohorts and the law (rounds 15 and 18). Supplementary Figure
   n is scripts/figS<n>.py: Supplementary Figure 1 the library's correlation structure, 2-10 the earlier
   supplementary figures (python figS3.py; python run_restyled.py figS2.py and figS4.py to figS8.py, whose scripts
   save under their earlier names; python figS9.py, round 10; python figS10.py, round 12), 11 the first-pass sample of
   published models (round 13). figprep_example.py, figprep_r.R, figprep_roc.py and figprep_km.R write the data of
   the example panels (icb/v13_fig/); cns.py, cns10.py, run_restyled.py, fig67_lib.py, figstyle.py and v4kit.py are
   the style of the earlier supplementary figures.
13 The model check of the package (decoyfdr 0.3.0: check_models(), validation_fdr(), cohort_fvr()). Acceptance
   criteria written and hashed before the run (validation/TOOL_check_models_acceptance.md); applications/
   v24_check_models_repro.R <TYPE> reruns round 15 part A and the published models of part B5 through the function
   (icb/v24_tool/), decoys/v24_check_models_repro.py compares the output with round 15, and the result is
   validation/TOOL_check_models_result.md. decoys/v24_cohort_table.py writes the cohort table the package ships
   (decoyfdr/inst/extdata/validation_cohorts_r15.csv) and decoys/v24_cohort_page.py the searchable page
   docs/validation_cohorts.html from it.
