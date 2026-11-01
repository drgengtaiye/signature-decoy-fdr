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
7  Figures (validation/PLAN_round8_figure_displays.md and its four amendments: the descriptive displays added last;
   Figure 5c, j and Additional file 1 Fig. S2: validation/PLAN_round9_standard_fdr_comparison.md; Figure 8:
   validation/PLAN_round10_content_and_structure.md, and panels f and g validation/PLAN_round11_external_replication.md).
   figprep_example.py, figprep_r.R, figprep_roc.py and figprep_km.R write the data of the example panels
   (icb/v13_fig/); cns.py holds the shared style and cns10.py its palette and the forms of the last redraw. Figure 1
   scripts/fig1.py; Figure 2 fig2.py; Figure 3 fig3.py; Figure 4 fig4.py; Figure 5 fig5.py; Figure 6 fig6.py; Figure 7
   fig7.py; Figure 8 fig8.py. Additional file 1: python figS2.py; python run_restyled.py figS1.py (which also draws the figure set's
   earlier Figure 4, not used) and figS3.py to figS7.py, whose scripts save under their earlier names (FigS2_calib to
   FigS6_scale); python figS8.py (round 10: single-cohort estimates and libraries); run_restyled.py applies the palette of cns10.py (fig67_lib.py, figstyle.py and v4kit.py are their
   shared panels and style).
