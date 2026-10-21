"""Additional file 3 (analysis code) for the Genome Biology split.

usage: build_addfile3.py <Additional_file_5_analysis_code.zip from the GM package> <outdir>

Starts from the GM code archive, keeps the scripts and deposited plans behind what the GB
paper still reports, drops those that served only the attribution / triage / locked-prediction
/ TME-library / journal-build work, and adds the target-decoy pipeline, the GB figure and
manuscript scripts, and the producers of kept tables that the GM archive lacked. Kept tables that
also carry results of the separate paper (the locked-signature library of the I-SPY2 interaction
tables, the locked predictions of the registered cohorts) are cut to this paper's rows and
columns (FILTERS). The archive is mirrored by the public repository REPO. Optional files
(the GB Figure 2 and I-SPY2 supplementary figure scripts, the reference and Additional file 1
builders) are taken when they exist, so rerun this script after they change.
"""
import csv
import io
import sys
import tempfile
import zipfile
from pathlib import Path

SRC, OUTD = Path(sys.argv[1]), Path(sys.argv[2])
OUT = OUTD / "Additional_file_3_analysis_code.zip"
G = Path("/groups/adv2105_gp/yichen/Yi/com/gastric")
FS, ICB, PAN, GB = G / "figscripts", G / "icb", G / "pancan", G / "revision_GB" / "scripts"

# --- from the GM archive, unchanged ------------------------------------------------------
KEEP_SCRIPTS = [
    # benchmarks, common scale, dependence
    "null_matched_v2.R", "null_matched_prog.R", "build_dissociation_v3.py", "perm_inference.py",
    "v5_clusters.py", "cf1_cohort_stats.R", "cf2_fit.py", "positive2.R",
    # the three pre-specified tests of PLAN_three_tests.md
    "ta1_spread.R", "ta2_test.py", "tb1_axis.R", "tb2_test.py", "tc1_gene_arms.py", "tc2_perm.py",
    # survival under checkpoint blockade, within-arm analyses, pooled interaction, gastric random axes
    "null_matched_icbsurv.R", "analyse_icbsurv.py", "v8_within_arm_null.py", "v8b_arm_power_check.py",
    "wa1_perm.py", "pi1_pooled.py", "arm_split.py", "ispy2_within_arm_null.py", "v2_random_logrank.R",
    # figures (GM numbering in the file name; see README) and the shared style
    "fig1_v6.py", "fig3_v6.py", "fig4_v6.py", "fig5_v6.py", "fig6_v6.py",
    "supp01_scoring.py", "supp_gastric.py", "supp_gastric_cox.py", "supp_trials.py", "supp_ceiling.py",
    "figstyle.py", "v4kit.py",
]
KEEP_VALIDATION = [
    # deposited plans still behind GB analyses, with their hash files
    "PREREGISTRATION_external_validation.md", "AMENDMENT_1_before_scoring.md", "PREREGISTRATION.sha256",
    "PLAN_closed_form_null.md", "PLAN_icb_survival.md", "PLAN_icb_survival.sha256",
    "PLAN_pooled_interaction.md", "PLAN_three_tests.md", "PLAN_within_arm_inference.md",
    # the GEO search and the I-SPY2 analysis
    "geo_icb_search.csv", "ispy2_prep.py", "score_ispy2_native.R", "ispy2_interaction.py",
    "ispy2_within_arm_null.py", "ispy2_direction_posthoc.py", "ispy2_interaction_by_signature.csv",
    "ispy2_interaction_summary.csv", "ispy2_within_arm_null.csv", "ispy2_direction_posthoc.csv",
    # the seven registered cohorts, read by tests A and B
    "cohort_manifest.csv", "prep_cohorts.py", "score_validation.R", "validation_scored.csv",
    "validation_by_signature.csv",
]
REPO = "https://github.com/drgengtaiye/signature-decoy-fdr"


def _rows(keep_row=lambda r: True, drop_cols=()):
    """A CSV filter: keep the rows keep_row accepts, without the columns in drop_cols."""
    def f(text):
        rd = csv.DictReader(io.StringIO(text))
        cols = [c for c in rd.fieldnames if c not in drop_cols]
        out = io.StringIO()
        w = csv.DictWriter(out, cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(r for r in rd if keep_row(r))
        return out.getvalue()
    return f


# archive path -> filter applied to the file's text before it is written
FILTERS = {
    "validation/ispy2_interaction_by_signature.csv": _rows(lambda r: r["library"] == "trial library"),
    "validation/ispy2_interaction_summary.csv": _rows(lambda r: r["library"] == "trial library"),
    "validation/validation_by_signature.csv": _rows(drop_cols=("pred",)),
}

# --- added from the live tree --------------------------------------------------------------
ADD = {
    "decoys": [ICB / f for f in [
        "null_matched_v3_coexpr.R", "run_null_v2.sh", "agg_null_v3.R",
        "calib_decoys_v3.R", "calib_score_v3.R", "run_calib_v3.sh",
        "calib_v3_fdr.py", "calib_v3_robust.py", "calib_v3_rcurve.py",
        "null_natpool_v4.R", "calib_read_natbank.R", "run_null_natpool_v4.sh",
        "calib_ownbank_v3.R", "run_calib_v3_own.sh", "run_calib_v3_ownnat.sh", "calib_v3_ownbank_compare.py",
        "pool_shift_check.R", "export_decoy_genes_v3.R"]],
    "scripts": [ICB / "standardised_effects.py", FS / "fig_targetdecoy.py"],
    "trials": [PAN / f for f in [
        "interaction2.py", "arm_specific.py", "v4_arm_split.py", "v7_lib.py",
        "v7_interaction_null_valid.py", "v7_fig2_draws.py", "v7_calibration.py",
        "v7_calibration_grid.py", "v7_calibration_real.py", "v7_observed_check.py"]],
    "manuscript": [GB / f for f in [
        "build_background.py", "build_targetdecoy.py", "build_front_v2.py", "build_manuscript_GB.py",
        "build_manuscript_GB_v2.py", "build_manuscript_GB_v3.py",
        "audit_refs.py", "build_addfile2.py", "build_addfile3.py", "build_addfile4.py", "build_package.py"]],
}
OPTIONAL = {
    "scripts": [FS / "fig2_gb.py", *sorted(FS.glob("*ispy2*.py"))],
    "manuscript": [GB / "build_refs_GB.py", *sorted(GB.glob("build_addfile1*.py"))],
}
# without a GB Figure 2 or an I-SPY2 figure script, fall back to the GM scripts that draw them
if not (FS / "fig2_gb.py").exists():
    KEEP_SCRIPTS.append("fig2_v6.py")
if not list(FS.glob("*ispy2*.py")):
    KEEP_SCRIPTS.append("fig8_v6.py")

README = """\
Analysis code for: Matched decoy modules give immune gene-expression signatures a false discovery rate

Scripts were run on one HPC system and carry its absolute paths; set the data root at the top of
each script (and D in scripts/v4kit.py) before running elsewhere. The same files are at
{repo}.

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
6  Figures. Figure 1 fig1_v6.py; Figure 2 {fig2}; Figure 3 fig3_v6.py; Figure 4 fig_targetdecoy.py;
   Figure 5 fig4_v6.py; Figure 6 fig5_v6.py; Figure 7 fig6_v6.py (file names keep the previous
   numbering). Additional file 1: supp01_scoring.py, supp_gastric.py, supp_gastric_cox.py,
   supp_trials.py, supp_ceiling.py, {ispy2}.
7  Manuscript. manuscript/build_manuscript_GB.py rebuilds the manuscript from the previous version,
   taking its new sections from build_background.py, build_targetdecoy.py and build_front_v2.py;
   build_manuscript_GB_v2.py applies the journal format (abstract, section order, tables) and the
   text fixes of the split, build_refs_GB.py then renumbers the references, and
   build_manuscript_GB_v3.py sets the code-availability statement. build_addfile1.py,
   build_addfile2.py, build_addfile3.py and build_addfile4.py build Additional files 1-4, and
   build_package.py assembles the submission.
"""


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=OUTD) as tmp:
        old = Path(tmp) / "gm_code"          # the GM archive, extracted on its own
        old.mkdir()
        with zipfile.ZipFile(SRC) as z:
            z.extractall(old)
        files = [(old / "scripts" / f, f"scripts/{f}") for f in KEEP_SCRIPTS]
        files += [(old / "validation" / f, f"validation/{f}") for f in KEEP_VALIDATION]
        for d, paths in ADD.items():
            files += [(p, f"{d}/{p.name}") for p in paths]
        for d, paths in OPTIONAL.items():
            files += [(p, f"{d}/{p.name}") for p in paths if p.exists()]
        missing = [str(p) for p, _ in files if not p.exists()]
        assert not missing, f"missing: {missing}"
        names = [a for _, a in files]
        assert len(names) == len(set(names)), "duplicate archive path"
        fig2 = "fig2_gb.py" if (FS / "fig2_gb.py").exists() else "fig2_v6.py (panel e and two rings of a belong to the separate paper)"
        isp = ", ".join(p.name for p in FS.glob("*ispy2*.py")) or "fig8_v6.py panels g-i (I-SPY2)"
        tmp_out = OUT.with_suffix(".tmp")
        with zipfile.ZipFile(tmp_out, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("README.txt", README.format(fig2=fig2, ispy2=isp, repo=REPO))
            for p, a in sorted(files, key=lambda x: x[1]):
                if a in FILTERS:
                    z.writestr(a, FILTERS[a](p.read_text(encoding="utf-8")))
                else:
                    z.write(p, a)
        tmp_out.replace(OUT)
    print(f"{OUT}: {len(files) + 1} files")


if __name__ == "__main__":
    main()
