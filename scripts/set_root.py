"""Point every script at a local copy of the data.

    python scripts/set_root.py /path/to/data_root [--dry-run]

The scripts name their inputs and outputs under one root, written /path/to/analysis_root, under which the analysis
keeps icb/ (checkpoint cohorts, null banks, score stores, package runs and most tables), pancan/ (TCGA), breast/ (the
breast cancer cohorts and library), external/ (the independent survival cohorts), figs_pub/ (figures) and this
repository as signature-decoy-fdr/. Lay out the inputs listed in INPUTS.csv under a local root with the same
sub-directories and run this once: it writes the local root into every .R, .py, .sh, .sbatch and .csv file of the
repository that names the placeholder (the .csv files are the tables of built cohorts, which give each cohort's
directory). The figure scripts also need the four Helvetica faces (Helvetica-Regular.ttf, -Bold, -Oblique,
-BoldOblique) in the directory named by the environment variable HELVETICA_FACES.
"""
import pathlib, sys

OLD = "/path/to/analysis_root"
new = sys.argv[1].rstrip("/")
dry = "--dry-run" in sys.argv
repo = pathlib.Path(__file__).resolve().parents[1]
n = 0
for p in sorted(repo.rglob("*")):
    if p.suffix in (".R", ".py", ".sh", ".sbatch", ".csv") and p.name != "set_root.py" and ".git" not in p.parts:
        s = p.read_text(errors="surrogateescape")
        if OLD in s:
            n += 1
            print(f"{p.relative_to(repo)}: {s.count(OLD)} path(s)")
            if not dry:
                p.write_text(s.replace(OLD, new), errors="surrogateescape")
print(f"{n} files {'would be ' if dry else ''}rewritten to {new}")
