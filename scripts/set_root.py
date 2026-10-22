"""Point every script at a local copy of the data.

    python scripts/set_root.py /path/to/data_root [--dry-run]

The scripts were run on one HPC system and carry its absolute paths. They all hang off one root,
/groups/adv2105_gp/yichen/Yi/com/gastric, under which the analysis keeps icb/ (checkpoint cohorts, null banks, score
stores, package runs and most tables), pancan/ (TCGA), breast/ (the breast cancer cohorts and library) and figs_pub/
(figures). Lay out the inputs listed in INPUTS.csv under a new root with the same sub-directories and run this once;
it rewrites the root in every .R, .py and .sh file of the repository.
"""
import pathlib, sys

OLD = "/groups/adv2105_gp/yichen/Yi/com/gastric"
new = sys.argv[1].rstrip("/")
dry = "--dry-run" in sys.argv
repo = pathlib.Path(__file__).resolve().parents[1]
n = 0
for p in sorted(repo.rglob("*")):
    if p.suffix in (".R", ".py", ".sh") and p.name != "set_root.py":
        s = p.read_text()
        if OLD in s:
            n += 1
            print(f"{p.relative_to(repo)}: {s.count(OLD)} path(s)")
            if not dry:
                p.write_text(s.replace(OLD, new))
print(f"{n} files {'would be ' if dry else ''}rewritten to {new}")
