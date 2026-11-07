"""The validation-cohort table shipped with decoyfdr 0.3.0 (cohort_fvr()), from icb/v20_round15/cohorts_fvr.csv
(round 15 B4). Adds the accession and the count of null models validated; every other value is copied unchanged.

    python3 -I v24_cohort_table.py <cohorts_fvr.csv> <out csv>
"""
import sys
import pandas as pd

F = pd.read_csv(sys.argv[1])
F.insert(2, "accession", F.cohort.str.split("_", n=1).str[1])
nv = F.null_models * F.fvr
assert (nv - nv.round()).abs().max() < 1e-9
F.insert(F.columns.get_loc("fvr"), "null_validated", nv.round().astype(int))
for c in ("patients", "events"):
    assert (F[c] == F[c].round()).all(); F[c] = F[c].astype(int)
assert len(F) == 43 and F.cohort.is_unique
F.to_csv(sys.argv[2], index=False)
print(len(F), "cohorts ->", sys.argv[2])
