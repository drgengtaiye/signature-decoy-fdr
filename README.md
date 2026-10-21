# signature-decoy-fdr

Analysis code for **Matched decoy modules give immune gene-expression signatures a false discovery rate**
(manuscript under review).

52 published immune gene-expression signatures are each compared with random gene sets matched on
size and scoring engine, in 12 checkpoint-treated cohorts and 20 TCGA cancer types. 2,912 decoy
modules, grown from genes in no published signature and matched to each signature on size, engine,
expression level and internal co-expression, then attach a false discovery rate to the library's
claims.

This repository holds the same files as Additional file 3 of the article. `README.txt` describes
the layout and the order in which the scripts run. The scripts carry the absolute paths of the HPC
system they were run on; set the data root at the top of each script (and `D` in
`scripts/v4kit.py`) before running elsewhere. All input data are public; the article's data
availability statement lists the sources.

| Folder | Contents |
|---|---|
| `scripts/` | benchmarks, dependence, randomized-trial analyses and figures |
| `decoys/` | the target-decoy calibration |
| `trials/` | treatment-interaction models and the no-interaction simulation |
| `validation/` | deposited analysis plans with their SHA-256 hashes, the GEO search, I-SPY2, the seven registered cohorts of tests A and B |
| `manuscript/` | build the manuscript and Additional files |

License: MIT (see `LICENSE`).
