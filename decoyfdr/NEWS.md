# decoyfdr 0.2.2

* `dfdr_cohort()` accepts a `SummarizedExperiment` (an assay and the outcome columns of its `colData`); the matrix
  input and every computation are unchanged.
* `dfdr_toy()` simulates a small library; every exported function has a runnable example.
* DESCRIPTION gains `biocViews` and `BugReports`; the vignette labels its chunks and ends with `sessionInfo()`;
  this file is added.

# decoyfdr 0.2.1

* A parallel worker that dies, for instance for lack of memory, stops the run instead of dropping its cohort; memory
  is collected before forking; examples added.

# decoyfdr 0.2.0

* `grow_decoys()` gains `n_sets`, independent decoy sets drawn from substreams of each target's stream (set 1
  reproduces a 0.1.0 run with the same seed); `run_decoy_fdr()` defaults to ten sets.
* `decoy_fdr()` reports the range of the per-set estimates and the pi0-corrected estimate.

# decoyfdr 0.1.0

* First version.
