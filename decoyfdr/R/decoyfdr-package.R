#' decoyfdr: matched decoy modules for false discovery rates of signature libraries
#'
#' Start with [dfdr_cohort()] and [dfdr_library()], then [run_decoy_fdr()]; [calibrate_null()] checks the
#' estimate under the global null. See `vignette("decoyfdr")`.
#'
#' @keywords internal
#' @useDynLib decoyfdr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
