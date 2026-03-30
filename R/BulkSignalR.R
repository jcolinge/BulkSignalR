#' BulkSignalR: A package to infer ligand-receptor
#' interactions from bulk transcriptomics or proteomics
#'
#' Inference of ligand-receptor interactions from bulk
#' (transcriptomic or proteomic) data. BulkSignalR bases its inferences
#' on the LRdb database. It relies on a statistical model that
#' is specific to bulk data sets. Different visualization and data
#' summary functions are proposed to help navigating results.
#'
#' @keywords internal 
"_PACKAGE"

utils::globalVariables(c(
    "BulkSignalR_LRdb",
    "%v%", "%>%", "%do%", "%dopar%", "ComplexHeatmap"
))
