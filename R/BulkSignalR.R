#' BulkSignalR: A package to infer ligand-receptor interactions from bulk transcriptomics or proteomics
#'
#' Inference of ligand-receptor (LR) interactions from bulk
#' (transcriptomic or proteomic) data. BulkSignalR bases its inferences
#' on the LRdb database included in our other package, SingleCellSignalR
#' available from Bioconductor. It relies on a statistical model that
#' is specific to bulk data sets. Different visualization and data
#' summary functions are proposed to help navigating results.
#'
#' @docType package
#' @name BulkSignalR
NULL

globalVariables(c("reactome","gobp","%v%","%>%","%do%","%dopar%","ComplexHeatmap","LRdb"))
