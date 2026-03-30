#' BulkSignalR Cluster Comparison Object
#'
#' An S4 class to represent the comparison of two clusters of samples to
#' infer LR interactions based on the resulting P-values,
#' log-fold-changes (logFC), and expression values.
#'
#' @slot col.clusterA   Column indices for the samples in cluster A.
#' @slot col.clusterB   Column indices for the samples in cluster B.
#' @slot differential.stats  Comparison statistics A versus B 
#' as a data.frame and
#' containing at least 3 columns named 'pval', 'logFC', and 'expr'.
#'
#' @export
#' @examples
#' new("BSRClusterComp")
setClass("BSRClusterComp",
    slots = c(
        col.clusterA = "integer",
        col.clusterB = "integer",
        differential.stats = "data.frame"
    ),
    prototype = list(
        col.clusterA = as.integer(c(1, 2)),
        col.clusterB = as.integer(c(3, 4)),
        differential.stats = data.frame(
            pval = c(0.01, 0.01),
            logFC = c(1, -1), expr = c(1, 2)
        )
    )
)

setValidity(
    "BSRClusterComp",
    function(object) {
        if (!is.integer(object@col.clusterA)) {
            return("col.clusterA indices are not all integers")
        }
        if (length(object@col.clusterA) < 1) {
            return("col.clusterA empty")
        }
        if (!is.integer(object@col.clusterB)) {
            return("col.clusterB indices are not all integers")
        }
        if (length(object@col.clusterB) < 1) {
            return("col.clusterB empty")
        }
        if (length(intersect(object@col.clusterA, 
        object@col.clusterB)) > 0) {
            return("col.cluster1 and col.clusterB are not disjoint")
        }
        if (!is.data.frame(object@differential.stats)) {
            return("specified stats are not a data.frame")
        }

        TRUE
    }
)

setMethod(
    "show", "BSRClusterComp",
    function(object) {
        if (length(object@col.clusterA) > 5) {
            cat("Cluster A columns:", 
            object@col.clusterA[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster A columns:",
                object@col.clusterA[
                    seq_len(length(object@col.clusterA))], "\n"
            )
        }
        if (length(object@col.clusterB) > 5) {
            cat("Cluster B columns:", object@col.clusterB[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster B columns:",
                object@col.clusterB[
                    seq_len(length(object@col.clusterB))], "\n"
            )
        }
        print(utils::head(object@differential.stats))
    }
)

# Constructor ====================================

#' Definition of the comparison between two clusters of samples
#'
#' Define the columns of the expression matrix that belong to each cluster,
#' and store the result of the cluster differences statistical analysis
#' obtained by an external tool such as edgeR or DESeq2 in a dedicated
#' data frame.
#'
#' @name BSRClusterComp
#'
#' @param obj    A BSRDataModelComp object.
#' @param col.clusterA   Cluster A column indices.
#' @param col.clusterB   Cluster B column indices.
#' @param differential.stats  A data.frame containing statistics about
#' the differential
#' analysis cluster A versus B. \code{differentialStats} must contain 
#' at least the
#' columns 'pval' (for P-values), 'logFC' for log-fold-changes A/B, and
#' 'expr' for the expression of the genes in cluster A.
#'
#' @details Create a BSRClusterComp object describing a comparison
#' of two clusters of columns taken from the expression matrix
#' in the BSRDataModelComp object \code{obj}. Such a cluster comparison
#' description is the basis for inferring LRIs from differential
#' expression P-values instead of correlation analysis.
#'
#' The rows of \code{differentialStats} must be in the same order 
#' as those of the count
#' matrix in \code{obj}. Alternatively, \code{differentialStats}
#' rows can be named and a 1-1 correspondence must exist between these names
#' and those of the count matrix.
#'
#' @return A BSRClusterComp object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc, package = "BulkSignalR")
#' normal <- grep("^N", names(sdc))
#' bsrdm <- BSRDataModel(sdc[, -normal])
#'
#' # define the comparison
#' bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
#' colA <- as.integer(1:3)
#' colB <- as.integer(12:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = runif(n, 0, 10)
#' )
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- BSRClusterComp(bsrdm.comp, colA, colB, stats)
#'
#' @importFrom methods new
BSRClusterComp <- function(obj,
    col.clusterA, 
    col.clusterB, 
    differential.stats) {
    if (!is.integer(col.clusterA)) {
        stop("col.clusterA must contain integer indices")
    }
    if (!is.integer(col.clusterB)) {
        stop("col.clusterB must contain integer indices")
    }
    if (length(intersect(col.clusterA, col.clusterB)) > 0) {
        stop("col.clusterA and col.clusterB must be disjoint")
    }
    if (any(col.clusterA < 1 | col.clusterA > ncol(ncounts(obj)))) {
        stop("col.clusterA indices must fall in [1; ncol(ncounts)]")
    }
    if (any(col.clusterB < 1 | col.clusterB > ncol(ncounts(obj)))) {
        stop("col.clusterB indices must fall in [1; ncol(ncounts)]")
    }
    if (!is.data.frame(differential.stats)) {
        stop("differential.stats must be a data.frame")
    }
    if (!all(c("pval", "logFC", "expr") %in% names(differential.stats))) {
        stop("differential.stats data.frame must contain",
            " columns named 'pval', 'logFC', and 'expr'")
    }
    if (nrow(differential.stats) != nrow(ncounts(obj))) {
        stop("differential.stats and ncounts(obj) number of rows differ")
    }
    if (!is.null(rownames(differential.stats)) &&
        (sum(rownames(differential.stats) %in% 
        rownames(ncounts(obj))) != nrow(differential.stats))) {
        stop("differential.stats rownames defined",
            " but do not all match ncounts(obj)")
    }
    if (is.null(rownames(differential.stats))) {
        rownames(differential.stats) <- rownames(ncounts(obj))
    }

    new("BSRClusterComp", col.clusterA = col.clusterA, 
    col.clusterB = col.clusterB, 
    differential.stats = differential.stats)
} # BSRClusterComp


# Accessors & setters ====================================

setGeneric("colClusterA", signature="x",
    function(x) standardGeneric("colClusterA")
)
#' Cluster A columns accessor
#'
#' @name colClusterA
#' @aliases colClusterA,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return col.clusterA
#' @examples
#' bsrcc <- new("BSRClusterComp")
#' colClusterA(bsrcc)
#' @export
setMethod("colClusterA", "BSRClusterComp", function(x) x@col.clusterA)

setGeneric("colClusterA<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colClusterA<-")
)
#' Cluster A columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colClusterA<-", "BSRClusterComp", function(x, value) {
    x@col.clusterA <- value
    methods::validObject(x)
    x
})

setGeneric("colClusterB", signature="x",
    function(x) standardGeneric("colClusterB")
)
#' Cluster B columns accessor
#'
#' @name colClusterB
#' @aliases colClusterB,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return col.clusterB
#' @examples
#' bsrcc <- new("BSRClusterComp")
#' colClusterB(bsrcc)
#' @export
setMethod("colClusterB", "BSRClusterComp", function(x) x@col.clusterB)

setGeneric("colClusterB<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colClusterB<-")
)
#' Cluster B columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colClusterB<-", "BSRClusterComp", function(x, value) {
    x@col.clusterB <- value
    methods::validObject(x)
    x
})


setGeneric("differentialStats", signature="x",
    function(x) standardGeneric("differentialStats")
)
#' Cluster comparison statistics accessor
#'
#' @name differentialStats
#' @aliases differentialStats,BSRClusterComp-method
#' @param x BSRClusterComp object
#' @return diffferential.stats
#' @examples
#' bsrcc <- new("BSRClusterComp")
#' differentialStats(bsrcc)
#' @export
setMethod("differentialStats", "BSRClusterComp", 
function(x) x@differential.stats)

setGeneric("differentialStats<-", signature=c("x", "value"),
    function(x, value) standardGeneric("differentialStats<-")
)
#' Cluster comparison statistics setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("differentialStats<-", "BSRClusterComp", function(x, value) {
    x@differential.stats <- value
    methods::validObject(x)
    x
})
