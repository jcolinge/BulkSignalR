library(methods)

#' BulkSignalR Cluster Comparison Object
#'
#' An S4 class to represent the comparison of two clusters of samples to
#' infer LR interactions based on the resulting P-values and
#' log-fold-changes (logFC).
#'
#' @slot colA   Column indices for the samples in cluster A.
#' @slot colB   Column indices for the samples in cluster B.
#' @slot stats  Comparison statistics A versus B as a data.frame and
#' containing at least two columns named 'pval' and 'logFC'.
#' 
#' @export
#' @examples
#' bsrdm <- new("BSRDataModel", ncounts=matrix(1.5, nrow=2, ncol=4,
#'              dimnames=list(c("A","B"), c("E","F","G","H"))),
#'              log.transformed=TRUE, normalization="TC")
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' n <- nrow(ncounts(bsrdm.comp))
#' edger.stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(edger.stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, edger.stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "my_comparison")
#' 
setClass("BSRClusterComp",
         slots=c(colA="integer",
                 colB="integer",
                 stats="data.frame"
         ),
         prototype=list(
           colA=as.integer(1:2),
           colB=as.integer(3:4),
           stats=data.frame(pval=c(0.01,0.01),logFC=c(1,-1))
         ))

setValidity("BSRClusterComp",
            function(object) {
              if(!is.integer(object@colA))
                return("colA indices are not all integers")
              if (length(object@colA) < 1)
                return("colA empty")
              if(!is.integer(object@colB))
                return("colB indices are not all integers")
              if (length(object@colB) < 1)
                return("colB empty")
              if (length(intersect(object@colA, object@colB)) > 0)
                return("colA and colB are not disjoint")
              if(!is.data.frame(object@stats))
                return("specified stats are not a data.frame")

              TRUE
            }
)

setMethod("show", "BSRClusterComp",
          function(object) {
            if (length(object@colA) > 5)
              cat("Cluster A columns:", object@colA[1:5], "...\n")
            else
              cat("Cluster A columns:", object@colA[1:length(object@colA)],"\n")
            if (length(object@colB) > 5)
              cat("Cluster B columns:", object@colB[1:5], "...\n")
            else
              cat("Cluster B columns:", object@colB[1:length(object@colB)],"\n")
            print(utils::head(object@stats))
          }
)


# Accessors & setters ========================================================

if (!isGeneric("colA")) {
  if (is.function("colA"))
    fun <- colA
  else
    fun <- function(x) standardGeneric("colA")
  setGeneric("colA", fun)
}
#' Cluster A columns accessor
#'
#' @name colA
#' @aliases colA,BSRClusterComp-method
#' @param x object BSRClusterComp 
#' @export
setMethod("colA", "BSRClusterComp", function(x) x@colA)

if (!isGeneric("colA<-")) {
  if (is.function("colA<-"))
    fun <- `colA<-`
  else
    fun <- function(x, value) standardGeneric("colA<-")
  setGeneric("colA<-", fun)
}
#' Cluster A columns setter (internal use only)
#'
#' @param x object BSRClusterComp 
#' @param value value to be set for BSRClusterComp
#' @keywords internal 
setMethod("colA<-", "BSRClusterComp", function(x,value){
  x@colA <- value
  methods::validObject(x)
  x
})


if (!isGeneric("colB")) {
  if (is.function("colB"))
    fun <- colB
  else
    fun <- function(x) standardGeneric("colB")
  setGeneric("colB", fun)
}
#' Cluster B columns accessor
#'
#' @name colB
#' @aliases colB,BSRClusterComp-method
#' @param x object BSRClusterComp 
#' @export
setMethod("colB", "BSRClusterComp", function(x) x@colB)

if (!isGeneric("colB<-")) {
  if (is.function("colB<-"))
    fun <- `colB<-`
  else
    fun <- function(x, value) standardGeneric("colB<-")
  setGeneric("colB<-", fun)
}
#' Cluster B columns setter (internal use only)
#'
#' @param x object BSRClusterComp 
#' @param value value to be set for BSRClusterComp
#' @keywords internal 
setMethod("colB<-", "BSRClusterComp", function(x,value){
  x@colB <- value
  methods::validObject(x)
  x
})


if (!isGeneric("stats")) {
  if (is.function("stats"))
    fun <- stats
  else
    fun <- function(x) standardGeneric("stats")
  setGeneric("stats", fun)
}
#' Cluster comparison statistics accessor
#'
#' @name stats
#' @aliases stats,BSRClusterComp-method
#' @param x BSRClusterComp object
#' @export
setMethod("stats", "BSRClusterComp", function(x) x@stats)

if (!isGeneric("stats<-")) {
  if (is.function("stats<-"))
    fun <- `stats<-`
  else
    fun <- function(x, value) standardGeneric("stats<-")
  setGeneric("stats<-", fun)
}
#' Cluster comparison statistics setter (internal use only)
#'
#' @param x object BSRClusterComp 
#' @param value value to be set for BSRClusterComp
#' @keywords internal 
setMethod("stats<-", "BSRClusterComp", function(x,value){
  x@stats <- value
  methods::validObject(x)
  x
})
