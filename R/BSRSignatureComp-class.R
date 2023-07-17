library(methods)

#' BulkSignalR ligand-receptor signature object for cluster comparisons
#'
#' S4 class to represent gene signatures associated with ligand-receptor
#' interactions that were inferred from the comparison of two clusters
#' of samples. This class inherits from BSRSignature.
#'
#' @slot cmp.name   The name of the comparison.
#' @slot tg.pval    A list of target genes P-values.
#'
#' @export
#' @examples
#' new("BSRSignatureComp")
#'
setClass("BSRSignatureComp",
         contains = c("BSRSignature"),
         slots=c(cmp.name="character",
                 tg.pval="list"),
         prototype=list(
           cmp.name="happy",
           pathways="path 1",
           ligands=list("A"),
           receptors=list("B"),
           t.genes=list(c("a","b","c")),
           tg.corr=list(c(0.1,0.2,0.3)),
           tg.pval=list(c(0.05,0.1,0.008))
         ))

setValidity("BSRSignatureComp",
            function(object) {
              if(!is.character(object@cmp.name))
                return("cmp.name is not character")
              if (length(object@cmp.name) == 0)
                return("cmp.name must have a length > 0")
              if(!is.list(object@tg.pval))
                return("tg.pval is not a list")
              
              TRUE
            }
)

setMethod("show", "BSRSignatureComp", function(object) {
  callNextMethod()
  cat("Cluster comparison name:", object@cmp.name, "\n")
})


# Accessors & setters ========================================================

if (!isGeneric("tgPval")) {
  if (is.function("tgPval"))
    fun <- tgPval
  else
    fun <- function(x) standardGeneric("tgPval")
  setGeneric("tgPval", fun)
}
#' Target genes accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRSignatureComp-method
#' @param x BSRSignature
#' @export
setMethod("tgPval", "BSRSignatureComp", function(x) x@tg.pval)

if (!isGeneric("cmpName")) {
  if (is.function("cmpName"))
    fun <- tgcmpName
  else
    fun <- function(x) standardGeneric("cmpName")
  setGeneric("cmpName", fun)
}
#' Target genes accessor
#'
#' @name cmpName
#' @aliases cmpName,BSRSignatureComp-method
#' @param x BSRSignature
#' @export
setMethod("cmpName", "BSRSignatureComp", function(x) x@cmp.name)
