library(methods)

#' BulkSignalR ligand-receptor signature Object
#'
#' S4 class to represent gene signatures
#' of inferred ligand-receptor interactions, including
#' their reduced versions.
#'
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot t.genes  A list of target genes, one entry per LR interaction.
#' @slot pathways  An atomic vector of pathway names, one per interaction.
#' @slot tg.corr  A list of target genes correlation.
#'
#' @export
#' @examples
#' new("BSRSignature")
#'
setClass("BSRSignature",
         slots=c(pathways="character",
                 ligands="list",
                 receptors="list",
                 t.genes="list",
                 tg.corr="list"),
         prototype=list(
             pathways="path 1",
             ligands=list("A"),
             receptors=list("B"),
             t.genes=list(c("a","b","c")),
             tg.corr=list(c(0.1,0.2,0.3))
         ))

setValidity("BSRSignature",
            function(object) {
                if(!is.character(object@pathways))
                    return("pathways is not character")
                if(!is.list(object@ligands))
                    return("ligands is not a list")
                if(!is.list(object@receptors))
                    return("receptors is not a list")
                if(!is.list(object@t.genes))
                    return("t.genes is not a list")
                if(!is.list(object@tg.corr))
                    return("tg.corr is not a list")
                TRUE
            }
)

setMethod("show", "BSRSignature", function(object) {
    print(data.frame(L=sapply(object@ligands, function(x) paste(x,collapse=";")),
                   R=sapply(object@receptors, function(x) paste(x,collapse=";")),
                   pathways=object@pathways,
                   tGenes=sapply(object@t.genes, function(x) paste(x,collapse=";"))

        )[seq_len(min(5, length(object@ligands))),]
    )
})


# Accessors & setters ========================================================

if (!isGeneric("pathways")) {
    if (is.function("pathways"))
        fun <- pathways
    else
        fun <- function(x) standardGeneric("pathways")
    setGeneric("pathways", fun)
}
#' pathways accessor
#'
#' @name pathways
#' @aliases pathways,BSRSignature-method
#' @param x BSRSignature
#' @export
setMethod("pathways", "BSRSignature", function(x) x@pathways)

if (!isGeneric("ligands")) {
    if (is.function("ligands"))
        fun <- ligands
    else
        fun <- function(x) standardGeneric("ligands")
    setGeneric("ligands", fun)
}
#' ligands accessor
#'
#' @name ligands
#' @aliases ligands,BSRSignature-method
#' @param x BSRSignature
#' @export
setMethod("ligands", "BSRSignature", function(x) x@ligands)

if (!isGeneric("receptors")) {
    if (is.function("receptors"))
        fun <- receptors
    else
        fun <- function(x) standardGeneric("receptors")
    setGeneric("receptors", fun)
}
#' receptors accessor
#'
#' @name receptors
#' @aliases receptors,BSRSignature-method
#' @param x BSRSignature
#' @export
setMethod("receptors", "BSRSignature", function(x) x@receptors)

if (!isGeneric("tGenes")) {
    if (is.function("tGenes"))
        fun <- tGenes
    else
        fun <- function(x) standardGeneric("tGenes")
    setGeneric("tGenes", fun)
}
#' Target genes accessor
#'
#' @name tGenes
#' @aliases tGenes,BSRSignature-method
#' @param x BSRSignature
#' @export
setMethod("tGenes", "BSRSignature", function(x) x@t.genes)

if (!isGeneric("tgCorr")) {
    if (is.function("tgCorr"))
        fun <- tgCorr
    else
        fun <- function(x) standardGeneric("tgCorr")
    setGeneric("tgCorr", fun)
}
#' Target gene correlations accessor
#'
#' @name tgCorr
#' @aliases tgCorr,BSRSignature-method
#' @param x BSRSignature
#' @export
setMethod("tgCorr", "BSRSignature", function(x) x@tg.corr)
