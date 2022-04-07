library(methods)

#' @title BulkSignalR ligand-receptor signature Object
#' @description An S4 class to represent gene signatures
#' of inferred ligand-receptor interactions, including
#' their reduced versions.
#'
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot t.genes  A list of target genes, one entry per LR interaction.
#' @slot pathways  An atomic vector of pathway names, one per interaction.
#'
#' @export
#' @examples
#' new("BSRSignatures")
#'
setClass("BSRSignature",
         slots=c(pathways="character",
                 ligands="list",
                 receptors="list",
                 t.genes="list"),
         prototype=list(
             pathways="path 1",
             ligands=list("A"),
             receptors=list("B"),
             t.genes=list(c("a","b","c"))
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

                TRUE
            }
)

setMethod("show", "BSRSignature", function(object) {
    head(
        data.frame(L=sapply(object@ligands, function(x) paste(x,collapse=";")),
                   R=sapply(object@receptors, function(x) paste(x,collapse=";")),
                   pathway=object@pathways,
                   tGenes=sapply(object@t.genes, function(x) paste(x,collapse=";"))
        )
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
#' @title pathways accessor
#' @export
setMethod("pathways", "BSRSignature", function(x) x@pathways)

if (!isGeneric("ligands")) {
    if (is.function("ligands"))
        fun <- ligands
    else
        fun <- function(x) standardGeneric("ligands")
    setGeneric("ligands", fun)
}
#' @title ligands accessor
#' @export
setMethod("ligands", "BSRSignature", function(x) x@ligands)

if (!isGeneric("receptors")) {
    if (is.function("receptors"))
        fun <- receptors
    else
        fun <- function(x) standardGeneric("receptors")
    setGeneric("receptors", fun)
}
#' @title receptors accessor
#' @export
setMethod("receptors", "BSRSignature", function(x) x@receptors)

if (!isGeneric("tGenes")) {
    if (is.function("tGenes"))
        fun <- tGenes
    else
        fun <- function(x) standardGeneric("tGenes")
    setGeneric("tGenes", fun)
}
#' @title Target genes accessor
#' @export
setMethod("tGenes", "BSRSignature", function(x) x@t.genes)

