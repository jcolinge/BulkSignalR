#' BulkSignalR ligand-receptor signature Object
#'
#' S4 class to represent gene signatures
#' of inferred ligand-receptor interactions, including
#' their reduced versions.
#'
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot tg.genes  A list of target genes, one entry per LR interaction.
#' @slot pathways  An atomic vector of pathway names, one per interaction.
#' @slot tg.corr  A list of target genes correlation.
#'
#' @export
#' @examples
#' new("BSRSignature")
setClass("BSRSignature",
    slots = c(
        pathways = "character",
        ligands = "list",
        receptors = "list",
        tg.genes = "list",
        tg.corr = "list"
    ),
    prototype = list(
        pathways = "path 1",
        ligands = list("A"),
        receptors = list("B"),
        tg.genes = list(c("a", "b", "c")),
        tg.corr = list(c(0.1, 0.2, 0.3))
    )
)

setValidity(
    "BSRSignature",
    function(object) {
        if (!is.character(object@pathways)) {
            return("pathways is not character")
        }
        if (!is.list(object@ligands)) {
            return("ligands is not a list")
        }
        if (!is.list(object@receptors)) {
            return("receptors is not a list")
        }
        if (!is.list(object@tg.genes)) {
            return("tg.genes is not a list")
        }
        if (!is.list(object@tg.corr)) {
            return("tg.corr is not a list")
        }
        TRUE
    }
)
setMethod("show", "BSRSignature", function(object) {
    print(data.frame(
        L = vapply(
            object@ligands, function(x) paste(x, collapse = ";"),
            character(1)
        ),
        R = vapply(
            object@receptors, function(x) paste(x, collapse = ";"),
            character(1)
        ),
        pathways = object@pathways,
        tgGenes = vapply(
            object@tg.genes, function(x) paste(x, collapse = ";"),
            character(1)
        )
    )[seq_len(min(5, length(object@ligands))), ])
})

# Constructor ========================================================

#' Extract gene signatures of LR pair activity
#'
#' Obtain gene signatures reflecting ligand-receptor as well as
#' receptor downstream activity to
#' score ligand-receptor pairs across samples subsequently with
#' \code{"\link[=BSRDataModel-class]{scoreLRGeneSignatures}"}
#'
#' @name BSRSignature
#'
#' @param obj    BSRinference object.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param with.pw.id    A logical indicating whether the ID of a pathway
#' should be concatenated to its name.
#' @return A BSRSignature object containing a gene signature for each triple
#' ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is
#' comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}.
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redP <- reduceToPathway(bsrinf)
#' bsrsig.redP <- BSRSignature(bsrinf, qval.thres = 0.001)
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom methods new
BSRSignature <- function(obj,
    pval.thres = NULL, qval.thres = NULL, with.pw.id = FALSE) {
    if (is.null(pval.thres) && is.null(qval.thres)) {
        stop("Either a P- or a Q-value threshold must be provided")
    }

    # reduce and select
    obj <- reduceToBestPathway(obj)
    pairs <- LRinter(obj)
    if (!is.null(pval.thres)) {
        selected <- pairs$pval <= pval.thres
    } else {
        selected <- pairs$qval <= qval.thres
    }

    # obtain the signature object
    pairs <- pairs[selected, ]
    ligands <- ligands(obj)[selected]
    receptors <- receptors(obj)[selected]
    if (with.pw.id) {
        pathways <- paste0(pairs$pw.id, ": ", pairs$pw.name)
    } else {
        pathways <- pairs$pw.name
    }
    tg.genes <- tgGenes(obj)[selected]
    tg.corrs <- tgCorr(obj)[selected]

    for (i in seq_len(nrow(pairs))) {
        tg <- tg.genes[[i]]
        tg.genes[[i]] <- tg[pairs$rank[i]:length(tg)]

        tc <- tg.corrs[[i]]
        tg.corrs[[i]] <- tc[pairs$rank[i]:length(tc)]
    }

    new("BSRSignature",
        ligands = ligands,
        receptors = receptors, 
        tg.genes = tg.genes, 
        tg.corr = tg.corrs,
        pathways =  pathways
    )
} # BSRSignature

# Accessors & setters ========================================================

setGeneric("pathways", signature="x",
    function(x) standardGeneric("pathways")
)
#' pathways accessor
#'
#' @name pathways
#' @aliases pathways,BSRSignature-method
#' @param x BSRSignature
#' @return pathways
#' @examples
#' bsr.sig <- new("BSRSignature")
#' pathways(bsr.sig)
#' @export
setMethod("pathways", "BSRSignature", function(x) x@pathways)

#' ligands accessor
#'
#' @name ligands
#' @aliases ligands,BSRSignature-method
#' @param x BSRSignature
#' @return ligands
#' @examples
#' bsr.sig <- new("BSRSignature")
#' ligands(bsr.sig)
#' @export
setMethod("ligands", "BSRSignature", function(x) x@ligands)

#' receptors accessor
#'
#' @name receptors
#' @aliases receptors,BSRSignature-method
#' @param x BSRSignature
#' @return receptors
#' @examples
#' bsr.sig <- new("BSRSignature")
#' ligands(bsr.sig)
#' @export
setMethod("receptors", "BSRSignature", function(x) x@receptors)

#' Target genes accessor
#'
#' @name tgGenes
#' @aliases tgGenes,BSRSignature-method
#' @param x BSRSignature
#' @return tgGenes
#' @examples
#' bsr.sig <- new("BSRSignature")
#' tgGenes(bsr.sig)
#' @export
setMethod("tgGenes", "BSRSignature", function(x) x@tg.genes)

#' Target gene correlations accessor
#'
#' @name tgCorr
#' @aliases tgCorr,BSRSignature-method
#' @param x BSRSignature
#' @return tgCorr
#' @examples
#' bsr.sig <- new("BSRSignature")
#' tgCorr(bsr.sig)
#' @export
setMethod("tgCorr", "BSRSignature", function(x) x@tg.corr)

