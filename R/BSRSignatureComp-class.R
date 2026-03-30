#' BulkSignalR ligand-receptor signature object for cluster comparisons
#'
#' S4 class to represent gene signatures associated with ligand-receptor
#' interactions that were inferred from the comparison of two clusters
#' of samples. This class inherits from BSRSignature.
#'
#' @slot cmp.name   The name of the comparison.
#' @slot tg.pval    A list of target genes P-values.
#' @slot tg.logFC   A list of target genes logFC.
#' @slot tg.expr   A list of target genes expression
#'
#' @export
#' @examples
#' new("BSRSignatureComp")
#'
setClass("BSRSignatureComp",
    contains = c("BSRSignature"),
    slots = c(
        cmp.name = "character",
        tg.expr = "list",
        tg.pval = "list",
        tg.logFC = "list"
    ),
    prototype = list(
        cmp.name = "myComparison",
        tg.expr = list(c(1, 2, 3)),
        tg.pval = list(c(0.05, 0.1, 0.008)),
        tg.logFC = list(c(-1, 0, 2))
    )
)

setValidity(
    "BSRSignatureComp",
    function(object) {
        if (!is.character(object@cmp.name)) {
            return("cmp.name is not character")
        }
        if (length(object@cmp.name) == 0) {
            return("cmp.name must have a length > 0")
        }
        if (!is.list(object@tg.expr)) {
            return("tg.expr is not a list")
        }
        if (!is.list(object@tg.pval)) {
            return("tg.pval is not a list")
        }
        if (!is.list(object@tg.logFC)) {
            return("tg.logFC is not a list")
        }

        TRUE
    }
)

setMethod("show", "BSRSignatureComp", function(object) {
    callNextMethod()
    cat("Cluster comparison name:", object@cmp.name, "\n")
})

# Constructor ========================================================

# Obtain gene signatures from a BSRInference object

#' Extract gene signatures of LR pair activity
#'
#' Obtains gene signatures reflecting ligand-receptor as well as
#' receptor downstream activity to
#' score ligand-receptor pairs across samples subsequently with
#' \code{"\link[=BSRInferenceComp-class]{scoreLRGeneSignatures}"}
#'
#' @name BSRSignatureComp
#'
#' @param obj    BSRInferenceComp object.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param with.pw.id    A logical indicating whether the ID of a pathway
#' should be concatenated to its name.
#' @return A BSRSignatureComp object containing a gene signature for each triple
#' ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is
#' comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}.
#' @export
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' 
#' bsrinf.redP <- reduceToPathway(bsrinf.comp)
#' bsrsig.redP <- BSRSignatureComp(bsrinf.redP, qval.thres = 0.001)
#' @importFrom foreach %do% %dopar%
#' @importFrom methods new
BSRSignatureComp <- function(obj,
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
    t.corrs <- tgCorr(obj)[selected]
    t.pvals <- tgPval(obj)[selected]
    t.logFCs <- tgLogFC(obj)[selected]
    t.exprs <- tgExpr(obj)[selected]

    for (i in seq_len(nrow(pairs))) {
        tg <- tg.genes[[i]]
        tg.genes[[i]] <- tg[pairs$rank[i]:length(tg)]
        tc <- t.corrs[[i]]
        t.corrs[[i]] <- tc[pairs$rank[i]:length(tc)]
        tp <- t.pvals[[i]]
        t.pvals[[i]] <- tp[pairs$rank[i]:length(tp)]
        tl <- t.logFCs[[i]]
        t.logFCs[[i]] <- tl[pairs$rank[i]:length(tl)]
        te <- t.exprs[[i]]
        t.exprs[[i]] <- te[pairs$rank[i]:length(te)]

    }

    new("BSRSignatureComp",
        ligands = ligands,
        receptors = receptors, 
        tg.genes = tg.genes, 
        tg.corr = t.corrs,
        pathways = pathways, 
        cmp.name = comparisonName(obj),
        tg.pval = t.pvals, 
        tg.logFC = t.logFCs,        
        tg.expr = t.exprs
    )
} # BSRSignatureComp

# Accessors & setters ========================================================

#' Target gene expression accessor
#'
#' @name tgExpr
#' @aliases tgExpr,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.expr
#' @export
setMethod("tgExpr", "BSRSignatureComp", function(x) x@tg.expr)

#' Target gene P-values accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.pval
#' @export
setMethod("tgPval", "BSRSignatureComp", function(x) x@tg.pval)


#' Target gene logFC accessor
#'
#' @name tgLogFC
#' @aliases tgLogFC,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.logFC
#' @export
setMethod("tgLogFC", "BSRSignatureComp", function(x) x@tg.logFC)

#' Comparison name accessor
#'
#' @name comparisonName
#' @aliases comparisonName,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return cmp.name
#' 
#' @export
setMethod("comparisonName", "BSRSignatureComp", function(x) x@cmp.name)
