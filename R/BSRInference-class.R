#' BulkSignalR Inference Object
#'
#' An S4 class to represent inferred ligand-receptor interactions.
#'
#' @slot LRinter  A data frame describing the (ligand,receptor,pathway) triples
#' with P- and Q-values.
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot tg.genes  A list of target genes, one entry per LR interaction.
#' @slot tg.corr  A list of target gene correlations to the receptor, one
#' entry per interaction
#' @slot inf.param  The parameters used for the inference.
#'
#' @details This class contains inferred LR interactions along with
#' their statistical significance. Data representation supports subsequent
#' reductions to pathways, etc. See reduction functions
#' \code{"\link[=BSRInference-class]{reduceToBestPathway}"},
#' \code{"\link[=BSRInference-class]{reduceToLigand}"},
#' \code{"\link[=BSRInference-class]{reduceToReceptor}"} and
#' \code{"\link[=BSRInference-class]{reduceToPathway}"}.
#' @export
#' @examples
#' new("BSRInference")
setClass("BSRInference",
    slots = c(
        LRinter = "data.frame",
        ligands = "list",
        receptors = "list",
        tg.genes = "list",
        tg.corr = "list",
        inf.param = "list"
    ),
    prototype = list(
        LRinter = data.frame(
            L = "A", R = "B", pw.id = "123", pw.name = "one pw",
            pval = 1.0, qval = 1.0, LR.corr = 0.6, rank = 2,
            len = 50, rank.corr = 0.6,
            stringsAsFactors = FALSE
        ),
        ligands = list("A"),
        receptors = list("B"),
        tg.genes = list(c("a", "b", "c")),
        tg.corr = list(c(-0.5, 0.1, 0.8)),
        inf.param = list()
    )
)


setValidity(
    "BSRInference",
    function(object) {
        if (!is.data.frame(object@LRinter)) {
            return("LRinter is not a data frame")
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
        if (!is.list(object@inf.param)) {
            return("inf.param is not a list")
        }

        TRUE
    }
)


setMethod(
    "show", "BSRInference",
    function(object) {
        cat("Reference database: ", object@inf.param$reference, "\n", sep = "")
        print(utils::head(object@LRinter[
            order(object@LRinter$qval),
            c("L", "R", "pval", "qval", "pw.id", "pw.name"),
        ]))[5, ]
        cat("Inference parameters:\n")
        utils::str(object@inf.param)
    }
)


# Constructor ========================================================

#' Inference of ligand-receptor interactions
#'
#' Computes putative LR interactions along with their statistical confidence.
#' In this initial inference, all the relevant pathways are reported,
#' see reduction functions to reduce this list.
#'
#' @name BSRInference
#'
#' @param obj         A BSRDataModel output by \code{\link{BSRDataModel}} with
#' statistical model parameters trained by the
#' \code{"\link[=BSRDataModel-class]{learnParameters}"} method.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#' considered target genes.
#' @param min.cor         The minimum Spearman correlation required between
#' the ligand and the receptor.
#' @param fdr.proc      The procedure for adjusting P-values according to
#' \code{\link[multtest]{mt.rawp2adjp}}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#' @param use.full.network  A logical to avoid limiting the reference network
#'   to the detected genes and use the whole reference network.
#'
#' @details Perform the initial ligand-receptor inference. Initial means that
#' no reduction is applied. All the (ligand, receptor, downstream pathway)
#' triples are reported, i.e., a given LR pair may appear multiple times
#' with different pathways downstream the receptor. Specific reduction
#' functions are available from the package to operate subsequent
#' simplifications based on the BSRInference object created by the initial
#' inference.
#'
#' Parameters defining minimum/maximum pathway sizes, etc. are set to NULL
#' by default, meaning that their values will be taken from what was set
#' during the training of the statistical model with
#' \code{"\link[=BSRDataModel-class]{learnParameters}"}
#'
#' To use different
#' values at the time of inference sounds like a bad idea, although this
#' could be used to explore without retraining the underlying model.
#' Retraining of the model with adjusted parameters is advised following
#' such an exploration.
#'
#' @return A BSRInference object with initial inferences set.
#'
#' @export
#'
#' @examples
#' data(bsrdm, package = "BulkSignalR")
#' data(immune.signatures, package = "BulkSignalR")
#' 
#' # We use a subset of the reference to speed up
#' # inference in the context of the example.
#' 
#' reactSubset <- getResource(resourceName = "Reactome",
#' cache = FALSE)
#' 
#' subset <- c("REACTOME_BASIGIN_INTERACTIONS",
#' "REACTOME_SYNDECAN_INTERACTIONS",
#' "REACTOME_ECM_PROTEOGLYCANS",
#' "REACTOME_CELL_JUNCTION_ORGANIZATION")
#' 
#' reactSubset <- reactSubset[
#' reactSubset$`Reactome name` %in% subset,]
#' 
#' resetPathways(dataframe = reactSubset,
#' resourceName = "Reactome")
#' 
#' bsrinf <- BSRInference(bsrdm,
#'     min.cor = 0.2,restrict.genes=immune.signatures$gene,
#'     reference="REACTOME")
#' @importFrom methods new
BSRInference <- function(obj, rank.p = 0.55,
    min.cor = 0.25, restrict.genes = NULL,
    reference = c("REACTOME-GOBP", "REACTOME", "GOBP"),
    max.pw.size = NULL, min.pw.size = NULL, min.positive = NULL,
    use.full.network = FALSE, restrict.pw = NULL,
    with.complex = NULL, fdr.proc = c(
        "BH", "Bonferroni", "Holm", "Hochberg",
        "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {

    if (is.null(max.pw.size)) {
        max.pw.size <- parameters(obj)$max.pw.size
    }
    if (is.null(min.pw.size)) {
        min.pw.size <- parameters(obj)$min.pw.size
    }
    if (is.null(min.positive)) {
        min.positive <- parameters(obj)$min.positive
    }
    if (is.null(with.complex)) {
        with.complex <- parameters(obj)$with.complex
    }
    if (is.null(use.full.network)) {
        use.full.network <- parameters(obj)$use.full.network
    }
    reference <- match.arg(reference)
    fdr.proc <- match.arg(fdr.proc)
    if (rank.p < 0 || rank.p > 1) {
        stop("rank.p must lie in [0;1]")
    }

    inf.param <- list()
    inf.param$min.corr <- min.cor
    inf.param$restrict.genes <- restrict.genes
    lr <- .getCorrelatedLR(obj, min.cor = min.cor, 
        restrict.genes = restrict.genes)

    inf.param$reference <- reference
    inf.param$min.pw.size <- min.pw.size
    inf.param$max.pw.size <- max.pw.size
    inf.param$with.complex <- with.complex
    inf.param$min.positive <- min.positive
    inf.param$use.full.network <- use.full.network
    inf.param$restrict.pw <- restrict.pw
    pairs <- .checkReceptorSignaling(obj, lr,
        reference = reference,
        min.pw.size = min.pw.size, max.pw.size = max.pw.size,
        min.positive = min.positive, with.complex = with.complex,
        use.full.network = use.full.network, restrict.pw = restrict.pw
    )

    inf.param$fdr.proc <- fdr.proc
    inf.param$rank.p <- rank.p
    inter <- .pValuesLR(pairs, 
        parameters(obj), 
        rank.p = rank.p, 
        fdr.proc = fdr.proc)

    ligands <- strsplit(inter$L, ";")
    receptors <- strsplit(inter$R, ";")
    tg <- strsplit(inter$target.genes, ";")
    tgcorr <- lapply(
        strsplit(inter$target.corr, ";"),
        function(x) as.numeric(x)
    )
    inf.param$ligand.reduced <- FALSE
    inf.param$receptor.reduced <- FALSE
    inf.param$pathway.reduced <- FALSE

    new("BSRInference",
        LRinter = inter[, c(
            "L", "R", "pw.id", "pw.name", "pval", "qval",
            "LR.corr", "rank", "len", "rank.corr"
        )], ligands = ligands,
        receptors = receptors, tg.genes = tg, tg.corr = tgcorr,
        inf.param = inf.param
    )
    
} # BSRInference


# Accessors & setters ========================================================

setGeneric("LRinter", signature="x",
    function(x) standardGeneric("LRinter")
)
#' LRinter accessor
#'
#' @name LRinter
#' @aliases LRinter,BSRInference-method
#' @param x BSRInference object
#' @return LRinter
#' @examples
#' bsrinf <- new ("BSRInference")
#' LRinter(bsrinf)
#' @export
setMethod("LRinter", "BSRInference", function(x) x@LRinter)


setGeneric("LRinter<-", signature=c("x", "value"),
    function(x, value) standardGeneric("LRinter<-")
)
#' LRinter setter (internal use only)
#'
#' @param x BSRInference object
#' @param value value to be set to BSRInference
#' @return returns \code{NULL}
#' @keywords internal
setMethod("LRinter<-", "BSRInference", function(x, value) {
    x@LRinter <- value
    methods::validObject(x)
    x
})


setGeneric("ligands", signature="x",
    function(x) standardGeneric("ligands")
)
#' ligands accessor
#'
#' @name ligands
#' @aliases ligands,BSRInference-method
#' @param x BSRInference object
#' @return ligands
#' @export
setMethod("ligands", "BSRInference", function(x) x@ligands)


setGeneric("ligands<-", signature=c("x", "value"),
    function(x, value) standardGeneric("ligands<-")
)
#' ligands setter (internal use only)
#' @param x BRSInference object
#' @param value Value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("ligands<-", "BSRInference", function(x, value) {
    x@ligands <- value
    methods::validObject(x)
    x
})


setGeneric("receptors", signature="x",
    function(x) standardGeneric("receptors")
)
#' receptors accessor
#'
#' @name receptors
#' @aliases receptors,BSRInference-method
#' @param x BRSInference object
#' @return receptors
#' @export
setMethod("receptors", "BSRInference", function(x) x@receptors)


setGeneric("receptors<-", signature=c("x", "value"),
    function(x, value) standardGeneric("receptors<-")
)
#' receptors setter (internal use only)
#'
#' @param x BRSInference object
#' @param value value to be set for BRSInference
#' @return returns \code{NULL}
#' @keywords internal
setMethod("receptors<-", "BSRInference", function(x, value) {
    x@receptors <- value
    methods::validObject(x)
    x
})


setGeneric("tgGenes", signature="x",
    function(x) standardGeneric("tgGenes")
)
#' Target genes accessor
#'
#' @name tgGenes
#' @aliases tgGenes,BSRInference-method
#' @param x BSRInferance object
#' @return tgGenes
#' @export
setMethod("tgGenes", "BSRInference", function(x) x@tg.genes)


setGeneric("tgGenes<-", signature=c("x", "value"),
    function(x, value) standardGeneric("tgGenes<-")
)
#' Target genes setter (internal use only)
#' 
#' @param x BSRInference object
#' @param value value to be set BSRInference
#' @return returns \code{NULL}
#' @keywords internal
setMethod("tgGenes<-", "BSRInference", function(x, value) {
    x@tg.genes <- value
    methods::validObject(x)
    x
})


setGeneric("tgCorr", signature="x",
    function(x) standardGeneric("tgCorr")
)
#' Target gene correlations accessor
#'
#' @name tgCorr
#' @aliases tgCorr,BSRInference-method
#' @param x BSRInference object
#' @export
setMethod("tgCorr", "BSRInference", function(x) x@tg.corr)

setGeneric("tgCorr<-", signature=c("x", "value"),
    function(x, value) standardGeneric("tgCorr<-")
)
#' Target gene correlations setter (internal use only)
#' @param x BSRInference object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' 
#' @keywords internal
setMethod("tgCorr<-", "BSRInference", function(x, value) {
    x@tg.corr <- value
    methods::validObject(x)
    x
})


setGeneric("inferenceParameters", signature="x",
    function(x) standardGeneric("inferenceParameters")
)
#' Inference parameters accessor
#'
#' @name inferenceParameters
#' @aliases inferenceParameters,BSRInference-method
#' @param x BRSInference object.
#' @return inf.param
#' @examples
#' bsrinf <- new ("BSRInference")
#' inferenceParameters(bsrinf)
#' @export
setMethod("inferenceParameters", "BSRInference", function(x) x@inf.param)


setGeneric("inferenceParameters<-", signature=c("x", "value"),
    function(x, value) standardGeneric("inferenceParameters<-")
)
#' Inference parameters setter (internal use only)
#' @param x BRSInference object.
#' @param value value to be set.
#' @return returns \code{NULL}
#' 
#' @keywords internal
setMethod("inferenceParameters<-", "BSRInference", function(x, value) {
    x@inf.param <- value
    methods::validObject(x)
    x
})


# simplified table view ========================================================

setGeneric("LRinterShort", signature="x",
    function(x) standardGeneric("LRinterShort")
)
#' Simplified LRinter accessor reporting the essential columns
#'
#' @name LRinterShort
#' @aliases LRinterShort,BSRInference-method
#' @param x BSRInference object
#' @return LRinterShort
#' @export
setMethod(
    "LRinterShort", "BSRInference",
    function(x) {
        x@LRinter[, c(
            "L", "R", "pw.id", "pw.name",
            "qval", "LR.corr", "len"
        )]
    }
)


# Rescoring ====================================================================

setGeneric("rescoreInference", signature="obj",
    function(obj,...) standardGeneric("rescoreInference")
)
#' Inference re-scoring
#'
#' A method to re-score an existing BSRInference object
#' (P- and Q-value estimations).
#'
#' @name rescoreInference
#' @aliases rescoreInference,BSRInference-method
#'
#' @param obj BSRInferecence object.
#' @param param The matching BSRDataModel parameters
#'              (can be obtained with the accessor \code{param})
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @details A BSRInference object should be created by calling
#' \code{"\link[=BSRDataModel-class]{BSRInference}"}
#'
#' Parameters controlling the estimation
#' of the statistical significance of the ligand/receptor/pathway triples
#' (\code{param}) are provided at the time of calling the latter method.
#'
#' Nonetheless, it
#' might be useful to change the initially-provided parameters, in
#' which case this method should not be called.
#'
#' @return A BSRInference object.
#'
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' data(bsrdm, package = "BulkSignalR")
#' 
#' bsrinf.new <- rescoreInference(bsrinf,
#'                       param = parameters(bsrdm)
#'               )
setMethod("rescoreInference", "BSRInference", function(obj, param, 
    rank.p = 0.55, fdr.proc = c("BH", "Bonferroni", "Holm",
        "Hochberg", "SidakSS", "SidakSD", "BY", "ABH", "TSBH"
        )) {

    if (rank.p < 0 || rank.p > 1) {
        stop("rank.p must lie in [0;1]")
    }
    fdr.proc <- match.arg(fdr.proc)

    # extract the necessary data from the BSRInference object
    pairs <- LRinter(obj)
    tg.genes <- tgGenes(obj)
    tg.corr <- tgCorr(obj)

    # prepare the chosen model CDF
    LR.par <- param$LR.0$model
    RT.par <- param$RT.0$model
    if (LR.par$distrib != RT.par$distrib) {
        stop("Distinct statistical models for LR and RT nulls are not allowed")
    }
    if (LR.par$distrib == "censored_normal") {
        cdf <- .cdfGaussian
    } else if (LR.par$distrib == "censored_mixed_normal") {
        cdf <- .cdfMixedGaussian
    } else if (LR.par$distrib == "empirical") {
        cdf <- .cdfEmpirical
    } else if (LR.par$distrib == "kernel_empirical") {
        cdf <- .cdfKernelEmpirical
    } else if (LR.par$distrib == "censored_stable") {
        cdf <- .cdfAlphaStable
    } else {
        stop("Unknown statistical model: ", LR.par$LR.0$model$distrib)
    }

    # recompute P-values
    for (i in seq_len(nrow(pairs))) {
        tg <- tg.genes[[i]]
        spears <- tg.corr[[i]]

        # estimate the LR correlation P-value
        if (pairs$LR.corr[i] >= 0) {
            # normal case
            p.lr <- 1 - cdf(pairs$LR.corr[i], LR.par)
        } else {
            # to enable searching for inhibitory L-R interactions
            p.lr <- cdf(pairs$LR.corr[i], LR.par)
        }

        # estimate the target gene correlation P-value based on rank statistics
        # for the individual correlation Gaussian model
        len <- pairs$len[i]
        r <- min(max(1, trunc(rank.p * len)), len)
        pairs$rank[i] <- r
        rank.corr <- spears[r]
        p.rt <- stats::pbinom(r - 1, len, cdf(rank.corr, RT.par))
        pairs$pval[i] <- p.lr * p.rt
        pairs$rank.corr[i] <- rank.corr
    }

    # recompute the Q-values
    rawp <- pairs$pval
    if (length(rawp) > 1){
        adj <- multtest::mt.rawp2adjp(rawp, fdr.proc)
        pairs$qval <- adj$adjp[order(adj$index), fdr.proc]
    }
    else {
        pairs$qval <- rawp
    }

    # update the BSRInference object
    inf.param <- inferenceParameters(obj)
    inf.param$fdr.proc <- fdr.proc
    inf.param$rank.p <- rank.p
    inferenceParameters(obj) <- inf.param
    LRinter(obj) <- pairs

    obj
    
}) # rescoreInference


# Reduction and pathway stat methods ===========================================

setGeneric("getPathwayStats", signature="obj",
    function(obj,...) standardGeneric("getPathwayStats")
)
#' Basic statistics about hit pathways
#'
#' @name getPathwayStats
#' @aliases getPathwayStats,BSRInference-method
#'
#' @param obj    BSRinf object.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @return A table with the pathways selected after the chosen threshold was
#'   applied to rows in \code{LRinter(obj)}.
#' Each pathway is reported along with various statistics:
#' the number of selected receptors in this pathway, the total number of
#' receptors described this pathway,
#' the number of selected ligand-receptor pairs hitting this pathway,
#' and the total number of
#' ligand-receptor pairs described that could hit this pathway.
#'
#' Obviously, one could imagine computing enrichment in receptors or
#' ligand-receptor pairs
#' based on such statistics, but the actual meaning of such an analysis
#' would be ambiguous since
#' the pathways were already selected as significantly regulated by the
#' receptor. We thus did not implement
#' this (hypergeometric test) computation.
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#'
#' pw.stat <- getPathwayStats(bsrinf)
#' 
#'
#' @importFrom foreach %do% %dopar%
setMethod("getPathwayStats", "BSRInference", function(obj,
    pval.thres = NULL, qval.thres = NULL) {
    if (inferenceParameters(obj)$ligand.reduced || 
        inferenceParameters(obj)$receptor.reduced) {
        stop(
            "Cannot be applied to interactions involving",
            " reduced receptors or ligands"
        )
    }

    # local binding
    id <- NULL

    pairs <- LRinter(obj)
    if (!is.null(pval.thres)) {
        pairs <- pairs[pairs$pval <= pval.thres, ]
    } else {
        pairs <- pairs[pairs$qval <= qval.thres, ]
    }

    # pairs reduced to the receptor & pathway names
    pairs.R <- unique(pairs[, c("R", "pw.id", "pw.name")])
    upw <- unique(pairs[, c("pw.id", "pw.name")])
    pw.names <- stats::setNames(upw$pw.name, upw$pw.id)

    # number of "hits" on a pathway with or without the ligand combinatorics
    pw.ids <- table(pairs$pw.id)
    pw.ids.R <- table(pairs.R$pw.id)

    # number of ligands for each receptor
    R.n.comb <- table(.SignalR$BulkSignalR_LRdb$receptor)

    foreach::foreach(id = names(pw.ids), .combine = rbind) %do% {
        # number of receptors that are in the current pathway,
        # depending on whether it is a GOBP or Reactome pathway
        if (regexpr("^R-", id) != -1) {
            Rs <- intersect(
                .SignalR$BulkSignalR_LRdb$receptor,
                .SignalR$BulkSignalR_Reactome[
                .SignalR$BulkSignalR_Reactome$`Reactome ID` == id, "Gene name"]
            )
        } else {
            Rs <- intersect(
                .SignalR$BulkSignalR_LRdb$receptor,
                .SignalR$BulkSignalR_Gobp[
                .SignalR$BulkSignalR_Gobp$`GO ID` == id, "Gene name"]
            )
        }

        # non-combinatorial version (ignore ligands)
        tot.R <- length(Rs)

        # ligand combinatorics included
        tot.LR <- sum(R.n.comb[Rs])

        data.frame(
            pw.id = id, pw.name = pw.names[id], n.R = pw.ids.R[id],
            tot.R = tot.R, n.LR = pw.ids[id], tot.LR = tot.LR,
            stringsAsFactors = FALSE
        )
    }
    
}) # getPathwayStats


setGeneric("reduceToBestPathway", signature="obj",
    function(obj,...) standardGeneric("reduceToBestPathway")
)
#' Keep one pathway per ligand-receptor pair
#'
#' @name reduceToBestPathway
#' @aliases reduceToBestPathway,BSRInference-method
#'
#' @param obj BSRInference object
#'
#' @return A BSRInference object reduced to only report one pathway per
#' ligand-receptor pair. The pathway with the
#' smallest P-value is selected.
#'
#' @details Ligand-receptor pairs
#' are evaluated in relation with pathways that allow checking receptor
#' downstream correlations. It is thus possible
#' that several pathways are reported for a same LR pair.
#'
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' bsrinf.redBP <- reduceToBestPathway(bsrinf)
#'
#' @importFrom rlang .data
setMethod("reduceToBestPathway", "BSRInference", function(obj) {
    # Here we access the object slots directly as this procedure
    # is dependent of actual data representation

    # get best p-value pathway per LR pair
    ligands <- list()
    receptors <- list()
    tg.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    LR <- unique(pairs[, c("L", "R")])
    for (i in seq_len(nrow(LR))) {
        L <- LR$L[i]
        R <- LR$R[i]

        pwr <- pairs[pairs$L == L & pairs$R == R, ]
        k <- which.min(pwr$pval)
        j <- which(pairs$L == L & pairs$R == R & pairs$pw.id == pwr$pw.id[k])
        ligands <- c(ligands, obj@ligands[j])
        receptors <- c(receptors, obj@receptors[j])
        tg.genes <- c(tg.genes, obj@tg.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        LRinter <- rbind(LRinter, pairs[j, ])
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE

    obj
    
}) # reduceToBestPathway


setGeneric("reduceToReceptor", signature="obj",
    function(obj,...) standardGeneric("reduceToReceptor")
)
#' Aggregate the ligands of a same receptor
#'
#' Simplifies a ligand-receptor table to focus on the receptors.
#'
#' @name reduceToReceptor
#' @aliases reduceToReceptor,BSRInference-method
#'
#' @return BSRInference object reduced to one row per receptor.
#' All the ligands are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the line with the
#' pathway featuring the smallest P-value.
#' @param obj BRSInference object
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redR <- reduceToReceptor(bsrinf)
#'
#' @importFrom rlang .data
setMethod("reduceToReceptor", "BSRInference", function(obj) {
    # Here we access the object slots directly as this procedure
    # is dependent of actual data representation

    if (inferenceParameters(obj)$ligand.reduced) {
        stop("Already reduced to receptor")
    } # because ligands were reduced

    # pool the ligands
    ligands <- list()
    receptors <- list()
    tg.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (R in unique(pairs$R)) {
        lig <- pairs[pairs$R == R, ]
        k <- which.min(lig$pval)
        j <- which(pairs$R == lig$R[k] & pairs$pw.id == lig$pw.id[k])[1]
        ligands <- c(ligands, list(unique(lig$L)))
        receptors <- c(receptors, list(R))
        tg.genes <- c(tg.genes, obj@tg.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j, ]
        to.add[1, "L"] <- paste0("{", paste(unique(lig$L), collapse = ";"), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$ligand.reduced <- TRUE

    obj
    
}) # reduceToReceptor


setGeneric("reduceToLigand", signature="obj",
    function(obj,...) standardGeneric("reduceToLigand")
)
#' Aggregate the receptors of a same ligand
#'
#' Simplifies a ligand-receptor table to focus on the ligands.
#'
#' @name reduceToLigand
#' @aliases reduceToLigand,BSRInference-method
#'
#' @return A BSRInference object but reduced to one row per ligand.
#' All the receptors are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the pathway with
#' the smallest P-value.
#' @param obj BSRInference object
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redL <- reduceToLigand(bsrinf)
#'
#' @importFrom rlang .data
setMethod("reduceToLigand", "BSRInference", function(obj) {
    # Here we access the object slots directly as this procedure
    # is dependent of actual representation

    if (inferenceParameters(obj)$receptor.reduced) {
        stop("Already reduced to ligand")
    } # because receptors were reduced

    # pool the receptors
    ligands <- list()
    receptors <- list()
    tg.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (L in unique(pairs$L)) {
        rec <- pairs[pairs$L == L, ]
        k <- which.min(rec$pval)
        j <- which(pairs$L == L & 
            pairs$R == rec$R[k] & 
            pairs$pw.id == rec$pw.id[k])
        ligands <- c(ligands, list(L))
        receptors <- c(receptors, list(unique(rec$R)))
        tg.genes <- c(tg.genes, obj@tg.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j, ]
        to.add[1, "R"] <- paste0("{", paste(unique(rec$R), collapse = ";"), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj
    
}) # reduceToLigand


setGeneric("reduceToPathway", signature="obj",
    function(obj,...) standardGeneric("reduceToPathway")
)
#' Aggregate ligands and receptors at the pathway level
#'
#' Simplifies a ligand-receptor inference object to focus on
#' the pathways.
#'
#' @name reduceToPathway
#' @aliases reduceToPathway,BSRInference-method
#'
#' @return A BSRInference object reduced to only report one row per pathway.
#' The information of which ligand interacted with which receptor is lost as
#' all the ligands and all the receptors forming pairs related to a certain
#' pathway are combined.
#' For a given pathway, the reported P-values and target genes are those of
#' the best ligand-receptor pair that
#' was in this pathway.
#' Receptors and ligands are combined in two semi-colon-separated
#' lists surrounded by curly brackets in the tabular slot \code{LRinter},
#' while the list representation slots (\code{ligands} and
#' \code{receptors}) are update accordingly.
#'
#' @param obj BSRInference object
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redP <- reduceToPathway(bsrinf)
#' @importFrom rlang .data
setMethod("reduceToPathway", "BSRInference", function(obj) {
    # Here we access the object slots directly as this procedure
    # is dependent of actual representation

    if (inferenceParameters(obj)$receptor.reduced) {
        stop("Already reduced to ligand")
    } # because receptors were reduced
    if (inferenceParameters(obj)$ligand.reduced) {
        stop("Already reduced to receptor")
    } # because ligands were reduced

    # reduce to unique pathways
    ligands <- list()
    receptors <- list()
    tg.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (p in unique(pairs$pw.id)) {
        j <- which(pairs$pw.id == p)[1]
        ligands <- c(ligands, list(unique(pairs$L[pairs$pw.id == p])))
        receptors <- c(receptors, list(unique(pairs$R[pairs$pw.id == p])))
        tg.genes <- c(tg.genes, obj@tg.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j, ]
        to.add[1, "L"] <- paste0("{", paste(unique(pairs$L[pairs$pw.id == p]),
            collapse = ";"
        ), "}")
        to.add[1, "R"] <- paste0("{", paste(unique(pairs$R[pairs$pw.id == p]),
            collapse = ";"
        ), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$ligand.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj
    
}) # reduceToPathway


# Reset gene names to initial organism =========================================

setGeneric("resetToInitialOrganism", signature="obj",
    function(obj,...) standardGeneric("resetToInitialOrganism")
)
#'  Reset gene names to initial organism provided in the first instance
#'
#' @name resetToInitialOrganism
#' @aliases resetToInitialOrganism,BSRInference-method
#'
#' @param obj  BSRInference object
#' @param conversion.dict   A dictionnary
#'
#' @return An BSRInference object updated for gene names.
#' The gene names are replaced by the ones from
#' the organism provided in the first instance.
#'
#' @export
#' @examples
#' data(bodyMap.mouse, package = "BulkSignalR")
#' data(bsrinf.mouse, package = "BulkSignalR")
#' data(ortholog.dict, package = "BulkSignalR")
#' 
#' #idx <- sample(nrow(bodyMap.mouse), 7500)
#' 
#' #bodyMap.mouse <- bodyMap.mouse[idx,1:3]
#' 
#' #ortholog.dict <- findOrthoGenes(
#' #    from_organism = "mmusculus",
#' #    from_values = rownames(bodyMap.mouse)
#' #)
#'
#' #matrix.expression.human <- convertToHuman(
#' #    counts = bodyMap.mouse,
#' #    dictionary = ortholog.dict
#' #)
#'
#' #bsrdm <- BSRDataModel(
#' #    counts = matrix.expression.human,
#' #    species = "mmusculus",
#' #    conversion.dict = ortholog.dict
#' #)
#' 
#' #bsrdm <- learnParameters(bsrdm,
#' #    quick = TRUE  
#' #)
#' 
#' #reactSubset <- getResource(resourceName = "Reactome",
#' #cache = TRUE)
#' 
#' #subset <- c("REACTOME_BASIGIN_INTERACTIONS",
#' #"REACTOME_SYNDECAN_INTERACTIONS",
#' #"REACTOME_ECM_PROTEOGLYCANS",
#' #"REACTOME_CELL_JUNCTION_ORGANIZATION")
#' 
#' #reactSubset <- reactSubset[
#' #reactSubset$`Reactome name` %in% subset,]
#' 
#' #bsrinf.mouse <- BSRInference(bsrdm,reference="REACTOME")
#' 
#' bsrinf <- resetToInitialOrganism(bsrinf.mouse, 
#' conversion.dict = ortholog.dict)
setMethod("resetToInitialOrganism", "BSRInference", function(obj,
    conversion.dict) {
    # Need to check conversion.dict format

    conversion.dict$human.gene.name <- rownames(conversion.dict)

    LRinter(obj)$L <- .geneNameConversion(LRinter(obj)$L, conversion.dict)
    LRinter(obj)$R <- .geneNameConversion(LRinter(obj)$R, conversion.dict)

    ligands(obj) <- .geneNameConversion(ligands(obj), conversion.dict)
    receptors(obj) <- .geneNameConversion(receptors(obj), conversion.dict)
    tgGenes(obj) <- .geneNameConversion(tgGenes(obj), conversion.dict)

    obj
    
}) # resetToInitialOrganism


#' @title Convert gene symbols to another organism
#'
#' @description Convert gene symbols to another organism
#' based on a dictionary with human and orthologs in the other species.
#'
#' @param genes The genes you want to convert.
#' @param conversion.dict A data frame containing
#' gene names for the source species and Homo sapiens.
#'
#' @return Depend on the type of input genes
#' LRinter return a vector of genes
#' tgGenes receptors ligands: return list of list of genes
#' @keywords internal
.geneNameConversion <- function(genes, conversion.dict) {
    # print(".geneNameConversion")
    if (typeof(genes) == "character") {
        genes.df <- data.frame(human.gene.name = genes)

        genes.df$id <- seq_len(nrow(genes.df))

        genes.converted <- merge(genes.df, conversion.dict,
            by.x = "human.gene.name", sort = FALSE, all = FALSE
        )
        genes.converted <- genes.converted[order(genes.converted$id), ]

        genes.converted$human.gene.name <- NULL
        genes.converted$id <- NULL

        as.vector(unlist(genes.converted))
    } else if (typeof(genes) == "list") {
        list <- list()

        for (i in seq_len(length(genes))) {
            genes.df <- data.frame(human.gene.name = genes[[i]])
            genes.df$id <- seq_len(nrow(genes.df))
            genes.converted <- merge(genes.df, conversion.dict,
                by.x = "human.gene.name",
                sort = FALSE, all = FALSE
            )
            genes.converted <- genes.converted[order(genes.converted$id), ]

            genes.converted$human.gene.name <- NULL
            genes.converted$id <- NULL

            list[[i]] <- as.vector(unlist(genes.converted))
            rm(genes.df)
            rm(genes.converted)
        }

        list
    } else {
        stop("Something went wrong during gene conversion.", call. = FALSE)
    }
  
} # .geneNameConversion
