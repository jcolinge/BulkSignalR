#' BulkSignalR cluster comparison-based inference object
#'
#' An S4 class to represent ligand-receptor interactions inferred from
#' a comparison between two clusters of samples. This class inherits from
#' BSRInference.
#'
#' @slot cmp.name  The name of the BSRClusterComp object in a BSRDataModelComp
#' object comp list.
#' @slot src.cmp.name  The name of an optional BSRClusterComp object in a
#' BSRDataModelComp object comp list in case paracrine inferences were
#' performed.
#' @slot tg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot tg.logFC  A list of target gene logFC, one
#' entry per interaction
#' @slot tg.expr  A list of target gene expression, one
#' entry per interaction
#'
#' @details This class is contains inferred LR interactions along with
#' their statistical significance. Data representation supports subsequent
#' reductions to pathways, etc. See reduction functions
#' \code{"\link[=BSRInferenceComp-class]{reduceToBestPathway}"},
#' \code{"\link[=BSRInferenceComp-class]{reduceToLigand}"},
#' \code{"\link[=BSRInferenceComp-class]{reduceToReceptor}"} and
#' \code{"\link[=BSRInferenceComp-class]{reduceToPathway}"}.
#' @export
#' @examples
#' new("BSRInferenceComp")
setClass("BSRInferenceComp",
    contains = "BSRInference",
    slots = c(
        cmp.name = "character",
        src.cmp.name = "character",
        tg.pval = "list",
        tg.logFC = "list",
        tg.expr = "list"
    ),
    prototype = list(
        cmp.name = "myComparison",
        src.cmp.name = "",
        tg.pval = list(c(0.05, 0.1, 0.008)),
        tg.logFC = list(c(-1, 0, 2)),
        tg.expr = list(c(1, 2, 3))
    )
)


setValidity(
    "BSRInferenceComp",
    function(object) {
        if (!is.character(object@cmp.name)) {
            return("cmp.name is not of character type")
        }
        if (!is.character(object@src.cmp.name)) {
            return("src.cmp.name is not of character type")
        }
        if (length(object@cmp.name) == 0) {
            return("cmp.name must have a length > 0")
        }
        if (!is.list(object@tg.pval)) {
            return("tg.pval is not a list")
        }
        if (!is.list(object@tg.logFC)) {
            return("tg.logFC is not a list")
        }
        if (!is.list(object@tg.expr)) {
            return("tg.expr is not a list")
        }

        TRUE
    }
)


setMethod(
    "show", "BSRInferenceComp",
    function(object) {
        callNextMethod()
        cat("Cluster comparison name:", object@cmp.name, "\n")
        cat("Source cluster comparison name:", object@src.cmp.name, "\n")
    }
)


# Constructor ========================================================

#' Inference of ligand-receptor interactions based on regulation
#'
#' This method supports two configurations that we refer to
#' as paracrine and autocrine.
#'
#' In the autocrine case, a single cluster comparison name is provided.
#' In the corresponding cluster comparison, a group of samples A was
#' compared to a group of samples B to determine fold-changes and associated
#' P-values. The inferred ligand-receptor interactions take place in the
#' samples of group A. A typical single-cell example would be a population of
#' macrophages (group A) compared to all the other populations (group B) to
#' represent specific increased or decreased expression in macrophages. The
#' resulting ligand-receptor interactions will be autocrine interactions
#' that are exacerbated (or reduced depending on the chosen parameters) in
#' macrophages.
#'
#' In the paracrine case, two cluster comparison names must be provided.
#' For instance, a first comparison could involve macrophages versus all
#' the other cell populations as above. The second comparison could be
#' B-cells against all the other populations. Now, calling
#' \code{BSRInferenceComp()}
#' with comparison macrophages *versus* the rest and, as source comparison,
#' B-cells *versus* the rest, will result in inferring interactions between
#' B-cells (ligands) and macrophages (receptors and downstream pathways). To
#' obtain macrophages to B-cells paracrine interactions, it is necessary to call
#' the method a second time with permuted cluster comparison names. Another
#' example in spatial transcriptomics could be two thin bands at the boundary of
#' two tissue regions, one emitting the ligand and the other one expressing the
#' receptor.
#'
#' In this initial inference, all the receptor-containing pathways are reported,
#' see reduction functions to reduce this list.
#'
#' @name BSRInferenceComp
#'
#' @param obj       A BSRDataModelComp object.
#' @param cmp.name        The name of the cluster comparison that should be used
#' for the inference. Autocrine interactions if only this comparison name is
#' provided, paracrine if a source comparison name is provided as well.
#' @param src.cmp.name    The name of the source cluster comparison that should
#' be used for paracrine interaction inferences.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param max.pval        The maximum P-value imposed to both the ligand
#'   and the receptor.
#' @param min.logFC       The minimum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
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
#' @param min.t.logFC     The minimum log2 fold-change allowed for
#'   targets in case pos.targets or neg.targets are used.
#' @param pos.targets   A logical imposing that all the network targets must
#'   display positive logFC, i.e. logFC >= min.t.logFC.
#' @param neg.targets   A logical imposing that all the network targets must
#'   display negative logFC, i.e. logFC <= - min.t.logFC.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#' @param use.full.network  A logical to avoid limiting the reference network
#' to the detected genes and use the whole reference network.
#'
#' @details Perform the initial ligand-receptor inference. Initial means that
#' no reduction is applied. All the (ligand, receptor, downstream pathway)
#' triples are reported, i.e., a given LR pair may appear multiple times
#' with different pathways downstream the receptor. Specific reduction
#' functions are available from the package to operate subsequent
#' simplifications based on the BSRInferenceComp object created by this method.
#'
#' Here, ligand-receptor interactions are inferred based on gene or protein
#' regulation-associated P-values when comparing two clusters of samples. Since
#' a BSRDataModelComp object can contain several such comparisons, the name
#' of the comparison to use must be specified (parameter \code{cmp.name}).
#'
#' Note that since the introduction of the use.full.network parameter
#' (April 29, 2024), the pathway sizes are always computed before potential
#' intersection with the observed data (use.full.network set to FALSE) for
#' consistency. Accordingly, the minimum and maximum pathway default values
#' have been raised from 5 & 200 to 5 & 400 respectively. By default,
#' use.full.network is set to FALSE.
#'
#' In addition to statistical significance estimated according to BulkSignalR
#' statistical model, we compute SingleCellSignalR original LR-score,
#' based on L and R cluster average expression. 
#' In the paracrine case, L average expression
#' is taken from the source cluster.
#'
#' @return A BSRInferenceComp object with initial inferences set.
#'
#' @export
#'
#' @examples
#' data(bsrdm.comp, package = "BulkSignalR")
#' data(immune.signatures, package = "BulkSignalR")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf.comp <- BSRInferenceComp(bsrdm.comp, max.pval = 1, 
#' reference="REACTOME",
#' "random.example")
#' 
#' @importFrom methods new
BSRInferenceComp <- function(obj, cmp.name, 
    src.cmp.name = NULL, rank.p = 0.55,
    max.pval = 0.01, min.logFC = 1, neg.receptors = FALSE,
    pos.targets = FALSE, neg.targets = FALSE,
    min.t.logFC = 0.5, restrict.genes = NULL,
    use.full.network = FALSE,
    reference = c("REACTOME-GOBP", "REACTOME", "GOBP"),
    max.pw.size = 400, min.pw.size = 5, min.positive = 2,
    restrict.pw = NULL, with.complex = TRUE,
    fdr.proc = c(
        "BH", "Bonferroni", "Holm", "Hochberg",
        "SidakSS", "SidakSD", "BY", "ABH", "TSBH"
        )) {


    if (!(cmp.name %in% names(comparison(obj)))) {
        stop("cmp.name must exist in the names ",
            "of comparisons contained in obj")
    }
    if (!is.null(src.cmp.name) && !(src.cmp.name %in% names(comparison(obj)))) {
        stop("src.cmp.name must exist in the names",
        " of comparisons contained in obj")
    }
    reference <- match.arg(reference)
    fdr.proc <- match.arg(fdr.proc)
    if (min.logFC <= 0) {
        stop("min.logFC must be >0")
    }
    if (min.t.logFC <= 0) {
        stop("min.t.logFC must be >0")
    }
    if (rank.p < 0 || rank.p > 1) {
        stop("rank.p must lie in [0;1]")
    }
    if (neg.targets && pos.targets) {
        stop("neg.targets and pos.targets cannot be TRUE simultaneously")
    }

    # retrieve the BSRClusterComp object(s)
    cc <- comparison(obj)[[cmp.name]]
    if (!is.null(src.cmp.name)) {
        scc <- comparison(obj)[[src.cmp.name]]
    } else {
        scc <- NULL
    }

    # store inference parameters and retrieve relevant L & R
    inf.param <- list()
    inf.param$log.transformed.data <- logTransformed(obj)
    inf.param$mu <- mu(obj)
    inf.param$col.clusterA <- colClusterA(cc)
    inf.param$col.clusterB <- colClusterB(cc)
    if (is.null(src.cmp.name)) {
        inf.param$inference.type <- "autocrine"
    } else {
        inf.param$inference.type <- "paracrine"
        inf.param$src.colA <- colClusterA(scc)
        inf.param$src.colB <- colClusterB(scc)
    }
    inf.param$max.pval <- max.pval
    inf.param$min.logFC <- min.logFC
    inf.param$neg.receptors <- neg.receptors
    inf.param$pos.targets <- pos.targets
    inf.param$neg.targets <- neg.targets
    inf.param$min.t.logFC <- min.t.logFC
    inf.param$restrict.genes <- restrict.genes
    inf.param$use.full.network <- use.full.network
    lr <- .getRegulatedLR(obj, cc, scc,
        max.pval = max.pval, min.logFC = min.logFC,
        neg.receptors = neg.receptors, restrict.genes = restrict.genes
    )

    # apply BSR model on the targets
    inf.param$reference <- reference
    inf.param$min.pw.size <- min.pw.size
    inf.param$max.pw.size <- max.pw.size
    inf.param$with.complex <- with.complex
    inf.param$min.positive <- min.positive
    inf.param$restrict.pw <- restrict.pw
    pairs <- .checkRegulatedReceptorSignaling(obj, cc, lr,
        reference = reference,
        pos.targets = pos.targets, neg.targets = neg.targets,
        min.t.logFC = min.logFC, use.full.network = use.full.network,
        min.pw.size = min.pw.size, max.pw.size = max.pw.size,
        min.positive = min.positive, with.complex = with.complex,
        restrict.pw = restrict.pw
    )
    inf.param$fdr.proc <- fdr.proc
    inf.param$rank.p <- rank.p

    # compute P-values
    inter <- .pValuesRegulatedLR(pairs, 
        parameters(obj), rank.p = rank.p, fdr.proc = fdr.proc)

    # compute LR-score for compatibility with SingleCellSignalR version 1
    if (is.null(scc)) {
        inter$L.expr <- differentialStats(cc)[inter$L, "expr"]
    } else {
        inter$L.expr <- differentialStats(scc)[inter$L, "expr"]
    }
    inter$R.expr <- differentialStats(cc)[inter$R, "expr"]
	if (sum(inter$L.expr<0) + sum(inter$R.expr<0) > 0){
	    # impossible to compute LR-score with negative values,
		# this may happen for scProt-MS data
		inter$LR.score <- 0
	}
	else{
	    # normal case
        if (inf.param$log.transformed.data) {
            sq <- sqrt(inter$L.expr * inter$R.expr)
        } else {
            sq <- sqrt(log1p(inter$L.expr) / log(2) * log1p(inter$R.expr) / log(2))
        }
        inter$LR.score <- sq / (inf.param$mu + sq)
	}

    # prepare the accompanying lists
    ligands <- strsplit(inter$L, ";")
    receptors <- strsplit(inter$R, ";")
    tg <- strsplit(inter$target.genes, ";")
    tgpval <- lapply(
        strsplit(inter$target.pval, ";"),
        function(x) as.numeric(x)
    )
    tglogFC <- lapply(
        strsplit(inter$target.logFC, ";"),
        function(x) as.numeric(x)
    )
    tgcorr <- lapply(
        strsplit(inter$target.corr, ";"),
        function(x) as.numeric(x)
    )
    tgexpr <- lapply(
        strsplit(inter$target.expr, ";"),
        function(x) as.numeric(x)
    )
    inf.param$ligand.reduced <- FALSE
    inf.param$receptor.reduced <- FALSE
    inf.param$pathway.reduced <- FALSE

    # instantiate the object
    if (is.null(src.cmp.name)) {
        src.cmp.name.char <- ""
    } else {
        src.cmp.name.char <- src.cmp.name
    }

    new("BSRInferenceComp",
        LRinter = inter[, c(
            "L", "R", "pw.id", "pw.name", "pval", "qval",
            "L.logFC", "R.logFC", "LR.pval", "LR.corr",
            "rank", "len", "rank.pval", "rank.corr",
            "LR.score", "L.expr", "R.expr"
        )],
        ligands = ligands, 
        receptors = receptors, 
        tg.genes = tg, 
        tg.corr = tgcorr, inf.param = inf.param,
        tg.pval = tgpval, tg.logFC = tglogFC, tg.expr = tgexpr, 
        cmp.name = cmp.name, src.cmp.name = src.cmp.name.char
    )
    
} # BSRInferenceComp


# Accessors & setters ========================================================

setGeneric("comparisonName", signature="x",
    function(x) standardGeneric("comparisonName")
)
#' Comparison name accessor
#'
#' @name comparisonName
#' @aliases comparisonName,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return cmp.name
#' @examples
#' bsrinf <- new("BSRInferenceComp")
#' comparisonName(bsrinf)
#' @export
setMethod("comparisonName", "BSRInferenceComp", function(x) x@cmp.name)


setGeneric("comparisonName<-", signature=c("x", "value"),
    function(x, value) standardGeneric("comparisonName<-")
)
#' Comparison name setter (internal use only)
#' 
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("comparisonName<-", "BSRInferenceComp", function(x, value) {
    x@cmp.name <- value
    methods::validObject(x)
    x
})


setGeneric("sourceComparisonName", signature="x",
    function(x) standardGeneric("sourceComparisonName")
)
#' Source comparison name accessor
#'
#' @name sourceComparisonName
#' @aliases sourceComparisonName,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return src.comp.name
#' @examples
#' bsrinf <- new("BSRInferenceComp")
#' sourceComparisonName(bsrinf)
#' @export
setMethod("sourceComparisonName", "BSRInferenceComp", 
    function(x) x@src.cmp.name)


setGeneric("sourceComparisonName<-", signature=c("x", "value"),
    function(x, value) standardGeneric("sourceComparisonName<-")
)
#' Source comparison name setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("sourceComparisonName<-", "BSRInferenceComp", function(x, value) {
    x@src.cmp.name <- value
    methods::validObject(x)
    x
})


setGeneric("tgPval", signature="x",
    function(x) standardGeneric("tgPval")
)

#' Target gene P-values accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return tgPval
#' @examples
#' bsrinf <- new("BSRInferenceComp")
#' tgPval(bsrinf)
#' @export
setMethod("tgPval", "BSRInferenceComp", function(x) x@tg.pval)


setGeneric("tgPval<-", signature=c("x", "value"),
    function(x, value) standardGeneric("tgPval<-")
)
#' Target gene P-values setter (internal use only)
#'
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("tgPval<-", "BSRInferenceComp", function(x, value) {
    x@tg.pval <- value
    methods::validObject(x)
    x
})


setGeneric("tgLogFC", signature="x",
    function(x) standardGeneric("tgLogFC")
)
#' Target gene logFC accessor
#'
#' @name tgLogFC
#' @aliases tgLogFC,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return tgLogFC
#' @examples
#' bsrinf <- new("BSRInferenceComp")
#' tgLogFC(bsrinf)
#' @export
setMethod("tgLogFC", "BSRInferenceComp", function(x) x@tg.logFC)


setGeneric("tgLogFC<-", signature=c("x", "value"),
    function(x, value) standardGeneric("tgLogFC<-")
)
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("tgLogFC<-", "BSRInferenceComp", function(x, value) {
    x@tg.logFC <- value
    methods::validObject(x)
    x
})


setGeneric("tgExpr", signature="x",
    function(x) standardGeneric("tgExpr")
)
#' Target gene expression accessor
#'
#' @name tgExpr
#' @aliases tgExpr,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return tgExpr
#' @examples
#' bsrinf <- new("BSRInferenceComp")
#' tgExpr(bsrinf)
#' @export
setMethod("tgExpr", "BSRInferenceComp", function(x) x@tg.expr)


setGeneric("tgExpr<-", signature=c("x", "value"),
    function(x, value) standardGeneric("tgExpr<-")
)
#' Target gene expression setter (internal use only)
#'
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @return returns \code{NULL}
#' @keywords internal
setMethod("tgExpr<-", "BSRInferenceComp", function(x, value) {
    x@tg.expr <- value
    methods::validObject(x)
    x
})


# simplified table views =======================================================

#' Simplified LRinter accessor reporting the essential columns
#'
#' @name LRinterShort
#' @aliases LRinterShort,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return LRinterShort
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' LRinterShort(bsrinf.comp)[5,]
#' @export
setMethod(
    "LRinterShort", "BSRInferenceComp",
    function(x) {
        x@LRinter[, c(
            "L", "R", "pw.id", "pw.name",
            "qval", "L.logFC", "R.logFC", "len"
        )]
    }
)


setGeneric("LRinterScore", signature="x",
    function(x) standardGeneric("LRinterScore")
)
#' Simplified LRinter accessor with focus on the LR-score
#'
#' @name LRinterScore
#' @aliases LRinterScore,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @return LRinterScore
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' LRinterScore(bsrinf.comp)[5,]
#' @export
setMethod(
    "LRinterScore", "BSRInferenceComp",
    function(x) {
        x@LRinter[, c(
            "L", "R", "pw.id", "pw.name",
            "LR.score", "L.expr", "R.expr", "len"
        )]
    }
)


# Rescoring & updating =========================================================

#' Inference re-scoring
#'
#' A method to re-score an existing BSRInferenceComp object
#' (P- and Q-value estimations).
#'
#' @name rescoreInference
#' @aliases rescoreInference,BSRInferenceComp-method
#'
#' @param obj BSRInferecenceComp object.
#' @param param NULL by default
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @details A BSRInferenceComp object should be created by calling
#' \code{"\link[=BSRInferenceComp-class]{BSRInferenceComp}"}
#'
#' @return A BSRInferenceComp object.
#'
#' @export
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#'
#' bsrinf.less <- rescoreInference(bsrinf.comp, 
#' rank.p = 0.75)
setMethod("rescoreInference", "BSRInferenceComp", function(obj,
    param = NULL, 
    rank.p = 0.55,
    fdr.proc = c("BH", "Bonferroni", "Holm",
    "Hochberg", "SidakSS", "SidakSD", "BY", "ABH", "TSBH"
    )) {

    if (rank.p < 0 || rank.p > 1) {
        stop("rank.p must lie in [0;1]")
    }
    fdr.proc <- match.arg(fdr.proc)

    # extract the necessary data from the BSRInferenceComp object
    pairs <- LRinter(obj)
    tg.genes <- tgGenes(obj)
    tg.pval <- tgPval(obj)
    tg.corr <- tgCorr(obj)

    # recompute P-values, LR- and LRT-scores
    for (i in seq_len(nrow(pairs))) {
        tg <- tg.genes[[i]]
        spvals <- tg.pval[[i]]

        # get the LR correlation P-value
        p.lr <- pairs$LR.pval[i]

        # estimate the target gene correlation P-value based on rank statistics
        # for the individual correlation Gaussian model
        len <- pairs$len[i]
        r <- min(max(1, trunc(rank.p * len)), len)
        pairs$rank[i] <- r
        rank.pval <- spvals[r]
        # r-1 P-values are > rank.pval, prob to have r-1 or less
        # P-values > rank.pval is given by a binomial with success rate
        # equal to the probability to get a P-value > rank.pval, i.e.,
        # 1-rank.pval. If rank.pval is low (i.e., highly significant),
        # it becomes difficult to get as little 
        # as r-1 P-values > rank.pval by chance!
        p.rt <- stats::pbinom(r - 1, len, 1 - rank.pval) # cdf is punif here!
        pairs$pval[i] <- p.lr * p.rt
        pairs$rank.pval[i] <- rank.pval
        pairs$rank.corr[i] <- tg.corr[[i]][r]
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


setGeneric("updateInference", signature="obj",
    function(obj,...) standardGeneric("updateInference")
)
#' Inference updating
#'
#' A method to update the data underlying statistical significance estimations
#' prior to rescoring for an existing BSRInferenceComp object
#' (P- and Q-value estimations as well as LR-score).
#'
#' @name updateInference
#' @aliases updateInference,BSRInferenceComp-method
#'
#' @param obj BSRInferenceComp object.
#' @param bsrcc BSRClusterComp object relative to target cells.
#' @param ncounts Matrix counts normalized.
#' @param src.bsrcc BSRClusterComp object relative to source cells.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param max.pval        The maximum P-value imposed to both the ligand
#'   and the receptor.
#' @param min.logFC       The minimum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param min.LR.score       The minimum LR-score allowed for the interaction.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
#' @param fdr.proc      The procedure for adjusting P-values according to
#' \code{\link[multtest]{mt.rawp2adjp}}.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param min.t.logFC     The minimum log2 fold-change allowed for
#'   targets in case pos.targets or neg.targets are used.
#' @param pos.targets   A logical imposing that all the network targets must
#'   display positive logFC, i.e. logFC >= min.t.logFC.
#' @param neg.targets   A logical imposing that all the network targets must
#'   display negative logFC, i.e. logFC <= - min.t.logFC.
#'
#' @details A BSRInferenceComp object should be created by calling
#' \code{"\link[=BSRInferenceComp-class]{BSRInferenceComp}"}
#'
#' @return A BSRInferenceComp object. The main application of this
#' method is to take a "universal" inference obtained by assigning
#' each gene to good logFC, P-values and expression levels whose role
#' is to find all the reachable targets per receptor/pathway, and
#' to update it by using actual logFC, P-values, and expression data.
#' The benefit is to save time when multiple sample comparisons are
#' performed, only one network exploration is necessary. Note that
#' if a restrictive logic such as \code{positive.targets=TRUE} is used,
#' the result will be correct provided all the targets were in the
#' initial BSRInferenceComp object. If a restriction on the targets
#' was applied, then the update is likely to miss some targets, i.e.,
#' the statistical analysis will be wrong.
#'
#' In case no L-R interaction is above the required thresholds, the
#' value `NULL` is returned.
#'
#' Note that correlations are set to 1 to avoid
#' lengthy computations with scRNA-seq data and multiple cell
#' populations.
#'
#' The main role of this method is to support our SingleCellSignalR Version 2
#' package.
#'
#' @export
#' @examples
#' data(bsrdm.comp, package = "BulkSignalR")
#' data(bsrinf.comp, package = "BulkSignalR")
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' 
#' #bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
#'
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval = runif(n),
#' logFC = rnorm(n, 0, 2),
#' expr = runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' 
#' 
#' # update
#' stats$pval <- stats$pval / 100
#' stats$logFC <- stats$logFC + 0.5
#' 
#' bsrcc.2 <- BSRClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrinf.updated <- updateInference(bsrinf.comp, bsrcc.2,
#' max.pval = 1, min.logFC = 0.1)
setMethod("updateInference", "BSRInferenceComp", function(obj, bsrcc,
    ncounts, src.bsrcc = NULL, rank.p = 0.55,
    max.pval = 0.01, min.logFC = 1, min.LR.score = 0,
    neg.receptors = FALSE,
    pos.targets = FALSE, neg.targets = FALSE,
    min.t.logFC = 0.5, min.positive = 2,
    fdr.proc = c(
        "BH", "Bonferroni", "Holm", "Hochberg",
        "SidakSS", "SidakSD", "BY", "ABH", "TSBH"
        )) {
    if (!is(obj, "BSRInferenceComp")) {
        stop("obj must be a BSRInferenceComp object")
    }
    if (!is(bsrcc, "BSRClusterComp")) {
        stop("bsrcc must be a BSRClusterComp object")
    }
    if (!is.null(src.bsrcc) && !is(src.bsrcc, "BSRClusterComp")) {
        stop("src.bsrcc must be a BSRClusterComp object")
    }

    # get BSRInferenceComp object details
    inter <- LRinter(obj)
    L <- ligands(obj)
    R <- receptors(obj)
    t <- tgGenes(obj)
    c <- tgCorr(obj)
    lfc <- tgLogFC(obj)
    p <- tgPval(obj)
    e <- tgExpr(obj)
    par <- inferenceParameters(obj)
    mu <- par$mu
    logTransf <- par$log.transformed.data
    R.stats <- differentialStats(bsrcc)
    if (!is.null(src.bsrcc)) {
        L.stats <- differentialStats(src.bsrcc)
    } # paracrine
    else {
        L.stats <- R.stats
    } # autocrine

    # update parameters
    par$col.clusterA <- colClusterA(bsrcc)
    par$col.clusterB <- colClusterB(bsrcc)
    if (is.null(src.bsrcc)) {
        par$inference.type <- "autocrine"
    } else {
        par$inference.type <- "paracrine"
        par$src.colA <- colClusterA(src.bsrcc)
        par$src.colB <- colClusterB(src.bsrcc)
    }
    par$max.pval <- max.pval
    par$min.logFC <- min.logFC
    par$min.LR.score <- min.LR.score
    par$neg.receptors <- neg.receptors
    par$pos.targets <- pos.targets
    par$neg.targets <- neg.targets
    par$min.t.logFC <- min.t.logFC
    par$min.positive <- min.positive
    par$fdr.proc <- fdr.proc
    par$rank.p <- rank.p

    # assign correct logFC, P-values, correlations, and expression to L and R
    inter$L.logFC <- L.stats[inter$L, "logFC"]
    inter$R.logFC <- R.stats[inter$R, "logFC"]
    inter$LR.pval <- L.stats[inter$L, "pval"] * R.stats[inter$R, "pval"]
    inter$L.expr <- L.stats[inter$L, "expr"]
    inter$R.expr <- R.stats[inter$R, "expr"]

    # if (is.null(src.bsrcc)){
    #   corlr <- stats::cor(t(ncounts[, c(colA(bsrcc),colB(bsrcc))]),
    # method = "spearman")
    # }
    # else {
    #   corlr <- stats::cor(t(ncounts[, c(colA(bsrcc),colA(src.bsrcc))]),
    # method = "spearman")
    #}
    # for (i in seq_len(nrow(inter)))
    #   inter$LR.corr[i] <- corlr[inter$L[i], inter$R[i]]
    inter$LR.corr <- 1

    # LR-score
	if (sum(inter$L.expr<0) + sum(inter$R.expr<0) > 0){
	    # impossible to compute LR-score with negative values,
		# this may happen for scProt-MS data
		inter$LR.score <- 0
	}
	else{
	    # normal case
        if (logTransf) {
            sq <- sqrt(inter$L.expr * inter$R.expr)
        } else {
            sq <- sqrt(log1p(inter$L.expr) / log(2) *
			           log1p(inter$R.expr) / log(2))
        }
        inter$LR.score <- sq / (mu + sq)
	}

    # select on L & R as well as LR-scores
    good <- L.stats[inter$L, "pval"] <= max.pval & inter$L.logFC >= min.logFC &
        R.stats[inter$R, "pval"] <= max.pval
    if (neg.receptors) {
        good <- good & abs(inter$R.logFC) >= min.logFC
    } else {
        good <- good & inter$R.logFC >= min.logFC
    }
    good <- good & inter$LR.score >= min.LR.score
    if (sum(good) == 0) {
        return(NULL)
    }
    inter <- inter[good, ]
    L <- L[good]
    R <- R[good]
    t <- t[good]
    c <- c[good]
    lfc <- lfc[good]
    p <- p[good]
    e <- e[good]

    # assign correct logFC, P-values, 
    # correlations, and expression to the targets
    # and select the targets
    keep <- NULL
    for (i in seq_len(nrow(inter))) {
        genes <- t[[i]]
        logfc <- R.stats[genes, "logFC"]
        if (pos.targets || neg.targets) {
            if (pos.targets) {
                genes <- genes[logfc >= min.t.logFC]
            } else {
                genes <- genes[logfc <= -min.t.logFC]
            }
        } else {
            genes <- genes[abs(logfc) >= min.t.logFC]
        }
        if (length(genes) < min.positive) {
            keep <- c(keep, FALSE)
        } else {
            # if all conditions are met, list all target genes with
            # their regulation P-values in a data frame
            # row. Target genes are sorted wrt P-values in decreasing
            # order to keep the compatibility with correlation analysis,
            # where the most significant values are at the end.
            pv <- R.stats[genes, "pval"]
            o <- order(pv, decreasing = TRUE)
            pv <- pv[o]
            logfc <- R.stats[genes, "logFC"]
            logfc <- logfc[o]
            expr <- R.stats[genes, "expr"]
            expr <- expr[o]
            genes <- genes[o]
            co <- rep(1, length(genes)) # corlr[inter[i, "R"], genes]
            t[[i]] <- genes
            c[[i]] <- co
            lfc[[i]] <- logfc
            p[[i]] <- pv
            e[[i]] <- expr
            len <- length(genes)
            inter[i, "len"] <- len
            rank <- min(max(1, trunc(rank.p * len)), len)
            inter[i, "rank"] <- rank
            # rank.corr and rank.pval are updated by
            #  calling rescoreInference() below
            keep <- c(keep, TRUE)
        }
    }
    if (sum(keep) == 0) {
        return(NULL)
    }
    inter <- inter[keep, ]
    L <- L[keep]
    R <- R[keep]
    t <- t[keep]
    c <- c[keep]
    lfc <- lfc[keep]
    p <- p[keep]
    e <- e[keep]

    # update object
    LRinter(obj) <- inter
    ligands(obj) <- L
    receptors(obj) <- R
    tgGenes(obj) <- t
    tgCorr(obj) <- c
    tgPval(obj) <- p
    tgLogFC(obj) <- lfc
    tgExpr(obj) <- e
    inferenceParameters(obj) <- par

    # rescore and return
    rescoreInference(obj, par, rank.p = rank.p, fdr.proc = fdr.proc)
    
}) # updateInferenceComp


# Reduction and pathway stat methods ===========================================

#' Keep one pathway per ligand-receptor pair
#'
#' @name reduceToBestPathway
#' @aliases reduceToBestPathway,BSRInferenceComp-method
#'
#' @param obj BSRInferenceComp object
#'
#' @return A BSRInferenceComp object reduced to only report one pathway per
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
#' data(bsrinf.comp, package = "BulkSignalR")
#' 
#' 
#' reduceToBestPathway(bsrinf.comp)
#'
#' @importFrom rlang .data
setMethod("reduceToBestPathway", "BSRInferenceComp", function(obj) {
    # Here we access the object slots directly as this procedure
    # is dependent of actual data representation

    # get best p-value pathway per LR pair
    ligands <- list()
    receptors <- list()
    tg.genes <- list()
    tg.corr <- list()
    tg.pval <- list()
    tg.logFC <- list()
    tg.expr <- list()

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
        tg.pval <- c(tg.pval, obj@tg.pval[j])
        tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
        tg.expr <- c(tg.expr, obj@tg.expr[j])
        LRinter <- rbind(LRinter, pairs[j, ])
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@tg.pval <- tg.pval
    obj@tg.logFC <- tg.logFC
    obj@tg.expr <- tg.expr
    obj@inf.param$pathway.reduced <- TRUE

    obj
    
}) # reduceToBestPathway


#' Aggregate the ligands of a same receptor
#'
#' Simplifies a ligand-receptor table to focus on the receptors.
#'
#' @name reduceToReceptor
#' @aliases reduceToReceptor,BSRInferenceComp-method
#'
#' @return BSRInferenceComp object reduced to one row per receptor.
#' All the ligands are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the line with the
#' pathway featuring the smallest P-value.
#' The same logic applies to the LR-score, and the ligand
#' expression.
#' @param obj BRSInferenceComp object
#' @export
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' # reduction
#' bsrinf.redR <- reduceToReceptor(bsrinf.comp)
#' @importFrom rlang .data
setMethod("reduceToReceptor", "BSRInferenceComp", function(obj) {
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
    tg.pval <- list()
    tg.logFC <- list()
    tg.expr <- list()
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
        tg.pval <- c(tg.pval, obj@tg.pval[j])
        tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
        tg.expr <- c(tg.expr, obj@tg.expr[j])
        to.add <- pairs[j, ]
        to.add[1, "L"] <- paste0("{", paste(unique(lig$L), collapse = ";"), "}")
        to.add[1, "LR.score"] <- max(lig$LR.score)
        to.add[1, "L.expr"] <- max(lig$L.expr)
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@tg.pval <- tg.pval
    obj@tg.logFC <- tg.logFC
    obj@tg.expr <- tg.expr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$ligand.reduced <- TRUE

    obj
    
}) # reduceToReceptor


#' Aggregate the receptors of a same ligand
#'
#' Simplifies a ligand-receptor table to focus on the ligands.
#'
#' @name reduceToLigand
#' @aliases reduceToLigand,BSRInferenceComp-method
#'
#' @return A BSRInferenceComp object but reduced to one row per ligand.
#' All the receptors are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the pathway with
#' the smallest P-value.
#' The same logic applies to the LR-score, and the receptor
#' expression.
#' @param obj BSRInferenceComp object
#' @export
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' 
#' bsrinf.redL <- reduceToLigand(bsrinf.comp)
#' @importFrom rlang .data
setMethod("reduceToLigand", "BSRInferenceComp", function(obj) {
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
    tg.pval <- list()
    tg.logFC <- list()
    tg.expr <- list()
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
        tg.pval <- c(tg.pval, obj@tg.pval[j])
        tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
        tg.expr <- c(tg.expr, obj@tg.expr[j])
        to.add <- pairs[j, ]
        to.add[1, "R"] <- paste0("{", paste(unique(rec$R), collapse = ";"), "}")
        to.add[1, "LR.score"] <- max(rec$LR.score)
        to.add[1, "R.expr"] <- max(rec$R.expr)
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@tg.pval <- tg.pval
    obj@tg.logFC <- tg.logFC
    obj@tg.expr <- tg.expr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj
    
}) # reduceToLigand


#' Aggregate ligands and receptors at the pathway level
#'
#' Simplifies a ligand-receptor inference object to focus on
#' the pathways.
#'
#' @name reduceToPathway
#' @aliases reduceToPathway,BSRInferenceComp-method
#'
#' @return A BSRInferenceComp object reduced to only report one row per pathway.
#' The information of which ligand interacted with which receptor is lost as
#' all the ligands and all the receptors forming pairs related to a certain
#' pathway are combined.
#' For a given pathway, the reported P-values and target genes are those of
#' the best ligand-receptor pair that was in this pathway.
#' The same logic applies to the LR-score, and the ligand and receptor
#' expression.
#' Receptors and ligands are combined in two semi-colon-separated
#' lists surrounded by curly brackets in the tabular slot \code{LRinter},
#' while the list representation slots (\code{ligands} and
#' \code{receptors}) are update accordingly.
#'
#' @param obj BSRInferenceComp object
#' @export
#' @examples
#' data(bsrinf.comp, package = "BulkSignalR")
#' 
#' bsrinf.redP <- reduceToPathway(bsrinf.comp)
#' @importFrom rlang .data
setMethod("reduceToPathway", "BSRInferenceComp", function(obj) {
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
    tg.pval <- list()
    tg.logFC <- list()
    tg.expr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (p in unique(pairs$pw.id)) {
        j <- which(pairs$pw.id == p)[1]
        ligands <- c(ligands, list(unique(pairs$L[pairs$pw.id == p])))
        receptors <- c(receptors, list(unique(pairs$R[pairs$pw.id == p])))
        tg.genes <- c(tg.genes, obj@tg.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        tg.pval <- c(tg.pval, obj@tg.pval[j])
        tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
        tg.expr <- c(tg.expr, obj@tg.expr[j])
        to.add <- pairs[j, ]
        to.add[1, "L"] <- paste0("{", paste(unique(pairs$L[pairs$pw.id == p]),
            collapse = ";"
        ), "}")
        to.add[1, "R"] <- paste0("{", paste(unique(pairs$R[pairs$pw.id == p]),
            collapse = ";"
        ), "}")
        to.add[1, "LR.score"] <- max(pairs$LR.score[pairs$pw.id == p])
        to.add[1, "R.expr"] <- max(pairs$R.expr[pairs$pw.id == p])
        to.add[1, "L.expr"] <- max(pairs$L.expr[pairs$pw.id == p])
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@tg.genes <- tg.genes
    obj@tg.corr <- tg.corr
    obj@tg.pval <- tg.pval
    obj@tg.logFC <- tg.logFC
    obj@tg.expr <- tg.expr
    obj@inf.param$ligand.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj
    
}) # reduceToPathway
