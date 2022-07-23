library(methods)

#' BulkSignalR Inference Object
#'
#' An S4 class to represent inferred ligand-receptor interactions.
#'
#' @slot LRinter  A data frame describing the (ligand,receptor,pathway) triples
#' with P- and Q-values.
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot t.genes  A list of target genes, one entry per LR interaction.
#' @slot tg.corr  A list of target gene correlations to the receptor, one
#' entry per interaction
#' @slot inf.param  The parameters used for the inference.
#'
#' @details This class is a container for inferred LR interactions along with
#' their statistical confidence. Data representation supports subsequent
#' reductions to pathways, etc. See reduction functions
#' \code{\link{reduceToBestPathway}}, \code{\link{reduceToLigand}},
#' \code{\link{reduceToReceptor}}, and \code{\link{reduceToPathway}}.
#' @export
#' @examples
#' new("BSRInference")
#'
setClass("BSRInference",
         slots=c(LRinter="data.frame",
                 ligands="list",
                 receptors="list",
                 t.genes="list",
                 tg.corr="list",
                 inf.param="list"),
         prototype=list(
             LRinter=data.frame(L="A", R="B", LR.corr=0.6, pw.id="123",
                                pw.name="one pw", rank=2, len=50,
                                rank.corr=0.6, target.genes="a;b;c",
                                target.corr="-0.5;0.1;0.8",
                                pval=1.0, qval=1.0,
                                stringsAsFactors=FALSE),
             ligands=list("A"),
             receptors=list("B"),
             t.genes=list(c("a","b","c")),
             tg.corr=list(c(-0.5,0.1,0.8)),
             inf.param=list()
         ))

setValidity("BSRInference",
            function(object) {
                if(!is.data.frame(object@LRinter))
                    return("LRinter is not a data frame")
                if(!is.list(object@ligands))
                    return("ligands is not a list")
                if(!is.list(object@receptors))
                    return("receptors is not a list")
                if(!is.list(object@t.genes))
                    return("t.genes is not a list")
                if(!is.list(object@tg.corr))
                    return("tg.corr is not a list")
                if (!is.list(object@inf.param))
                    return("inf.param is not a list")

                TRUE
            }
)

setMethod("show", "BSRInference",
          function(object) {
              cat("Reference database: ", object@inf.param$reference, "\n", sep="")
               print(as.data.frame(object@LRinter[order(object@LRinter$qval),
                        c("L", "R", "pval", "qval", "pw.id", "pw.name"),]))[5,]
              cat("Inference parameters:\n")
              print(object@inf.param)
          }
)


# Accessors & setters ========================================================

if (!isGeneric("LRinter")) {
    if (is.function("LRinter"))
        fun <- LRinter
    else
        fun <- function(x) standardGeneric("LRinter")
    setGeneric("LRinter", fun)
}
#' LRinter accessor
#' @export
setMethod("LRinter", "BSRInference", function(x) x@LRinter)

if (!isGeneric("LRinter<-")) {
    if (is.function("LRinter<-"))
        fun <- `LRinter<-`
    else
        fun <- function(x,value) standardGeneric("LRinter<-")
    setGeneric("LRinter<-", fun)
}
#' LRinter setter (internal use only)
setMethod("LRinter<-", "BSRInference", function(x, value){
    x@LRinter <- value
    methods::validObject(x)
    x
})

if (!isGeneric("ligands")) {
    if (is.function("ligands"))
        fun <- ligands
    else
        fun <- function(x) standardGeneric("ligands")
    setGeneric("ligands", fun)
}
#' ligands accessor
#' @export
setMethod("ligands", "BSRInference", function(x) x@ligands)

if (!isGeneric("ligands<-")) {
    if (is.function("ligands<-"))
        fun <- `ligands<-`
    else
        fun <- function(x,value) standardGeneric("ligands<-")
    setGeneric("ligands<-", fun)
}
#' ligands setter (internal use only)
setMethod("ligands<-", "BSRInference", function(x, value){
    x@ligands <- value
    methods::validObject(x)
    x
})

if (!isGeneric("receptors")) {
    if (is.function("receptors"))
        fun <- receptors
    else
        fun <- function(x) standardGeneric("receptors")
    setGeneric("receptors", fun)
}
#' receptors accessor
#' @export
setMethod("receptors", "BSRInference", function(x) x@receptors)

if (!isGeneric("receptors<-")) {
    if (is.function("receptors<-"))
        fun <- `receptors<-`
    else
        fun <- function(x, value) standardGeneric("receptors<-")
    setGeneric("receptors<-", fun)
}
#' receptors setter (internal use only)
setMethod("receptors<-", "BSRInference", function(x, value){
    x@receptors <- value
    methods::validObject(x)
    x
})

if (!isGeneric("tGenes")) {
    if (is.function("tGenes"))
        fun <- tGenes
    else
        fun <- function(x) standardGeneric("tGenes")
    setGeneric("tGenes", fun)
}
#' Target genes accessor
#' @export
setMethod("tGenes", "BSRInference", function(x) x@t.genes)

if (!isGeneric("tGenes<-")) {
    if (is.function("tGenes<-"))
        fun <- `tGenes<-`
    else
        fun <- function(x,value) standardGeneric("tGenes<-")
    setGeneric("tGenes<-", fun)
}
#' Target genes setter (internal use only)
setMethod("tGenes<-", "BSRInference", function(x, value){
    x@t.genes <- value
    methods::validObject(x)
    x
})

if (!isGeneric("tgCorr")) {
    if (is.function("tgCorr"))
        fun <- tgCorr
    else
        fun <- function(x) standardGeneric("tgCorr")
    setGeneric("tgCorr", fun)
}
#' Target gene correlations accessor
#' @export
setMethod("tgCorr", "BSRInference", function(x) x@tg.corr)

if (!isGeneric("tgCorr<-")) {
    if (is.function("tgCorr<-"))
        fun <- `tgCorr<-`
    else
        fun <- function(x,value) standardGeneric("tgCorr<-")
    setGeneric("tgCorr<-", fun)
}
#' Target gene correlations setter (internal use only)
setMethod("tgCorr<-", "BSRInference", function(x, value){
    x@tg.corr <- value
    methods::validObject(x)
    x
})

if (!isGeneric("infParam")) {
    if (is.function("infParam"))
        fun <- infParam
    else
        fun <- function(x) standardGeneric("infParam")
    setGeneric("infParam", fun)
}
#' Inference parameters accessor
#' @export
setMethod("infParam", "BSRInference", function(x) x@inf.param)
if (!isGeneric("infParam<-")) {
    if (is.function("infParam<-"))
        fun <- `infParam<-`
    else
        fun <- function(x,value) standardGeneric("infParam<-")
    setGeneric("infParam<-", fun)
}
#' Inference parameters setter (internal use only)
setMethod("infParam<-", "BSRInference", function(x, value){
    x@inf.param <- value
    methods::validObject(x)
    x
})


# Rescoring ====================================================================

if (!isGeneric("rescoreInference")) {
    if (is.function("rescoreInference"))
        fun <- rescoreInference
    else
        fun <- function(obj, ...) standardGeneric("rescoreInference")
    setGeneric("rescoreInference", fun)
}
#' Inference re-scoring
#'
#' A method to re-score an existing BSRInference object
#' (P- and Q-value estimations).
#'
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @details A BSRInference object should be created by calling
#' \code{\link{initialInference}}. Parameters controlling the estimation
#' of the statistical significance of the ligand/receptor/pathway triples
#' are provided at the time of calling the latter method.
#'
#' Nonetheless, it
#' might be useful to change the initially-provided parameters, in
#' which case this method should be called.
#'
#' @return A BSRInference object.
#'
#' @export
#' @examples
#' print('rescoreInference')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.new <- rescoreInference(bsrinf, param=param(bsrdm), rank.p=0.75)
#'
setMethod("rescoreInference", "BSRInference", function(obj, param, rank.p=0.55,
                    fdr.proc=c("BH", "Bonferroni", "Holm",
                    "Hochberg", "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {

    if (rank.p < 0 || rank.p > 1)
        stop("rank.p must lie in [0;1]")
    fdr.proc <- match.arg(fdr.proc)

    # extract the necessary data from the BSRInference object
    pairs <- LRinter(obj)
    t.genes <- tGenes(obj)
    tg.corr <- tgCorr(obj)

    # prepare the chosen model CDF
    LR.par <- param$LR.0$model
    RT.par <- param$RT.0$model
    if (LR.par$distrib != RT.par$distrib)
        stop("Distinct statistical models for LR and RT nulls are not allowed")
    if (LR.par$distrib == 'censored_normal')
        cdf <- .cdfGaussian
    else if (LR.par$distrib == 'censored_mixed_normal')
        cdf <- .cdfMixedGaussian
    else if (LR.par$distrib == 'censored_stable')
        cdf <- .cdfAlphaStable
    else
        stop(paste0("Unknown statistical model: ", LR.par$LR.0$model$distrib))

    # recompute P-values
    for (i in 1:nrow(pairs)){
        tg <- t.genes[[i]]
        spears <- tg.corr[[i]]

        # estimate the LR correlation P-value
        p.lr <- 1 - cdf(pairs$LR.corr[i], LR.par)

        # estimate the target gene correlation P-value based on rank statistics
        # for the individual correlation Gaussian model
        len <- pairs$len[i]
        r <- min(max(1, trunc(rank.p*len)), len)
        rank.corr <- spears[r]
        p.rt <- stats::pbinom(r-1, len, cdf(rank.corr, RT.par))
        pairs$pval[i] <- p.lr*p.rt
    }

    # recompute the Q-values
    rawp <- pairs$pval
    adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
    pairs$qval <- adj$adjp[order(adj$index),fdr.proc]

    # update the BSRInference object
    inf.param <- infParam(obj)
    inf.param$fdr.proc <- fdr.proc
    inf.param$rank.p <- rank.p
    infParam(obj) <- inf.param
    LRinter(obj) <- pairs

    obj

}) # rescoreInference


# Reduction and pathway stat methods ===========================================


if (!isGeneric("getPathwayStats")) {
    if (is.function("getPathwayStats"))
        fun <- getPathwayStats
    else
        fun <- function(obj, ...) standardGeneric("getPathwayStats")
    setGeneric("getPathwayStats", fun)
}
#' Basic statistics about hit pathways
#'
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
#' print('getPathwayStats')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' bsrinf <- initialInference(bsrdm)
#'
#' pw.stat <- getPathwayStats(bsrinf)
#' head(pw.stat)
#' 
#' @importFrom foreach %do% %dopar%
setMethod("getPathwayStats", "BSRInference", function(obj,
                                pval.thres=NULL, qval.thres=NULL){

    if (infParam(obj)$ligand.reduced || infParam(obj)$receptor.reduced)
        stop(paste0("Cannot be applied to interactions involving",
                    " reduced receptors or ligands"))

    # local binding
    id <- NULL

    pairs <- LRinter(obj)
    if (!is.null(pval.thres))
        pairs <- pairs[pairs$pval <= pval.thres,]
    else
        pairs <- pairs[pairs$qval <= qval.thres,]

    # pairs reduced to the receptor & pathway names
    pairs.R <- unique(pairs[, c("R","pw.id","pw.name")])
    upw <- unique(pairs[, c("pw.id","pw.name")])
    pw.names <- stats::setNames(upw$pw.name, upw$pw.id)

    # number of "hits" on a pathway with or without the ligand combinatorics
    pw.ids <- table(pairs$pw.id)
    pw.ids.R <- table(pairs.R$pw.id)

    # number of ligands for each receptor
    R.n.comb <- table(SingleCellSignalR::LRdb$receptor)

    foreach::foreach(id=names(pw.ids), .combine=rbind) %do% {

        # number of receptors that are in the current pathway,
        # depending on whether it is a GOBP or Reactome pathway
        if (regexpr("^R-",id) != -1)
            Rs <- intersect(SingleCellSignalR::LRdb$receptor,
                            reactome[reactome$`Reactome ID`==id,"Gene name"])
        else
            Rs <- intersect(SingleCellSignalR::LRdb$receptor,
                            gobp[gobp$`GO ID`==id,"Gene name"])

        # non-combinatorial version (ignore ligands)
        tot.R <- length(Rs)

        # ligand combinatorics included
        tot.LR <- sum(R.n.comb[Rs])

        data.frame(pw.id=id, pw.name=pw.names[id], n.R=pw.ids.R[id],
                   tot.R=tot.R, n.LR=pw.ids[id], tot.LR=tot.LR,
                   stringsAsFactors=FALSE)
    }

}) # getPathwayStats


if (!isGeneric("reduceToBestPathway")) {
    if (is.function("reduceToBestPathway"))
        fun <- reduceToBestPathway
    else
        fun <- function(obj, ...) standardGeneric("reduceToBestPathway")
    setGeneric("reduceToBestPathway", fun)
}
#' Keep one pathway per ligand-receptor pair
#'
#' @return A BSRInference object reduced to only report one pathway per
#' ligand-receptor pair. The pathway with the
#' smallest P-value is selected.
#'
#' @details During the execution of \code{pValuesLR}, ligand-receptor pairs
#' are evaluated in relation with pathways that allow checking receptor
#' downstream correlations. It is thus possible
#' that several pathways are reported for a same LR pair.
#'
#' @export
#' @examples
#' print('reduceToBestPathway')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redBP  <- reduceToBestPathway(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToBestPathway", "BSRInference", function(obj) {

    # Here we access the object slots directly as this procedure
    # is dependent of actual data representation

    # get best p-value pathway per LR pair
    ligands <- list()
    receptors <- list()
    t.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    LR <- unique(pairs[, c("L","R")])
    for (i in seq_len(nrow(LR))){
    
        L <- LR$L[i]
        R <- LR$R[i]
        
        pwr <- pairs[pairs$L==L & pairs$R==R,]
        k <- which.min(pwr$pval)
        j <- which(pairs$L==L & pairs$R==R & pairs$pw.id==pwr$pw.id[k])
        ligands <- c(ligands, obj@ligands[j])
        receptors <- c(receptors, obj@receptors[j])
        t.genes <- c(t.genes, obj@t.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        LRinter <- rbind(LRinter, pairs[j,])
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@t.genes <- t.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE


    obj

}) # reduceToBestPathway


if (!isGeneric("reduceToReceptor")) {
    if (is.function("reduceToReceptor"))
        fun <- reduceToReceptor
    else
        fun <- function(obj, ...) standardGeneric("reduceToReceptor")
    setGeneric("reduceToReceptor", fun)
}
#' Aggregate the ligands of a same receptor
#'
#' Simplifies a ligand-receptor table to focus on the receptors.
#'
#' @return BSRInference object reduced to one row per receptor.
#' All the ligands are combined in a
#' semi-colon-separated list surrounded by curly braces in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the line with the
#' pathway featuring the smallest P-value.
#' @export
#' @examples
#' print('reduceToReceptor')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redR  <- reduceToReceptor(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToReceptor", "BSRInference", function(obj){

    # Here we access the object slots directly as this procedure
    # is dependent of actual data representation

    if (infParam(obj)$ligand.reduced)
        stop("Already reduced to receptor") # because ligands were reduced

    # pool the ligands
    ligands <- list()
    receptors <- list()
    t.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (R in unique(pairs$R)){
        lig <- pairs[pairs$R==R,]
        k <- which.min(lig$pval)
        j <- which(pairs$R==lig$R[k] & pairs$pw.id==lig$pw.id[k])[1]
        ligands <- c(ligands, list(unique(lig$L)))
        receptors <- c(receptors, list(R))
        t.genes <- c(t.genes, obj@t.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j,]
        to.add[1, "L"] <- paste0("{", paste(unique(lig$L), collapse=";"), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@t.genes <- t.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$ligand.reduced <- TRUE

    obj

}) # reduceToReceptor


if (!isGeneric("reduceToLigand")) {
    if (is.function("reduceToLigand"))
        fun <- reduceToLigand
    else
        fun <- function(obj, ...) standardGeneric("reduceToLigand")
    setGeneric("reduceToLigand", fun)
}
#' Aggregate the receptors of a same ligand
#'
#' Simplifies a ligand-receptor table to focus on the ligands.
#'
#' @return A BSRInference object but reduced to one row per ligand.
#' All the receptors are combined in a
#' semi-colon-separated list surrounded by curly braces in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the pathway with
#' the smallest P-value.
#' @export
#' @examples
#' print('reduceToLigand')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redL  <- reduceToLigand(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToLigand", "BSRInference", function(obj){

    # Here we access the object slots directly as this procedure
    # is dependent of actual representation

    if (infParam(obj)$receptor.reduced)
        stop("Already reduced to ligand") # because receptors were reduced

    # pool the receptors
    ligands <- list()
    receptors <- list()
    t.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (L in unique(pairs$L)){
        rec <- pairs[pairs$L==L,]
        k <- which.min(rec$pval)
        j <- which(pairs$L==L & pairs$R==rec$R[k] & pairs$pw.id==rec$pw.id[k])
        ligands <- c(ligands, list(L))
        receptors <- c(receptors, list(unique(rec$R)))
        t.genes <- c(t.genes, obj@t.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j,]
        to.add[1, "R"] <- paste0("{", paste(unique(rec$R), collapse=";"), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@t.genes <- t.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$pathway.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj

}) # reduceToLigand


if (!isGeneric("reduceToPathway")) {
    if (is.function("reduceToPathway"))
        fun <- reduceToPathway
    else
        fun <- function(obj, ...) standardGeneric("reduceToPathway")
    setGeneric("reduceToPathway", fun)
}
#' Aggregate ligands and receptors at the pathway level
#'
#' Simplifies a ligand-receptor inference object to focus on
#' the pathways.
#'
#' @return A BSRInference object reduced to only report one row per pathway.
#' The information of which ligand interacted with which receptor is lost as
#' all the ligands and all the receptors forming pairs related to a certain
#' pathway are combined.
#' For a given pathway, the reported P-values and target genes are those of
#' the best ligand-receptor pair that
#' was in this pathway.
#' Receptors and ligands are combined in two semi-colon-separated
#' lists surrounded by curly braces in the tabular slot \code{LRinter},
#' while the list representation slots (\code{ligands} and
#' \code{receptors}) are update accordingly.
#'
#' @export
#' @examples
#' print('reduceToPathway')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' @importFrom rlang .data
setMethod("reduceToPathway", "BSRInference", function(obj){

    # Here we access the object slots directly as this procedure
    # is dependent of actual representation

    if (infParam(obj)$receptor.reduced)
        stop("Already reduced to ligand") # because receptors were reduced
    if (infParam(obj)$ligand.reduced)
        stop("Already reduced to receptor") # because ligands were reduced

    # reduce to unique pathways
    ligands <- list()
    receptors <- list()
    t.genes <- list()
    tg.corr <- list()
    LRinter <- NULL
    pairs <- obj@LRinter
    for (p in unique(pairs$pw.id)){
        j <- which(pairs$pw.id==p)[1]
        ligands <- c(ligands, list(unique(pairs$L[pairs$pw.id==p])))
        receptors <- c(receptors, list(unique(pairs$R[pairs$pw.id==p])))
        t.genes <- c(t.genes, obj@t.genes[j])
        tg.corr <- c(tg.corr, obj@tg.corr[j])
        to.add <- pairs[j,]
        to.add[1, "L"] <- paste0("{", paste(unique(pairs$L[pairs$pw.id==p]),
                                            collapse=";"), "}")
        to.add[1, "R"] <- paste0("{", paste(unique(pairs$R[pairs$pw.id==p]),
                                            collapse=";"), "}")
        LRinter <- rbind(LRinter, to.add)
    }

    # update the object
    obj@LRinter <- LRinter
    obj@ligands <- ligands
    obj@receptors <- receptors
    obj@t.genes <- t.genes
    obj@tg.corr <- tg.corr
    obj@inf.param$ligand.reduced <- TRUE
    obj@inf.param$receptor.reduced <- TRUE

    obj

}) # reduceToPathway


# Obtain gene signatures from a BSRInference object ============================

if (!isGeneric("getLRGeneSignatures")) {
    if (is.function("getLRGeneSignatures"))
        fun <- getLRGeneSignatures
    else
        fun <- function(obj, ...) standardGeneric("getLRGeneSignatures")
    setGeneric("getLRGeneSignatures", fun)
}
#' Extract gene signatures of LR pair activity
#'
#' Obtains gene signatures reflecting ligand-receptor as well as
#' receptor downstream activity to
#' score ligand-receptor pairs across samples subsequently with
#' \code{\link{scoreLRGeneSignatures}}.
#'
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @return A BSRSignature object containing a gene signature for each triple
#' ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is
#' comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}.
#' In case \code{signed==TRUE},
#' the rank is defined for correlation absolute values.
#' @export
#' @examples
#' print('getLRGeneSignatures')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' bsrsig.redP <- getLRGeneSignatures(bsrinf.redP,qval.thres=0.001)
#'
#' @importFrom foreach %do% %dopar%
setMethod("getLRGeneSignatures", "BSRInference", function(obj,
        pval.thres=NULL, qval.thres=NULL){

    if (is.null(pval.thres) && is.null(qval.thres))
        stop("Either a P- or a Q-value threshold must be provided")

    # reduce and select
    obj <- reduceToBestPathway(obj)
    pairs <- LRinter(obj)
    if (!is.null(pval.thres))
        selected <- pairs$pval <= pval.thres
    else
        selected <- pairs$qval <= qval.thres

    # obtain the signature object
    pairs     <- pairs[selected,]
    ligands   <- ligands(obj)[selected]
    receptors <- receptors(obj)[selected]
    pathways  <- pairs$pw.name
    #pathways  <- paste(pairs$pw.id, pairs$pw.name)
    t.genes   <- tGenes(obj)[selected]
    t.corrs   <- tgCorr(obj)[selected]
    
    for (i in seq_len(nrow(pairs))){
        tg <- t.genes[[i]]
            t.genes[[i]] <- tg[pairs$rank[i]:length(tg)] 

        tc <- t.corrs[[i]]
            t.corrs[[i]] <- tc[pairs$rank[i]:length(tc)] 
    }

    new("BSRSignature", pathways=pathways, ligands=ligands,
        receptors=receptors, t.genes=t.genes,tg.corr=t.corrs)

}) # getLRGeneSignatures


# Reset gene names to initial organism providen in first instance
# ==========================

if (!isGeneric("resetToInitialOrganism")) {
  if (is.function("resetToInitialOrganism"))
    fun <- resetToInitialOrganism
  else
    fun <- function(obj, ...) standardGeneric("resetToInitialOrganism")
  setGeneric("resetToInitialOrganism", fun)
}

#' # Reset gene names to initial organism providen in first instance
#'
#' @param conversion.dict   A dictionnary
#'
#' @return An BSRInference object updated for gene names.
#' The gene names are replaced by the ones from
#' the organism providen in first instance.
#'
#' @export
#' @examples
#'data(bodyMap.mouse)
#'ortholog.dict    <- findOrthoGenes (from_organism = "mmusculus", 
#'                                     from_values = rownames(bodyMap.mouse))
#' 
#'matrix.expression.human <- convertToHuman(counts = bodyMap.mouse, 
#'                                           dictionary = ortholog.dict)
#'
#'bsrdm  <- prepareDataset(counts = matrix.expression.human,
#'           species = "mmusculus",
#'           conversion.dict = ortholog.dict)
#'
#'bsrdm <- learnParameters(bsrdm, null.model="normal" , quick = FALSE, 
#'           plot.folder="../man/figures/",filename = "bodyMap.mouse",
#'           verbose = TRUE)
#' 
#'bsrinf  <- initialInference(bsrdm)
#' 
#'bsrinf  <- resetToInitialOrganism(bsrinf, conversion.dict=ortholog.dict)
#'
setMethod("resetToInitialOrganism", "BSRInference", function(obj,
                  conversion.dict=data.frame(Gene.name="A",row.names = "B") ){

     print("resetToInitialOrganism")
     # Need to check conversion.dict format
    
     conversion.dict$human.gene.name  <- rownames(conversion.dict) 

     LRinter(obj)$L <- .geneNameConversion(LRinter(obj)$L,conversion.dict)
     LRinter(obj)$R <- .geneNameConversion(LRinter(obj)$R,conversion.dict)

     ligands(obj)   <- .geneNameConversion(ligands(obj),conversion.dict)
     receptors(obj) <- .geneNameConversion(receptors(obj),conversion.dict)
     tGenes(obj)    <- .geneNameConversion(tGenes(obj),conversion.dict)

    obj

})

#' @title Convert gene symbol to another organism 
#'
#' @description Convert gene symbol to another organism 
#' based on a dictionnary with human and ortholog species.
#'
#' @param genes genes you want to convert
#' @param conversion.dict Dataframe containing
#' gene names for source species and Homo Sapiens.
#'
#' @return Depend type of input genes 
#' LRinter return a vector of genes 
#' tGenes receptors ligands : return list of list of genes
#'
.geneNameConversion <- function(genes,conversion.dict=data.frame(Gene.name="A",row.names = "B")){

    #print(".geneNameConversion")
    if(typeof(genes) == "character"){
        genes.df <- data.frame(human.gene.name = genes)
     
        genes.df$id <- 1:nrow(genes.df) 

        genes.converted <- merge(genes.df,conversion.dict,by.x='human.gene.name',sort=FALSE,all=FALSE)
        genes.converted <- genes.converted[order(genes.converted$id), ]

        genes.converted$human.gene.name <- NULL
        genes.converted$id <- NULL

        as.vector(unlist(genes.converted))
    }
    else if (typeof(genes) == "list") {
        list <- list()

        for (i in seq_len(length(genes))){
            genes.df <- data.frame(human.gene.name = genes[[i]])
            genes.df$id <- 1:nrow(genes.df) 
            genes.converted <- merge(genes.df,conversion.dict,by.x='human.gene.name',sort=FALSE,all=FALSE)
            genes.converted <- genes.converted[order(genes.converted$id), ]

            genes.converted$human.gene.name <- NULL
            genes.converted$id <- NULL

            list[[i]] <-  as.vector(unlist(genes.converted))    
            rm(genes.df)
            rm(genes.converted)
        }

        list
    }
    else {
        stop("Something went wrong during gene conversion.", call. = FALSE)
    }

    
}
#.geneNameConversion