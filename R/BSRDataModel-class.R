library(methods)

#' BulkSignalR Data Model Object
#'
#' An S4 class to represent the expression data used for inferring
#' ligand-receptor interactions.
#'
#' @slot ncounts   Normalized read count matrix. Row names must be set to HUGO
#' official gene symbols.
#' @slot log.transformed  Logical indicating whether values in
#' \code{ncounts} were log2-transformed.
#' @slot normalization    Name of the normalization method.
#' @slot param            List containing the statistical model parameters.
#' @slot initial.organism         Organism for which the data were obtained.
#' @slot initial.orthologs        List of genes for which human
#' orthologs exist.
#' @export
#' @examples
#' new("BSRDataModel", ncounts=matrix(1.5, nrow=2, ncol=2,
#'       dimnames=list(c("A","B"), c("C","D"))),
#'       log.transformed=TRUE,
#'       normalization="TC")
#'
setClass("BSRDataModel",
         slots=c(initial.organism="character",
                 initial.orthologs="list",
                 ncounts="matrix",
                 log.transformed="logical",
                 normalization="character",
                 param="list"
                ),
         prototype=list(
             initial.organism="hsapiens",
             initial.orthologs=list("A","B","C"),
             ncounts=matrix(1.0,nrow=2,ncol=1,dimnames=list(c("A","B"),"C")),
             log.transformed=FALSE,
             normalization="UQ",
             param=list()

         ))

setValidity("BSRDataModel",
    function(object) {
        if(!is.matrix(object@ncounts))
            return("specified normalized counts are not a matrix")
        if(!is.numeric(object@ncounts))
            return("specified normalized counts are not numeric")
        if (is.null(row.names(object@ncounts)))
            return("specified normalized counts have no row names set")
        if (nchar(object@normalization) == 0)
            return("Normalization method must be set")

        TRUE
    }
)

setMethod("show", "BSRDataModel",
    function(object) {
        cat("Expression values are log2-transformed: ", object@log.transformed,
            "\n", sep="")
        cat("Normalization method: ", object@normalization, "\n", sep="")
        cat("Organism : ", object@initial.organism,"\n", sep="")
        cat("Statistical model parameters:\n")
        utils::str(object@param)
        cat("Expression data:\n")
        if (ncol(object@ncounts) > 8)
            print(head(object@ncounts[,1:8]))
        else
            print(head(object@ncounts))
    }
)


# Accessors & setters ========================================================

if (!isGeneric("initialOrganism")) {
    if (is.function("initialOrganism"))
        fun <- initialOrganism
    else
        fun <- function(x) standardGeneric("initialOrganism")
    setGeneric("initialOrganism", fun)
}
#' organism accessor
#'
#' @name initialOrganism
#' @aliases initialOrganism,BSRDataModel-method
#' @param x Object BSRDataModel
#' @export
setMethod("initialOrganism", "BSRDataModel", function(x) x@initial.organism)


if (!isGeneric("initialOrthologs")) {
    if (is.function("initialOrthologs"))
        fun <- initialOrthologs
    else
        fun <- function(x) standardGeneric("initialOrthologs")
    setGeneric("initialOrthologs", fun)
}
#' Model parameter accessor
#'
#' @name initialOrthologs
#' @aliases initialOrthologs,BSRDataModel-method
#' @param x Object BSRDataModel
#' @export
setMethod("initialOrthologs", "BSRDataModel", function(x) x@initial.orthologs)


if (!isGeneric("ncounts")) {
    if (is.function("ncounts"))
        fun <- ncounts
    else
        fun <- function(x) standardGeneric("ncounts")
    setGeneric("ncounts", fun)
}
#' Normalized count matrix accessor
#'
#' @name ncounts
#' @aliases ncounts,BSRDataModel-method
#' @param x object BSRDataModel 
#' @export
setMethod("ncounts", "BSRDataModel", function(x) x@ncounts)

if (!isGeneric("ncounts<-")) {
    if (is.function("ncounts<-"))
        fun <- `ncounts<-`
    else
        fun <- function(x, value) standardGeneric("ncounts<-")
    setGeneric("ncounts<-", fun)
}
#' Normalized count matrix setter (internal use only)
#'
#' @param x object BSRDataModel 
#' @param value valut to be set for BSRDataModel 
setMethod("ncounts<-", "BSRDataModel", function(x,value){
    x@ncounts <- value
    methods::validObject(x)
    x
})


if (!isGeneric("param")) {
    if (is.function("param"))
        fun <- param
    else
        fun <- function(x) standardGeneric("param")
    setGeneric("param", fun)
}
#' Model parameter accessor
#'
#' @name param
#' @aliases param,BSRDataModel-method
#' @param x BSRDataModel oject
#' @export
setMethod("param", "BSRDataModel", function(x) x@param)

if (!isGeneric("logTransformed")) {
    if (is.function("logTransformed"))
        fun <- logTransformed
    else
        fun <- function(x) standardGeneric("logTransformed")
    setGeneric("logTransformed", fun)
}
#' log.transformed accessor
#'
#' @name logTransformed
#' @aliases logTransformed,BSRDataModel-method
#' @param x Object BRSDataModel
#' @export
setMethod("logTransformed", "BSRDataModel", function(x) x@log.transformed)

if (!isGeneric("normalization")) {
    if (is.function("normalization"))
        fun <- normalization
    else
        fun <- function(x) standardGeneric("normalization")
    setGeneric("normalization", fun)
}
#' Normalization accessor
#'
#' @name normalization
#' @aliases normalization,BSRDataModel-method
#' @param x oject BSRDatamModel 
#' @export
setMethod("normalization", "BSRDataModel", function(x) x@normalization)



if (!isGeneric("learnParameters")) {
    if (is.function("learnParameters"))
        fun <- learnParameters
    else
        fun <- function(obj, ...) standardGeneric("learnParameters")
    setGeneric("learnParameters", fun)
}
#' Training of BulkSignalR model parameters
#'
#' Unique entry point for training the parameters behind
#' BulkSignalR statistical models.
#'
#' @name learnParameters
#' @aliases learnParameters,BSRDataModel-method
#'
#' @param obj   A BSRDatamodel without learned paramaters.
#' @param plot.folder   A folder name for generating control plots.
#' @param filename      Name of the output plot.
#' @param verbose       A logical activating progress messages for the user.
#' @param n.rand.LR     The number of random expression matrices to use for
#'   learning the ligand-receptor correlation distribution.
#' @param n.rand.RT     The number of random expression matrices to use for
#'   learning the receptor-target genes correlation distribution.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param quick           A logical indicating whether approximate parameters
#'   for the receptor-target correlations should be used.
#' @param null.model      The null model to use for Spearman correlation
#'   null distributions.
#'
#' @details Estimates the model parameters that are stored in the
#'   slot \code{param}.
#'
#'   In a reference pathway, i.e., a Reactome pathway or the genes of a
#'   GOBP term, the target genes are the genes coding for proteins forming a
#'   complex with the receptor and the genes in the pathway downstream the
#'   receptor, which are given as regulated by the pathway. If
#'   \code{with.complex} is set to \code{FALSE}, then only the regulated genes
#'   are considered. Participation to a complex, being regulated, and
#'   pathway directed topologies are defined by Reactome and KEGG pathways
#'   as provided by PathwayCommons.
#'
#'   The \code{min.pw.size}, \code{max.pw.size}, and \code{min.positive}
#'   parameters should be identical to the
#'   values intended when searching for ligand-receptor pairs with
#'   \code{\link{.getCorrelatedLR}}) and  \code{\link{.checkReceptorSignaling}})
#'   Although the statistical distributions are rather robust, it is not
#'   advisable to use different parameters that could introduce unanticipated
#'   biases, but for saving compute time and exploring.
#'
#'   The maximum pathway size is used to limit the redundancy inherent to GOBP
#'   and Reactome. The minimum pathway size is used to avoid overspecific,
#'   noninformative results.
#'
#'   BulkSignalR approach relies on modeling (Spearman) correlations and
#'   different models of null distributions are available for this purpose
#'   (parameter \code{null.model}). By default, the "automatic" option is
#'   selected meaning that censored normal and mixed normal as well as
#'   an empirical model based on Gaussian kernels (R \code{density()} function)
#'   are compared to pick the one closest to the data. Preference is given
#'   to normal and then mixture of normal over the empirical version for
#'   comparable quality of fit. It is also to bypass the automatic selection.
#'   Fitting of an alpha-stable distribution is quite time consuming as the
#'   computation of its PDF is compute-intensive. Finally, in the automaic
#'   selection mode, the choice of the actual model will be done based on
#'   the L-R null assuming a similar shape for the R-T null (with
#'   different parameters though, unless \code{quick} was set to \code{TRUE}).
#'
#' @return A BSRDataModel object with trained model parameters
#'
#' @export
#'
#' @examples
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
#' print("learnParameters")
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
#' print("learnParameters")
#' bsrdm <- learnParameters(bsrdm,plot.folder="./")
#' bsrdm
#' @importFrom methods new
setMethod("learnParameters", "BSRDataModel", function(obj, plot.folder = NULL,
      verbose = FALSE, n.rand.LR = 5L, n.rand.RT = 2L, with.complex = TRUE,
      max.pw.size = 200, min.pw.size = 5, min.positive = 4, quick = FALSE,
      null.model = c("automatic", "mixedNormal", "normal", "kernelEmpirical",
                     "empirical", "stable"), filename = "distribution") {

    obj@param$n.rand.LR <- as.integer(n.rand.LR)
    if (obj@param$n.rand.LR < 1)
        stop("Parameter n.rand.LR must be an integer > 0")
    obj@param$n.rand.RT <- as.integer(n.rand.RT)
    if (obj@param$n.rand.RT < 1)
        stop("Parameter n.rand.RT must be an integer > 0")
    obj@param$plot.folder <- plot.folder
    if (!is.null(plot.folder) && !file.exists(plot.folder))
        stop("The provided plot.folder does not exist")
    obj@param$file.name <- filename
    obj@param$with.complex <- with.complex
    obj@param$max.pw.size <- trunc(max.pw.size)
    if (obj@param$max.pw.size < 1)
        stop("Parameter max.pw.size must be an integer > 0")
    obj@param$min.pw.size <- trunc(min.pw.size)
    if (obj@param$min.pw.size < 1 || obj@param$min.pw.size > obj@param$max.pw.size)
        stop("Parameter min.pw.size must be an integer > 0 and <= than max.pw.size")
    obj@param$min.positive <- trunc(min.positive)
    if (obj@param$min.positive < 1)
        stop("Parameter min.positive must be an integer > 0")
    obj@param$quick <- quick
    null.model <- match.arg(null.model)
    if (null.model == "normal")
        trainModel <- .getGaussianParam
    else if (null.model == "mixedNormal")
        trainModel <- .getMixedGaussianParam
    else if (null.model == "kernelEmpirical")
        trainModel <- .getKernelEmpiricalParam
    else if (null.model == "empirical")
        trainModel <- .getEmpiricalParam
    else if (null.model == "stable")
        trainModel <- .getAlphaStableParam
    else if (null.model == "automatic")
        trainModel <- NULL
    else
        stop("No valid null model specified")

    # LR correlation null ----------------

    if (verbose)
        cat("Learning ligand-receptor correlation null distribution...\n")
    obj@param$min.corr.LR <- -1.0
    ds.LR.null <- .getEmpiricalNullCorrLR(obj@ncounts,
                                          n.rand=obj@param$n.rand.LR,
                                          min.cor=obj@param$min.corr.LR)
    rc <- ds.LR.null[[1]]$corr
    if (length(ds.LR.null) > 1)
        for (i in 2:length(ds.LR.null)) rc <- c(rc, ds.LR.null[[i]]$corr)
    obj@param$LR.0$n <- length(rc)

    # Null distribution model
    if (!is.null(plot.folder))
        file.name <- paste0(plot.folder, "/",filename,"_LR-null.pdf")
    else
        file.name <- NULL

    if (is.null(trainModel)){
        # automatic selection of the model
        np <- try(.getGaussianParam(rc, "LR correlation (null)"), silent=TRUE)
        if (inherits(np, "try-error"))
            np <- NULL
        mp <- try(.getMixedGaussianParam(rc, "LR correlation (null)"),
                  silent=TRUE)
        if (inherits(mp, "try-error"))
            mp <- NULL
        kp <- .getKernelEmpiricalParam(rc, "LR correlation (null)")
        if (verbose){
            cat("Automatic null model choice:\n")
            if (is.null(np))
                cat("  Censored normal estimation did not converge\n")
            else
                cat("  Censored normal D=", np$D,
                    ", Chi2=", np$Chi2, "\n", sep="")
            if (is.null(mp))
                cat("  Censored Mixture of normals estimation did not converge\n")
            else
                cat("  Censored mixture D=", mp$D,
                    ", Chi2=", mp$Chi2, "\n", sep="")
            cat("  Gaussian kernel empirical D=", kp$D,
                ", Chi2=", kp$Chi2, "\n", sep="")
        }
        npchi <- ifelse(is.null(np), 100, sqrt(np$Chi2))
        mpchi <- ifelse(is.null(mp), 100, sqrt(mp$Chi2))
        kpchi <- sqrt(kp$Chi2)
        if ((npchi < 1.25*mpchi) && (npchi < 2*kpchi)){
            trainModel <- .getGaussianParam
            if (verbose)
                cat("  ==> select censored normal\n")
        }
        else if (mpchi < 2*kpchi){
            trainModel <- .getMixedGaussianParam
            if (verbose)
                cat("  ==> select censored mixture of 2 normals\n")
        }
        else{
            trainModel <- .getKernelEmpiricalParam
            if (verbose)
                cat("  ==> select Gaussian kernel-based empirical\n")
        }
    }
    # actual training with the chosen model
    gp <- trainModel(rc, "LR correlation (null)", verbose = verbose,
                     file.name = file.name)
    obj@param$LR.0$model <- gp

    # RT correlation null ------------------------------------

    if (obj@param$quick) {
        # RT correlations are assumed to be equal to LR correlations
        if (verbose)
            cat(paste0("Quick learning, receptor-target correlation null ","
                   distribution assumed to be equal to ligand-receptor...\n"))
        obj@param$RT.0$n <- obj@param$LR.0$n
        obj@param$RT.0$model <- obj@param$LR.0$model
    }
    else {
        # RT correlations are actually learnt
        if (verbose)
            cat("Learning receptor-target correlation null distribution...\n")
        ds.RT.null <- .getEmpiricalNull(obj@ncounts,
                                        n.rand = obj@param$n.rand.RT,
                                        with.complex = obj@param$with.complex,
                                        max.pw.size = obj@param$max.pw.size,
                                        min.pw.size = obj@param$min.pw.size,
                                        min.positive = obj@param$min.positive,
                                        min.cor = obj@param$min.corr.LR)
        t <- ds.RT.null[[1]]
        if (length(ds.RT.null) > 1)
            for (i in 2:length(ds.RT.null)) t <- rbind(t, ds.RT.null[[i]])
        above <- unlist(strsplit(t$target.corr, split = "\\|"))
        r.corrs <- NULL
        for (i in seq_len(length(above))) {
            corr <- as.numeric(strsplit(above[i], split = ";")[[1]])
            r.corrs <- c(r.corrs, corr)
        }
        if (null.model == "stable")
            # sub-sample randomized R-T correlations to limit compute time
            r.corrs <- sample(r.corrs, obj@param$LR.0$n)
        obj@param$RT.0$n <- length(r.corrs)

        # fit null model
        if (!is.null(plot.folder))
            file.name <- paste0(plot.folder,  "/",filename,"_RT-null.pdf")
        else
            file.name <- NULL
        gp <- trainModel(r.corrs, "RT correlation (null)",
                                verbose = verbose, file.name = file.name)
        obj@param$RT.0$model <- gp
    }

    if (verbose)
        cat("Learning of statistical model parameters completed\n")
    obj

})  # learnParameters


# performing initial inference ===========================================


if (!isGeneric("initialInference")) {
    if (is.function("initialInference"))
        fun <- initialInference
    else
        fun <- function(obj, ...) standardGeneric("initialInference")
    setGeneric("initialInference", fun)
}
#' Inference of ligand-receptor interactions
#'
#' Computes putative LR interactions along with their statistical confidence.
#' In this initial inference, all the relevant pathways are reported,
#' see reduction functions to reduce this list.
#'
#' @name initialInference
#' @aliases initialInference,BSRDataModel-method
#'
#' @param obj         A BSRDataModel output by \code{\link{prepareDataset}} with
#' statistical model parameters trained by 
#' \code{"\link[=BSRDataModel-class]{learnParameters}"}
#'
#' method.
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
#'
#' @details Perform the initial ligand-receptor inference. Initial means that
#' reduction is applied, all the (ligand, receptor, downstream pathway)
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
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
#' print("learnParameters")
#' bsrdm <- learnParameters(bsrdm)
#' bsrdm
#'
#' print('perform inference')
#' bsrinf <- initialInference(bsrdm)
#' bsrinf
#' @importFrom methods new
setMethod("initialInference", "BSRDataModel", function(obj, rank.p=0.55,
        min.cor = 0.25,
        restrict.genes = NULL, reference=c("REACTOME-GOBP","REACTOME","GOBP"),
        max.pw.size=NULL, min.pw.size=NULL, min.positive=NULL, restrict.pw=NULL,
        with.complex=NULL, fdr.proc=c("BH","Bonferroni","Holm","Hochberg",
                                      "SidakSS","SidakSD","BY","ABH","TSBH")){

    if (is.null(max.pw.size))
        max.pw.size <- param(obj)$max.pw.size
    if (is.null(min.pw.size))
        min.pw.size <- param(obj)$min.pw.size
    if (is.null(min.positive))
        min.positive <- param(obj)$min.positive
    if (is.null(with.complex))
        with.complex <- param(obj)$with.complex
    reference <- match.arg(reference)
    fdr.proc <- match.arg(fdr.proc)
    if (rank.p < 0 || rank.p > 1)
        stop("rank.p must lie in [0;1]")

    inf.param <- list()
    inf.param$min.corr <- min.cor
    inf.param$restrict.genes <- restrict.genes
    lr <- .getCorrelatedLR(obj, min.cor=min.cor, restrict.genes=restrict.genes)

    inf.param$reference <- reference
    inf.param$min.pw.size <- min.pw.size
    inf.param$max.pw.size <- max.pw.size
    inf.param$with.complex <- with.complex
    inf.param$min.positive <- min.positive
    inf.param$restrict.pw <- restrict.pw
    pairs <- .checkReceptorSignaling(obj, lr, reference=reference,
                         min.pw.size=min.pw.size, max.pw.size=max.pw.size,
                         min.positive=min.positive, with.complex=with.complex,
                         restrict.pw=restrict.pw)

    inf.param$fdr.proc <- fdr.proc
    inf.param$rank.p <- rank.p
    inter <- .pValuesLR(pairs, param(obj), rank.p=rank.p, fdr.proc=fdr.proc)

    ligands <- strsplit(inter$L, ";")
    receptors <- strsplit(inter$R, ";")
    tg <- strsplit(inter$target.genes, ";")
    tgcorr <- lapply(strsplit(inter$target.corr, ";"),
                     function(x) as.numeric(x))
    inf.param$ligand.reduced <- FALSE
    inf.param$receptor.reduced <- FALSE
    inf.param$pathway.reduced <- FALSE

    new("BSRInference", LRinter=inter[,c("L","R","LR.corr","pw.id","pw.name",
        "rank","len","rank.corr","pval","qval")], ligands=ligands,
        receptors=receptors, t.genes=tg, tg.corr=tgcorr,
        inf.param=inf.param)

}) # initialInference


# Scoring of gene signatures in a BSRSignature object ==========================

if (!isGeneric("scoreLRGeneSignatures")) {
  if (is.function("scoreLRGeneSignatures"))
    fun <- scoreLRGeneSignatures
  else
    fun <- function(obj, ...) standardGeneric("scoreLRGeneSignatures")
  setGeneric("scoreLRGeneSignatures", fun)
}
#' Score ligand-receptor gene signatures
#'
#' Compute ligand-receptor gene signature scores over a BSRDataModel.
#'
#' @name scoreLRGeneSignatures
#' @aliases scoreLRGeneSignatures,BSRDataModel-method
#'
#' @param obj           A BSRDataModel object.
#' @param sig           A BSRSignature object.
#' @param LR.weight    A number between 0 and 1 defining the relative weight
#' of the ligand and the receptor in the signature.
#' @param robust       A logical indicating that z-scores should be computed
#' with median and MAD instead of mean and standard deviation.
#' @param name.by.pathway     A logical indicating whether row names of the
#' resulting score matrix should be pathway names.
#' @param rownames.LRP A logical indicating, in case \code{name.by.pathway}
#' was set to TRUE, whether ligand and receptor names should be added on top.
#' No role if \code{name.by.pathway} was set to FALSE.
#' @param abs.z.score  A logical to use absolute z-scores (useful if the
#' activity of a paythway is reported by a mixture of up- and down-genes
#' whose z-score averages might hide actual activity).
#' @return A matrix containing the scores of each ligand-receptor gene
#' signature in each sample.
#'
#' @export
#' @examples
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
#' print("learnParameters")
#' bsrdm <- learnParameters(bsrdm)
#' bsrdm
#'
#' print('perform inference')
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redBP <- reduceToBestPathway(bsrinf) 
#' bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP,
#' qval.thres=0.001)
#'
#' scoresLR <- scoreLRGeneSignatures(bsrdm,bsrsig.redBP,
#'                        name.by.pathway=FALSE)
#' @importFrom foreach %do% %dopar%
#' @importFrom methods is
setMethod("scoreLRGeneSignatures", "BSRDataModel", function(obj,
                  sig, LR.weight=0.5, robust=FALSE,
                  name.by.pathway=FALSE, abs.z.score=FALSE,rownames.LRP=FALSE){

    if (!is(sig, "BSRSignature"))
        stop("sig must be a BSRSignature object")
    if (LR.weight<=0 || LR.weight>=1)
        stop("LRweight must reside in (0;1)")


    if( initialOrganism(obj)!="hsapiens" )
        all.genes <- unlist(initialOrthologs(obj))
    else
        all.genes <- rownames(ncounts(obj))

    # intersect signature gene names with RNA-seq data
    ncounts <- ncounts(obj)
    ligands <- sapply(ligands(sig), function(x) intersect(x, all.genes))
    receptors <- sapply(receptors(sig), function(x) intersect(x, all.genes))
    t.genes <- sapply(tGenes(sig), function(x) intersect(x, all.genes))

    good <- sapply(ligands, length) > 0 & sapply(receptors, length) > 0 &
            sapply(t.genes, length) > 0
    ligands   <- ligands[good]
    receptors <- receptors[good]
    t.genes   <- t.genes[good]
    pathways  <- pathways(sig)[good]

    # scale ncounts
    if (logTransformed(obj))
        ncounts <- 2**ncounts
    if (robust)
        z <- (ncounts-apply(ncounts,1,stats::median))/apply(ncounts,1,stats::mad)
    else
        z <- (ncounts-rowMeans(ncounts))/apply(ncounts,1,stats::sd)
    if (abs.z.score)
        z <- abs(z)

    if( initialOrganism(obj) != "hsapiens" )
        rownames(z) <- all.genes

    # compute the LR gene signatures
    i <- NULL
    pwn <- foreach::foreach(i=seq_len(length(pathways)), .combine=c) %do% {

        if (name.by.pathway){

            if(rownames.LRP){

                paste0("{",paste(ligands[[i]], collapse=";") ,"} / {",
                    paste(receptors[[i]], collapse=";"),"} | ",pathways[[i]])
            }
            else
                pathways[[i]]

        }
        else if (!name.by.pathway){
            paste0("{",paste(ligands[[i]], collapse=";") ,"} / {",
                    paste(receptors[[i]], collapse=";"),"}")
        }

    }

    res <- matrix(0,nrow=length(pathways),ncol=ncol(ncounts),dimnames=list(pwn,colnames(ncounts)))
    for (i in seq_len(length(pathways))){

        # average ligand z-score
        zz <- z[ligands[[i]],]
        if (is.matrix(zz))
            mL <- colSums(zz)/length(ligands[[i]])
        else
            mL <- zz

        # average receptor z-score
        zz <- z[receptors[[i]],]
        if (is.matrix(zz))
            mR <- colSums(zz)/length(receptors[[i]])
        else
            mR <- zz

        # average target gene z-score
        zz <- z[t.genes[[i]],]
        if (is.matrix(zz))
          mT <- colSums(zz)/length(t.genes[[i]])
        else
          mT <- zz

        res[i,] <- LR.weight*0.5*(mL+mR)+(1-LR.weight)*mT
    }

    res

}) # scoreLRGeneSignatures
