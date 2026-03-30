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
#' BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#'
setClass("BSRDataModel",
    slots = c(
        initial.organism = "character",
        initial.orthologs = "list",
        ncounts = "matrix",
        log.transformed = "logical",
        normalization = "character",
        param = "list"
    ),
    prototype = list(
        initial.organism = "hsapiens",
        initial.orthologs = list(), #list("A", "B", "C"),
        ncounts = matrix(1.0,
            nrow = 2, ncol = 1,
            dimnames = list(c("A", "B"), "C")
        ),
        log.transformed = FALSE,
        normalization = "UQ",
        param = list(spatial.smooth = FALSE)
    )
)

setValidity(
    "BSRDataModel",
    function(object) {
        if (!is.matrix(object@ncounts)) {
            return("specified normalized counts are not a matrix")
        }
        if (!is.numeric(object@ncounts)) {
            return("specified normalized counts are not numeric")
        }
        if (is.null(row.names(object@ncounts))) {
            return("specified normalized counts have no row names set")
        }
        if (nchar(object@normalization) == 0) {
            return("Normalization method must be set")
        }

        TRUE
    }
)

# NOTE: the constructor of BRSDataModel class is defined in the
# file dataPrepare.R for "historical" reasons

setMethod(
    "show", "BSRDataModel",
    function(object) {
        cat("Expression values are log2-transformed: ", object@log.transformed,
            "\n",
            sep = ""
        )
        cat("Normalization method: ", object@normalization, "\n", sep = "")
        cat("Organism: ", object@initial.organism, "\n", sep = "")
        cat("Statistical model parameters:\n")
        utils::str(object@param)
        cat("Expression data:\n")
        if (ncol(object@ncounts) > 8) {
            print(utils::head(object@ncounts[, seq_len(8)]))
        } else {
            print(utils::head(object@ncounts))
        }
    }
)


# Accessors & setters ========================================================

setGeneric("initialOrganism", signature="x",
    function(x) standardGeneric("initialOrganism")
)
#' organism accessor
#'
#' @name initialOrganism
#' @aliases initialOrganism,BSRDataModel-method
#' @param x Object BSRDataModel
#' @return initialOrganism
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' initialOrganism(bsrdm)
#' @export
setMethod("initialOrganism", "BSRDataModel", function(x) x@initial.organism)


setGeneric("initialOrthologs", signature="x",
    function(x) standardGeneric("initialOrthologs")
)
#' Model parameter accessor
#'
#' @name initialOrthologs
#' @aliases initialOrthologs,BSRDataModel-method
#' @param x Object BSRDataModel
#' @return initialOrthologs
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' initialOrthologs(bsrdm)
#' @export
setMethod("initialOrthologs", "BSRDataModel", function(x) x@initial.orthologs)


setGeneric("ncounts", signature="x",
    function(x) standardGeneric("ncounts")
)
#' Normalized count matrix accessor
#'
#' @name ncounts
#' @aliases ncounts,BSRDataModel-method
#' @param x object BSRDataModel
#' @return ncounts
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' ncounts(bsrdm)
#' @export
setMethod("ncounts", "BSRDataModel", function(x) x@ncounts)


setGeneric("ncounts<-", signature=c("x", "value"),
    function(x, value) standardGeneric("ncounts<-")
)
#' Normalized count matrix setter (internal use only)
#'
#' @param x object BSRDataModel
#' @param value value to be set for BSRDataModel
#' @return returns \code{NULL}
#' @keywords internal
setMethod("ncounts<-", "BSRDataModel", function(x, value) {
    x@ncounts <- value
    methods::validObject(x)
    x
})


setGeneric("parameters", signature="x",
    function(x) standardGeneric("parameters")
)
#' Model parameter accessor
#'
#' @name parameters
#' @aliases parameters,BSRDataModel-method
#' @param x BSRDataModel oject
#' @return param
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' parameters(bsrdm)
#' @export
setMethod("parameters", "BSRDataModel", function(x) x@param)


setGeneric("parameters<-", signature=c("x", "value"),
    function(x, value) standardGeneric("parameters<-")
)
#' Parameters dataModel setter (internal use only)
#'
#' @param x object BSRDataModel
#' @param value value to be set for BSRDataModel
#' @return returns \code{NULL}
#' @keywords internal
setMethod("parameters<-", "BSRDataModel", function(x, value) {
    x@param <- value
    methods::validObject(x)
    x
})


setGeneric("logTransformed", signature="x",
    function(x) standardGeneric("logTransformed")
)
#' log.transformed accessor
#'
#' @name logTransformed
#' @aliases logTransformed,BSRDataModel-method
#' @param x Object BRSDataModel
#' @return logTransformed
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' logTransformed(bsrdm)
#' @export
setMethod("logTransformed", "BSRDataModel", function(x) x@log.transformed)


setGeneric("normalization", signature="x",
    function(x) standardGeneric("normalization")
)
#' Normalization accessor
#'
#' @name normalization
#' @aliases normalization,BSRDataModel-method
#' @param x object BSRDatamModel
#' @return normalization
#' @examples
#' bsrdm <- BSRDataModel(
#'     counts = matrix(1.5,
#'                 nrow = 2, ncol = 2,
#'                 dimnames = list(c("A", "B"), c("C", "D"))
#'                    ),
#'     method = "TC",
#'     log.transformed = TRUE,
#'     normalize = FALSE,
#'     min.LR.found = 0
#' )
#' normalization(bsrdm)
#' @export
setMethod("normalization", "BSRDataModel", function(x) x@normalization)


setGeneric("learnParameters", signature="obj",
    function(obj, ...) standardGeneric("learnParameters")
)
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
#' @param min.corr.LR      The minimum ligand-receptor correlation required.
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
#' Note that since the introduction of the use.full.network parameter
#' (April 29, 2024) in the BSRInference method parameters,
#' the pathway sizes are always computed before potential
#' intersection with the observed data (use.full.network set to FALSE) for
#' consistency. Accordingly, the minimum and maximum pathway default values
#' have been raised from 5 & 200 to 5 & 400 respectively. By default,
#' use.full.network is set to TRUE, meaning no intersection and hence larger
#' pathways.
#'
#' @return A BSRDataModel object with trained model parameters
#'
#' @export
#'
#' @examples
#' data(sdc, package = "BulkSignalR")
#' idx <- sample(nrow(sdc), 4000)
#' bsrdm <- BSRDataModel(sdc[idx, c("N22","SDC17")],min.LR.found = 20)
#' bsrdm <- learnParameters(bsrdm, n.rand.LR = 1L,
#' verbose=FALSE,quick=TRUE)
#' @importFrom methods new
setMethod(
    "learnParameters", "BSRDataModel",
    function(obj, plot.folder = NULL,
        verbose = FALSE, n.rand.LR = 5L, n.rand.RT = 2L, with.complex = TRUE,
        max.pw.size = 400, min.pw.size = 5, min.positive = 4, quick = FALSE,
        null.model = c(
            "automatic", "mixedNormal", "normal", "kernelEmpirical",
            "empirical", "stable"
            ), 
        filename = NULL, min.corr.LR = -1) {


        parameters(obj)$n.rand.LR <- as.integer(n.rand.LR)
        if (parameters(obj)$n.rand.LR < 1) {
            stop("Parameter n.rand.LR must be an integer > 0")
        }
        parameters(obj)$n.rand.RT <- as.integer(n.rand.RT)
        if (parameters(obj)$n.rand.RT < 1) {
            stop("Parameter n.rand.RT must be an integer > 0")
        }
        parameters(obj)$plot.folder <- plot.folder
        if (!is.null(plot.folder) && !file.exists(plot.folder)) {
            stop("The provided plot.folder does not exist")
        }
        parameters(obj)$file.name <- filename
        parameters(obj)$with.complex <- with.complex
        parameters(obj)$max.pw.size <- trunc(max.pw.size)
        if (parameters(obj)$max.pw.size < 1) {
            stop("Parameter max.pw.size must be an integer > 0")
        }
        parameters(obj)$min.pw.size <- trunc(min.pw.size)
        if (parameters(obj)$min.pw.size < 1 || 
            parameters(obj)$min.pw.size > parameters(obj)$max.pw.size) {
            stop(
                "Parameter min.pw.size must be",
                "an integer > 0 and <= than max.pw.size"
            )
        }
        parameters(obj)$min.positive <- trunc(min.positive)
        if (parameters(obj)$min.positive < 1) {
            stop("Parameter min.positive must be an integer > 0")
        }

        if (min.corr.LR < -1 || min.corr.LR > 1) {
            stop("min.corr.LR must lie in [-1;1]")
        }

        parameters(obj)$quick <- quick
        null.model <- match.arg(null.model)
        if (null.model == "normal") {
            trainModel <- .getGaussianParam
        } else if (null.model == "mixedNormal") {
            trainModel <- .getMixedGaussianParam
        } else if (null.model == "kernelEmpirical") {
            trainModel <- .getKernelEmpiricalParam
        } else if (null.model == "empirical") {
            trainModel <- .getEmpiricalParam
        } else if (null.model == "stable") {
            trainModel <- .getAlphaStableParam
        } else if (null.model == "automatic") {
            trainModel <- NULL
        } else {
            stop("No valid null model specified")
        }

        # LR correlation null ----------------

        if (verbose) {
            message("Learning ligand-receptor correlation null distribution...")
        }

        parameters(obj)$min.corr.LR <- min.corr.LR
        
        ds.LR.null <- .getEmpiricalNullCorrLR(obj)

        rc <- ds.LR.null[[1]]$corr
        if (length(ds.LR.null) > 1) {
            for (i in 2:length(ds.LR.null)) rc <- c(rc, ds.LR.null[[i]]$corr)
        }
        parameters(obj)$LR.0$n <- length(rc)

        # Null distribution model
        if (!is.null(plot.folder)) {
            file.name <- paste0(plot.folder, "/", filename, "_LR-null.pdf")
        } else {
            file.name <- NULL
        }

        if (is.null(trainModel)) {
            # automatic selection of the model
            np <- try(.getGaussianParam(rc, 
                "LR correlation (null)"), silent = TRUE)
            if (inherits(np, "try-error")) {
                np <- NULL
            }
            mp <- try(.getMixedGaussianParam(rc, "LR correlation (null)"),
                silent = TRUE
            )
            if (inherits(mp, "try-error")) {
                mp <- NULL
            }
            kp <- .getKernelEmpiricalParam(rc, "LR correlation (null)")
            if (verbose) {
              message("Automatic null model choice:")
              if (is.null(np)) {
                message("  Censored normal estimation did not converge")
              } else {
                message(
                    "  Censored normal D=", np$D,
                    ", Chi2=", np$Chi2
                )
              }
              if (is.null(mp)) {
                message("  Censored Mixture of normals",
                    " estimation did not converge")
              } else {
                message(
                    "  Censored mixture D=", mp$D,
                    ", Chi2=", mp$Chi2
                )
              }
              message(
                  "  Gaussian kernel empirical D=", kp$D,
                  ", Chi2=", kp$Chi2
              )
            }
            npchi <- ifelse(is.null(np), 100, sqrt(np$Chi2))
            mpchi <- ifelse(is.null(mp), 100, sqrt(mp$Chi2))
            kpchi <- sqrt(kp$Chi2)
            if ((npchi < 1.25 * mpchi) && (npchi < 2 * kpchi)) {
              trainModel <- .getGaussianParam
              if (verbose) {
                  message("  ==> select censored normal")
              }
            } else if (mpchi < 2 * kpchi) {
              trainModel <- .getMixedGaussianParam
              if (verbose) {
                  message("  ==> select censored mixture of 2 normals")
              }
            } else {
              trainModel <- .getKernelEmpiricalParam
              if (verbose) {
                  message("  ==> select Gaussian kernel-based empirical")
              }
            }
        }
        # actual training with the chosen model
        gp <- trainModel(rc, "LR correlation (null)",
                         verbose = verbose,
                         file.name = file.name
        )
        parameters(obj)$LR.0$model <- gp
        
        # RT correlation null ------------------------------------
        
        if (parameters(obj)$quick) {
          # RT correlations are assumed to be equal to LR correlations
          if (verbose) {
              message("Quick learning, receptor-target correlation null",
                  " distribution assumed to be equal to ligand-receptor...")
          }
          parameters(obj)$RT.0$n <- parameters(obj)$LR.0$n
          parameters(obj)$RT.0$model <- parameters(obj)$LR.0$model
        } else {
          # RT correlations are actually learnt
          if (verbose) {
              message("Learning receptor-target",
                  " correlation null distribution...")
          }
          ds.RT.null <- .getEmpiricalNull(obj)
          
          t <- ds.RT.null[[1]]
          if (length(ds.RT.null) > 1) {
            for (i in 2:length(ds.RT.null)) t <- rbind(t, ds.RT.null[[i]])
          }
          above <- unlist(strsplit(t$target.corr, split = "\\|"))
          r.corrs <- NULL
          for (i in seq_len(length(above))) {
            corr <- as.numeric(strsplit(above[i], split = ";")[[1]])
            r.corrs <- c(r.corrs, corr)
          }
          if (null.model == "stable") {
            # sub-sample randomized R-T correlations to limit compute time
            r.corrs <- sample(r.corrs, parameters(obj)$LR.0$n)
          }
          parameters(obj)$RT.0$n <- length(r.corrs)
          
          # fit null model
          if (!is.null(plot.folder)) {
            file.name <- paste0(plot.folder, "/", filename, "_RT-null.pdf")
          } else {
            file.name <- NULL
          }
          gp <- trainModel(r.corrs, "RT correlation (null)",
                           verbose = verbose, file.name = file.name
          )
          parameters(obj)$RT.0$model <- gp
        }
        
        if (verbose) {
            message("Learning of statistical model parameters completed")
        }
        obj
    }
    
) # learnParameters


# Scoring of gene signatures in a BSRSignature object ==========================

setGeneric("scoreLRGeneSignatures", signature="obj",
    function(obj, ...) standardGeneric("scoreLRGeneSignatures")
)
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
#' data(bsrdm, package = "BulkSignalR")
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redBP <- reduceToBestPathway(bsrinf)
#' bsrsig.redBP <- BSRSignature(bsrinf.redBP, qval.thres = 0.001)
#' res <-scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
#'     name.by.pathway = FALSE
#' )
#' @importFrom foreach %do% %dopar%
#' @importFrom methods is
#' @importFrom matrixStats rowMeans2 colSums2
setMethod("scoreLRGeneSignatures", "BSRDataModel", function(obj,
    sig, LR.weight = 0.5, robust = FALSE,
    name.by.pathway = FALSE, abs.z.score = FALSE,
    rownames.LRP = FALSE) {

    if (!is(sig, "BSRSignature")) {
        stop("sig must be a BSRSignature object")
    }
    if (LR.weight <= 0 || LR.weight >= 1) {
        stop("LRweight must reside in (0;1)")
    }


    if (initialOrganism(obj) != "hsapiens") {
        all.genes <- unlist(initialOrthologs(obj))
    } else {
        all.genes <- rownames(ncounts(obj))
    }

    # intersect signature gene names with RNA-seq data
    ncounts <- ncounts(obj)

    ligands <- vector("list", length(ligands(sig)))
    receptors <- vector("list", length(receptors(sig)))
    tg.genes <- vector("list", length(tgGenes(sig)))


    for (i in seq_along(ligands(sig))) {
        ligands[[i]] <- intersect(ligands(sig)[[i]], all.genes)
    }
    for (i in seq_along(receptors(sig))) {
        receptors[[i]] <- intersect(receptors(sig)[[i]], all.genes)
    }
    for (i in seq_along(tgGenes(sig))) {
        tg.genes[[i]] <- intersect(tgGenes(sig)[[i]], all.genes)
    }

    good <- lengths(ligands) > 0 &
        lengths(receptors) > 0 &
        lengths(tg.genes)  > 0

    ligands <- ligands[good]
    receptors <- receptors[good]
    tg.genes <- tg.genes[good]
    pathways <- pathways(sig)[good]

    # scale ncounts
    if (logTransformed(obj)) {
        ncounts <- 2**ncounts
    }
    if (robust) {
        z <- (ncounts - apply(ncounts, 1, stats::median)) 
        z <- z / apply(ncounts, 1, stats::mad)
    } else {
        z <- (ncounts - matrixStats::rowMeans2(ncounts)) / 
        apply(ncounts, 1, stats::sd)
    }
    if (abs.z.score) {
        z <- abs(z)
    }

    if (initialOrganism(obj) != "hsapiens") {
        rownames(z) <- all.genes
    }

    # compute the LR gene signatures
    i <- NULL
    pwn <- foreach::foreach(i = seq_len(length(pathways)), .combine = c) %do% {
        if (name.by.pathway) {
            if (rownames.LRP) {
                paste0(
                    "{", paste(ligands[[i]], collapse = ";"), "} / {",
                    paste(receptors[[i]], collapse = ";"), "} | ", pathways[[i]]
                )
            } else {
                pathways[[i]]
            }
        } else if (!name.by.pathway) {
            paste0(
                "{", paste(ligands[[i]], collapse = ";"), "} / {",
                paste(receptors[[i]], collapse = ";"), "}"
            )
        }
    }

    res <- matrix(0, nrow = length(pathways), ncol = ncol(ncounts), 
        dimnames = list(pwn, colnames(ncounts)))
    for (i in seq_len(length(pathways))) {
        # average ligand z-score
        zz <- z[ligands[[i]], ]
        if (is.matrix(zz)) {
            mL <- matrixStats::colSums2(zz) / length(ligands[[i]])
        } else {
            mL <- zz
        }

        # average receptor z-score
        zz <- z[receptors[[i]], ]
        if (is.matrix(zz)) {
            mR <- matrixStats::colSums2(zz) / length(receptors[[i]])
        } else {
            mR <- zz
        }

        # average target gene z-score
        zz <- z[tg.genes[[i]], ]
        if (is.matrix(zz)) {
            mT <- matrixStats::colSums2(zz) / length(tg.genes[[i]])
        } else {
            mT <- zz
        }

        res[i, ] <- LR.weight * 0.5 * (mL + mR) + (1 - LR.weight) * mT
    }

    res
}) # scoreLRGeneSignatures
