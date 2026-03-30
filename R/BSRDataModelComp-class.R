#' BulkSignalR Data Model Compare Object
#'
#' An S4 class to represent the expression data used for inferring
#' ligand-receptor interactions based on sample cluster comparisons.
#'
#' @slot comp   A named list of BSRClusterComp objects, one per
#' comparison.
#' @slot mu     A number representing the average value in the normalized and
#' lop1p-transformed gene expression matrix. This value is used to compute
#' the LR-score (cf. SingleCellSignalR paper, Cabello-Aguilar, et al.,
#' Nucleic Acids Res, 2020)
#' @export
#' @examples
#' new("BSRDataModelComp")
#'
setClass("BSRDataModelComp",
    contains = c("BSRDataModel"),
    slots = c(
        comp = "list",
        mu = "numeric"
    ),
    prototype = list(
        initial.organism = "hsapiens",
        initial.orthologs = list("A", "B", "C"),
        ncounts = matrix(1.0, nrow = 2, ncol = 1,
            dimnames = list(c("A", "B"), "C")),
        log.transformed = FALSE,
        normalization = "UQ",
        param = list(spatial.smooth = FALSE),
        comp = list(),
        mu = 1.0
    )
)

setValidity(
    "BSRDataModelComp",
    function(object) {
        if (length(object@comp) > 0) {
            if (is.null(names(object@comp))) {
                return("comp names must be set")
            }
            for (c in names(object@comp)) {
                if (!is(object@comp[[c]], "BSRClusterComp")) {
                    return("comp must contain objects of class BSRClusterComp")
                }
            }
        }
        if (!is.numeric(object@mu) || object@mu <= 0) {
            return("mu must be numeric and >0")
        }

        TRUE
    }
)

setMethod(
    "show", "BSRDataModelComp",
    function(object) {
        callNextMethod()
      cat("mu: ", object@mu, "\n", sep = "")
      cat("Defined comparisons:\n")
      utils::str(object@comp)
    }
)


# Accessors & setters ========================================================

setGeneric("comparison", signature="x",
    function(x) standardGeneric("comparison")
)
#' Comparisons list accessor
#'
#' @name comparison
#' @aliases comparison,BSRDataModelComp-method
#' @param x object BSRDataModelComp
#' @return comp
#' @examples
#' bsrdm.comp <- new("BSRDataModelComp")
#' comparison(bsrdm.comp)
#' @export
setMethod("comparison", "BSRDataModelComp", function(x) x@comp)


setGeneric("comparison<-", signature=c("x", "value"),
    function(x, value) standardGeneric("comparison<-")
)
#' Comparisons list setter (internal use only, use addComparison() otherwise)
#'
#' @param x object BSRDataModelComp
#' @param value value to be set for BSRDataModelComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("comparison<-", "BSRDataModelComp", function(x, value) {
    x@comp <- value
    methods::validObject(x)
    x
})


setGeneric("mu", signature="x",
    function(x) standardGeneric("mu")
)
#' Mu accessor
#'
#' @name mu
#' @aliases mu,BSRDataModelComp-method
#' @param x object BSRDataModelComp
#' @return mu
#' @examples
#' bsrdm.comp <- new("BSRDataModelComp")
#' mu(bsrdm.comp)
#' @export
setMethod("mu", "BSRDataModelComp", function(x) x@mu)


setGeneric("mu<-", signature=c("x", "value"),
    function(x, value) standardGeneric("mu<-")
)
#' Mu setter (internal use only)
#'
#' @param x object BSRDataModelComp
#' @param value value to be set for BSRDataModelComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("mu<-", "BSRDataModelComp", function(x, value) {
    x@mu <- value
    methods::validObject(x)
    x
})


#' Convert BSRDataModel to BSRDataModelComp
#' 
#' Enable the promotion of a BSRDataModel object into
#' a BSRDataModelComp object by
#' coercion.
#' 
#' @name coerce
#' @aliases coerce,BSRDataModel,BSRDataModelComp-method
#' @param from  BSRDataModel object  
#' @return A BSRDataModelComp object
#' @examples
#' bsrdm <- new("BSRDataModel")
#' bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
#'
#' @importFrom methods is new
#' @exportMethod coerce
setAs("BSRDataModel", "BSRDataModelComp", function(from) {

    if (!is(from, "BSRDataModel")) {
        stop("bsrdm must be of class BSRDataModel")
    }
    m <- mean(ncounts(from))
    if (!logTransformed(from)) {
        m <- log1p(m) / log(2)
    } # approximate of mu on the log2-transformed matrix

    new("BSRDataModelComp", from, 
        comp = list(),
        mu = m) 
})


setGeneric("addClusterComp", signature="obj",
    function(obj, ...) standardGeneric("addClusterComp")
)
#' Add a comparison between two clusters of samples
#'
#' Add a comparison to a BSRDataModelComp object.
#'
#' @name addClusterComp
#' @aliases addClusterComp,BSRDataModelComp-method
#'
#' @param obj    A BSRDataModelComp object output by
#'   \code{\link{setAs}}.
#' @param cmp   A BSRClusterComp object to add.
#' @param cmp.name  The name of the comparison to add.
#'
#' @details Add \code{cmp} to the list of comparisons contained in
#' \code{obj}.
#'
#' @return A BSRDataModelComp object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc, package = "BulkSignalR")
#' normal <- grep("^N", names(sdc))
#' bsrdm <- BSRDataModel(sdc[, -normal])
#'
#' # define the comparison
#' bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
#' colA <- as.integer(1:3)
#' colB <- as.integer(12:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = runif(n, 0, 10)
#' )
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- BSRClusterComp(bsrdm.comp, colA, colB, stats)
#'
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' @importFrom methods new
setMethod("addClusterComp", "BSRDataModelComp", function(obj, cmp,
    cmp.name) {
    if (!is(cmp, "BSRClusterComp")) {
        stop("cmp must be of class BSRClusterComp")
    }
    if (!is.character(cmp.name)) {
        stop("cmp.name must be of type character")
    }
    if (length(cmp.name) == 0) {
        stop("cmp.name must have length > 0")
    }
    if (cmp.name %in% names(comparison(obj))) {
        stop("cmp.name is already in the list of comparisons")
    }

    tmp <- c(comparison(obj), list(cmp))
    names(tmp)[length(tmp)] <- cmp.name
    comparison(obj) <- tmp
    obj
    
}) # addClusterComp


setGeneric("removeClusterComp", signature="obj",
    function(obj, ...) standardGeneric("removeClusterComp")
)
#' Remove a comparison from a BSRDataModelComp object.
#'
#' @name removeClusterComp
#' @aliases removeClusterComp,BSRDataModelComp-method
#'
#' @param obj    A BSRDataModelComp object output by
#'   \code{\link{setAs}}.
#' @param cmp.name  The name of the comparison to remove.
#'
#' @details Remove the comparison with \code{cmp.name} from the list of
#' comparisons contained in \code{obj}.
#'
#' @return A BSRDataModelComp object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc, package = "BulkSignalR")
#' normal <- grep("^N", names(sdc))
#' bsrdm <- BSRDataModel(sdc[, -normal])
#'
#' # define the comparison
#' bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
#' colA <- as.integer(1:3)
#' colB <- as.integer(12:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = runif(n, 0, 10)
#' )
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- BSRClusterComp(bsrdm.comp, colA, colB, stats)
#'
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' bsrdm.comp <- removeClusterComp(bsrdm.comp, "random.example")
#'
setMethod("removeClusterComp", "BSRDataModelComp", function(obj, cmp.name) {
    if (!is.character(cmp.name)) {
        stop("cmp.name must be of type character")
    }
    if (length(cmp.name) == 0) {
        stop("cmp.name must have length > 0")
    }
    if (!(cmp.name %in% names(comparison(obj)))) {
        stop("cmp.name must be in the list of comparisons")
    }

    if (length(comparison(obj)) == 1) {
        comparison(obj) <- list()
    } else {
        comparison(obj) <- comparison(obj)[names(comparison(obj)) != cmp.name]
    }

    obj
    
}) # removeClusterComp


# Constructor ==================================================================

#' Constructor of the BSRDataModelComp class
#'
#' Take a matrix or data frame containing RNA sequencing,
#' microarray, or expression proteomics data as input parameter
#' and return a BSRDataModelComp
#' object ready for subsequent use.
#' 
#' @param counts     A table or matrix of read counts.
#' @param species    Data were obtained for this organism.
#' @param normalize  A logical indicating whether \code{counts} should be
#'   normalized according to \code{method} or if it was normalized beforehand.
#' @param symbol.col The index of the column containing the gene symbols in case
#'   those are not the row names of \code{counts} already.
#' @param min.count  The minimum read count of a gene to be considered expressed
#'   in a sample.
#' @param prop       The minimum proportion of samples where a gene must be
#'   expressed higher than \code{min.count} to keep that gene.
#' @param method     The normalization method ('UQ' for upper quartile or 'TC'
#'   for total count). If \code{normalize==FALSE}, then method must be
#'   used to document the name of the normalization method applied by the user.
#' @param UQ.pc      Percentile for upper-quartile normalization, number
#' between 0 and 1 (in case the default 0.75 - hence the name - is not
#' appropriate).
#' @param log.transformed  A logical indicating whether expression data were
#'   already log2-transformed, e.g., some microarray data.
#' @param min.LR.found  The minimum number of ligands or receptors found in
#'   \code{count} row names after eliminating the rows containing too many
#'   zeros according to \code{min.count} and \code{prop}.
#' @param conversion.dict  Correspondence table of HUGO gene symbols
#' human/nonhuman. Not used unless the organism is different from human.
#' @param x.col In a SpatialExperiment object, the index of the column
#' containing the x coordinates in the dafaframe returned by rowData(), usually 
#' named array_row.
#' @param y.col  In a SpatialExperiment object, the index of the column
#' containing the y coordinates in the data.frame returned by rowData(), usually 
#' named array_col.
#' @param barcodeID.col   In a SpatialExperiment object, the index of the column
#' containing the barcodeID in the dafaframe returned by colData(), usually
#' named barcode_id.
#'
#' @return A BSRModelDataComp object with empty model parameters.
#'
#' @details The \code{counts} matrix or table should be provided with expression
#'   levels of protein coding genes in each samples (column) and
#'   \code{rownames(counts)} set to HUGO official gene symbols.
#'   For commodity, it is also possible 
#'   to provide \code{counts} with the
#'   gene symbols stored in one of its columns. This column must be specified
#'   with \code{symbol.col}. In such a case, \code{BSRDataModel} will extract
#'   this column and use it to set the row names. Because row names must be
#'   unique, \code{BSRDataModel} will eliminate rows with duplicated gene
#'   symbols by keeping the rows with maximum average expression. Gene symbol
#'   duplication may occur in protein coding genes after genome alignment
#'   due to errors in genome feature annotation files (GTF/GFF), where a handful
#'   of deprecated gene annotations might remain, or
#'   some genes are not given their fully specific symbols. If your read count
#'   extraction pipeline does not take care of this phenomenon, the maximum mean
#'   expression selection strategy implemented here should solve this difficulty
#'   for the sake of inferring ligand-receptor interactions.
#'
#'   If \code{normalize} is \code{TRUE} then normalization is performed
#'   according to \code{method}. If those two simple methods are not satisfying,
#'   then it is possible to provide a pre-normalized matrix setting
#'   \code{normalize} to \code{FALSE}. In such a case, the parameter
#'   \code{method} must be used to document the name of the normalization
#'   algorithm used.
#'
#'   In case proteomic or microarray data are provided, \code{min.count} must be
#'   understood as its equivalent with respect to those data types.
#'
#' @importFrom matrixStats rowMeans2 rowSums2 colSums2
#' @export
#' @examples
#' data(sdc, package = "BulkSignalR")
#' idx <- sample(nrow(sdc), 4000)
#' bsrdm.comp <- BSRDataModelComp(sdc[idx, c("N22","SDC17")],
#'                                normalize=FALSE, method="UQ")
BSRDataModelComp <- function(
    counts, normalize = TRUE, symbol.col = NULL, min.count = 10,
    prop = 0.1, method = c("UQ", "TC"), 
    log.transformed = FALSE, min.LR.found = 80,
    species = "hsapiens", conversion.dict = NULL,
    UQ.pc = 0.75,x.col = NULL, y.col = NULL,
    barcodeID.col = NULL) {
  
  bsrdm <- BSRDataModel(counts = counts,
                        normalize = normalize,
                        symbol.col = symbol.col,
                        min.count = min.count,
                        prop = prop,
                        method = method,
                        log.transformed = log.transformed,
                        min.LR.found = min.LR.found,
                        species = species,
                        conversion.dict = conversion.dict,
                        UQ.pc = UQ.pc,
                        x.col = x.col,
                        y.col = y.col,
                        barcodeID.col = barcodeID.col)
  as(bsrdm, "BSRDataModelComp")
  
} # BSRDataModelComp


# Scoring of gene signatures in a BSRSignature object ==========================

#' Score ligand-receptor gene signatures
#'
#' Compute ligand-receptor gene signature scores over a BSRDataModelComp
#' specific comparison.
#'
#' @name scoreLRGeneSignatures
#' @aliases scoreLRGeneSignatures,BSRDataModelComp-method
#'
#' @param obj           A BSRDataModelComp object.
#' @param sig           A BSRSignatureComp object.
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
#' # prepare data
#' data(bsrdm.comp, package = "BulkSignalR")
#' data(bsrinf.comp, package = "BulkSignalR")
#'
#' # reduction
#' bsrinf.red <- reduceToBestPathway(bsrinf.comp)
#' # signature extraction and scoring
#' bsrsig.red <- BSRSignatureComp(bsrinf.red, qval.thres = 1e-6)
#' scores.red <- scoreLRGeneSignatures(bsrdm.comp, bsrsig.red,
#'     name.by.pathway = TRUE, rownames.LRP = TRUE
#' )
#' @importFrom foreach %do% %dopar%
#' @importFrom methods is
#' @importFrom matrixStats rowMeans2 colSums2
setMethod("scoreLRGeneSignatures", "BSRDataModelComp", function(obj,
    sig, LR.weight = 0.5, robust = FALSE,
    name.by.pathway = FALSE, abs.z.score = FALSE,
    rownames.LRP = FALSE) {
    
    if (!is(sig, "BSRSignatureComp")) {
        stop("sig must be a BSRSignature object")
    }
    if (LR.weight <= 0 || LR.weight >= 1) {
        stop("LRweight must reside in (0;1)")
    }

    # retrieve the BSRClusterComp object
    cmp.name <- comparisonName(sig)
    if (!(cmp.name %in% names(comparison(obj)))) {
        stop("The comparison name in sig is not present in obj")
    }
    cc <- comparison(obj)[[cmp.name]]

    # species management
    if (initialOrganism(obj) != "hsapiens") {
        all.genes <- unlist(initialOrthologs(obj))
    } else {
        all.genes <- rownames(ncounts(obj))
    }

    # get the ncount matrix with proper columns
    ncounts <- ncounts(obj)[, c(colClusterA(cc), colClusterB(cc))]

    # intersect signature gene names with RNA-seq data
    ligands <- list()
    receptors <- list()
    tg.genes <- list()

    for (i in seq_along(ligands(sig))) {
        ligands[[i]] <- intersect(ligands(sig)[[i]], all.genes)
    }
    for (i in seq_along(receptors(sig))) {
        receptors[[i]] <- intersect(receptors(sig)[[i]], all.genes)
    }
    for (i in seq_along(tgGenes(sig))) {
        tg.genes[[i]] <- intersect(tgGenes(sig)[[i]], all.genes)
    }

    good <- vapply(ligands, length, integer(1)) > 0 &
        vapply(receptors, length, integer(1)) > 0 &
        vapply(tg.genes, length, integer(1)) > 0

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

    res <- matrix(0, nrow = length(pathways),
        ncol = ncol(ncounts), dimnames = list(pwn, colnames(ncounts)))
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
