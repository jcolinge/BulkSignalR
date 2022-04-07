#' @title Prepare a BSRDataModel object from expression data
#'
#' @description Take a matrix or data frame containing RNA sequencing,
#'   microarray, or expression proteomics data and returns a BSRDataModel
#'   object ready for subsequent training. Normally, BSRDataModel objects
#'   are not instantiated directly, but through this function.
#'
#' @param counts     A table or matrix of read counts.
#' @param normalize  A logical indicating whether \code{counts} should be
#'   normalized according to \code{method} or if it was normalized beforehand.
#' @param symbol.col The index of the column containing the gene symbols in case
#'   those are not the row names of \code{counts} already.
#' @param min.count  The minimum read count of a gene to be considered expressed
#'   in a sample.
#' @param prop       The minimum proportion of samples where a gene must be
#'   expressed higher than \code{min.count} to keep that gene.
#' @param method     The normalization method ('UQ' for upper quartile or 'TC'
#'   for total count).
#' @param log.transformed  A logical indicating whether expression data were
#'   already log2-transformed, e.g., some microarray data.
#' @param min.LR.found  The minimum number of ligands or receptors found in
#'   \code{count} row names after eliminating the rows containing too many
#'   zeros according to \code{min.count} and \code{prop}.
#'
#' @return A BSRModelData object with empty model parameters.
#'
#' @details The \code{counts} matrix or table should be provided with expression
#'   levels of protein coding genes in each samples (column) and
#'   \code{rownames(counts)} set to HUGO offcial gene symbols. For commodity, it
#'   is also possible to provide \code{counts} with the
#'   gene symbols stored in one of its columns. This column must be specified
#'   with \code{symbol.col}. In such a case, \code{prepareDataset} will extract
#'   this column and use it to set the row names. Because row names must be
#'   unique, \code{prepareDataset} will eliminate rows with duplicated gene
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
#'   \code{normalize} to \code{FALSE}.
#'
#'   In case proteomic or microarray data are provided, \code{min.count} must be
#'   understood as its equivalent with respect to those data.
#' @export
#' @examples
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
prepareDataset <- function(counts, normalize = TRUE, symbol.col = NULL, min.count = 10,
    prop = 0.1, method = c("UQ", "TC"), log.transformed = FALSE, min.LR.found = 80) {

    if (prop < 0 || prop > 1)
        stop("prop must lie in [0;1]")
    if (min.count < 0)
        stop("min.count must be positive")
    if (normalize)
        method <- match.arg(method)
    else
        if (nchar(method) == 0)
            stop(paste0("In case of user-normalized counts, the name of the ",
            "normalization must be documented through the parameter 'method'"))

    if (!is.null(symbol.col)) {
        if (!is.numeric(symbol.col))
            stop("symbol.col must be the index of the column containing the gene symbols")

        # simple but desperately slow counts <-
        # aggregate(.~symbol,data=counts,FUN=max)

        # home-made but fast
        symbols <- as.character(counts[, symbol.col])
        d <- symbols[duplicated(symbols)]
        bad <- NULL
        for (s in d) {
            i <- which(symbols == s)
            t <- rowSums(counts[i, -symbol.col])
            bad <- c(bad, i[-which.max(t)])
        }

        # remove duplicates and the gene symbol column
        counts <- counts[-bad, -symbol.col]
        rownames(counts) <- symbols[-bad]
    }

    if (is.null(rownames(counts)) || typeof(rownames(counts)) != "character")
        stop("The read count matrix must be provided with gene symbols as row names")

    # as of now we ensure that counts is a matrix
    if (!is.matrix(counts))
        counts <- data.matrix(counts)

    # avoid empty rows even if no normalization is performed here
    counts <- counts[rowSums(counts) > 0, ]

    if (normalize) {
        good.c <- rowSums(counts >= min.count) >= prop * ncol(counts)
        counts <- counts[good.c, ]
        if (method == "UQ")
            tot <- apply(counts, 2, function(x) stats::quantile(x[x > 0], prob = 0.75))
        else
            tot <- colSums(counts)
        ncounts <- sweep(counts, 2, tot/stats::median(tot), "/")
    }
    else
        ncounts <- counts

    nLR <- length(intersect(
        c(SingleCellSignalR::LRdb$ligand, SingleCellSignalR::LRdb$receptor),
        rownames(ncounts)))
    if (nLR < min.LR.found)
        stop(paste0("Not enough LR genes (",nLR," < ", min.LR.found,
                    " were found).\n"))

    new("BSRDataModel", ncounts=ncounts, log.transformed=log.transformed,
        normalization=toupper(method))

}  # prepareDataset
