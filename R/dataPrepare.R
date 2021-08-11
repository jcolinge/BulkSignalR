#' Prepare a BulkSignalR data set
#'
#' Take a matrix or data frame containing read counts along with sample types and combine them
#' into a list for subsequent use with BulkSignalR. Expression proteomics data could be used as well.
#'
#' @param counts     A table or matrix of read counts.
#' @param types      A vector of characters defining sample types.
#' @param normalize  A logical indicating whether \code{counts} should be normalized according to \code{method} or it was normalized beforehand.
#' @param symbol.col The index of the column containing the gene symbols in case those are not the rownames of \code{counts} already.
#' @param min.count  The minimum read count of a gene to be considered expressed in a sample.
#' @param prop       The minimum proportion of samples where a must be expressed to keep that gene.
#' @param method     The normalization method ("UQ" for upper quartile or "TC" for total count).
#' @param log.transformed  A logical indicating whether expression data were already log-transformed, e.g., some microarray data.
#' @return A list containing the read counts and the sample types for further use in BulkSignalR.
#'
#' The \code{counts} matrix or table should be provided with expression levels of protein coding genes in each samples (column)
#' and \code{rownames(counts)} set to the gene symbols. For commodity, it is also possible to provide
#' the expression matrix or table \code{counts} with the gene symbols stored in one of the columns.
#' This column must be specified with \code{symbol.col}. In such a case, \code{prepareDataset} will extract
#' this column and use it to set the row names. Because row names must be unique, \code{prepareDataset}
#' will eliminate rows with duplicated gene symbols by keeping the rows with maximum average expression.
#' Gene symbol duplication typically occurs in protein coding genes after genome alignment due to
#' errors in genome feature annotation files (GTF/GFF), where a handful of deprecated gene annotations may remain
#' along with the current ones, or some genes are not given their fully specific symbols. If your
#' read count extraction pipeline does not take care of this phenomenon, the maximum mean expression
#' selection strategy implemented here should solve this difficulty for the sake of inferring ligand-receptor
#' interactions.
#'
#' If \code{normalize} is \code{TRUE} then normalization is performed according to \code{method}.
#' If those two simple methods are not satisfying, then it is possible to provide a pre-normalized
#' matrix, e.g., using edgeR TMM algorithm.
#'
#' In case proteomic data are provided,
#' \code{min.count} must be understood as its equivalent with respect to those data (spectral counts, iBAQ, etc.).
#' @export
#' @examples
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
prepareDataset <- function(counts,types,normalize=TRUE,symbol.col=NULL,min.count=10,prop=0.1,method=c("UQ","TC"),log.transformed=FALSE){

  if (prop<0 || prop>1)
    stop("prop must lie in [0;1]")
  if (min.count<0)
    stop("min.count must be positive")
  method <- match.arg(method)

  if (!is.null(symbol.col)){
    if (!is.numeric(symbol.col))
      stop("symbol.col must be the index of the column containing the gene symbols")

    # simple but desperately slow
    # counts <- aggregate(.~symbol,data=counts,FUN=max)

    # home-made but fast
    symbols <- as.character(counts[,symbol.col])
    d <- symbols[duplicated(symbols)]
    bad <- NULL
    for (s in d){
      i <- which(symbols==s)
      t <- rowSums(counts[i,-symbol.col])
      bad <- c(bad,i[-which.max(t)])
    }

    # remove duplicates and the gene symbol column
    counts <- counts[-bad,-symbol.col]
    rownames(counts) <- symbols[-bad]
  }

  if (length(types)!=ncol(counts))
    stop("The length of types does not match the read count matrix number of columns")
  if (is.null(rownames(counts)) || typeof(rownames(counts))!="character")
    stop("The read count matrix must be provided with gene symbols as row names")

  if (normalize){
    good.c <- apply(counts,1,function(x) sum(x>=min.count))>=prop*ncol(counts)
    counts <- counts[good.c,]
    if (method=="UQ")
      tot <- apply(counts,2,function(x) stats::quantile(x[x>0],prob=0.75))
    else
      tot <- colSums(counts)
    ncounts <- sweep(counts,2,tot/stats::median(tot),"/")
  }
  else
    ncounts <- counts

  list(ncounts=ncounts,types=types,log.transformed=log.transformed)

} # prepareDataset
