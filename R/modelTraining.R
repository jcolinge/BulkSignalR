#' Internal function to generate expression matrix permutation indices
#'
#' @param ncounts         A matrix of normalized read counts.
#' @param n.bins          Number of bins.
#'
#' @return A list containing a vectors of row indices. The vector at index
#'   \code{i} in the list contains the row indices of rows with mean normalized
#'   read count in bin \code{i}.
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom stats quantile
.buildPermutationIndices <- function(ncounts, n.bins = 20) {

    rm <-  rowMeans(ncounts, na.rm = TRUE)
    breaks <- stats::quantile(rm, prob = (seq(0,n.bins))/n.bins)
    breaks[1] <- 0

    lapply(seq(2,length(breaks)),
           function(i) which(rm>breaks[i-1] & rm<=breaks[i])
    )

}  # .buildPermutationIndices


#' Internal function to shuffle permutation indices
#'
#' @param pind      Permutation indices such as returned by
#'   \code{\link{.buildPermutationIndices}}.
#'
#' @return A list with same structure as \code{pind} with shuffled indices
#'   within each bin.
#' @importFrom foreach %do% %dopar%
#'
.shufflePermutationIndices <- function(pind) {

    lapply(pind,
           function(x) sample(x, length(x))
    )

} # .shufflePermutationIndices


#' Internal function to generate a randomized expression matrix
#'
#' @param ncounts    A matrix of normalized read counts.
#' @param pind      Permutation indices such as returned by
#'   \code{\link{.buildPermutationIndices}}.
#'
#' @return \code{ncount} with shuffled row names (gene symbols). Shuffling is
#'   performed within rows of comparable average expression.
#'
#'
.buildPermutatedCountMatrix <- function(ncounts, pind) {

    symbols <- rownames(ncounts)
    rind <- .shufflePermutationIndices(pind)

    for (i in seq_len(length(pind)))
        symbols[pind[[i]]] <- symbols[rind[[i]]]
    rownames(ncounts) <- symbols

    ncounts

}  # .buildPermutatedCountMatrix


#' Internal function to fit a Gaussian distribution
#'
#' @description Maximum-likelihood estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param file.name   The file name of a PDF file.
#'
#' @return A list with the mean (\code{mu}) and standard deviation
#'   (\code{sigma}) estimates.
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
.getGaussianParam <- function(d, title, file.name = NULL) {
    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq = FALSE, main = title,
                       xlab = "Spearman correlation")
    }

    mu <- mean(d)
    sigma <- stats::sd(d)
    if (!is.null(file.name)) {
        x <- seq(-1, 1, by = 0.01)
        graphics::lines(x = x, y = stats::dnorm(x, mu,
                                                sigma), col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Normal",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
    }
    list(mu = mu, sigma = sigma)

}  # .getGaussianParam


#' Sampling of correlations downstream the receptors null distribution
#'
#' Perform receptor downstream analysis with
#' \code{.checkReceptorSignaling} based on randomized expression data and
#' ligand-receptor pairs selected from the same randomized data.
#'
#' @param ncounts         A matrix or table of normalized read counts.
#' @param n.rand          The number of repetitions.
#' @param min.cor         The minimum ligand-receptor Spearman correlation
#'   required.
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#'
#' @return A list of \code{n.rand} tables such as output by
#'   \code{.checkReceptorSignaling}. Each table is computed from a randomized
#'   expression matrix (randomized \code{ncounts}).
#'
#' @details A large number of correlations (ligand-receptor and receptor-downstream
#'   target genes) is reported in each randomized matrix. Therefore,
#'   \code{n.rand} should be
#'   given a modest value to avoid unnecessarily long computations.
#'
#'   See \code{\link{.checkReceptorSignaling}} for more details about the
#'   parameters.
#'
#' @importFrom foreach %do% %dopar%
#'
.getEmpiricalNull <- function(ncounts, n.rand = 5, min.cor = -1,
                             with.complex = TRUE, max.pw.size = 200,
                             min.pw.size = 5, min.positive = 4) {

    pindices <- .buildPermutationIndices(ncounts)
    r.ds <- prepareDataset(ncounts, normalize = FALSE, method = "ALREADY")
    if (foreach::getDoParWorkers() > 1)
        foreach::foreach(k = seq_len(n.rand), .combine = c) %dopar% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices)
            r.LR <- .getCorrelatedLR(r.ds, min.cor = min.cor)
            list(.checkReceptorSignaling(r.ds, r.LR,
                        with.complex = with.complex, max.pw.size = max.pw.size,
                        min.pw.size = min.pw.size, min.positive = min.positive)
            )
        }
    else
        foreach::foreach(k = seq_len(n.rand), .combine = c) %do% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices)
            r.LR <- .getCorrelatedLR(r.ds, min.cor = min.cor)
            list(.checkReceptorSignaling(r.ds, r.LR,
                         with.complex = with.complex, max.pw.size = max.pw.size,
                         min.pw.size = min.pw.size, min.positive = min.positive)
            )
        }

}  # .getEmpiricalNull


#' Sampling of ligand-receptor correlation null distribution
#'
#' Perform a ligand-receptor Spearman correlation analysis based
#' on randomized expression data.
#'
#' @param ncounts         A matrix or table of normalized read counts.
#' @param n.rand          The number of repetitions.
#' @param min.cor         The minimum ligand-receptor correlation required.
#'
#' @return A list of \code{n.rand} tables such as output by
#'   \code{\link{.getCorrelatedLR}}. Each table is computed from a randomized
#'   expression matrix (randomized \code{ncounts}).
#'
#' @details A large number of correlations is reported in each randomized
#'   matrix. Therefore,
#'   \code{n.rand} should be given a modest value to avoid unnecessarily long
#'   computations.
#'
#'   See \code{\link{.getCorrelatedLR}} for more details about the parameters.
#'
#' @importFrom foreach %do% %dopar%
.getEmpiricalNullCorrLR <- function(ncounts, n.rand = 5, min.cor = -1) {

    pindices <- .buildPermutationIndices(ncounts)
    r.ds <- prepareDataset(ncounts, normalize = FALSE, method = "ALREADY")

    if (foreach::getDoParWorkers() > 1)
        foreach::foreach(k = seq_len(n.rand), .combine = 'c',
            .packages="BulkSignalR"
            ) %dopar% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices)
            list(.getCorrelatedLR(r.ds, min.cor = min.cor))
        }
    else
    foreach::foreach(k = seq_len(n.rand), .combine = c) %do% {
        ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices)
        list(.getCorrelatedLR(r.ds, min.cor = min.cor))
    }

}  # .getEmpiricalNullCorrLR

