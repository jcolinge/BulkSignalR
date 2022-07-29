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
#' @param seed      Seed to reproduce exact same sampling.
#'
#' @return A list with same structure as \code{pind} with shuffled indices
#'   within each bin.
#' @importFrom foreach %do% %dopar%
#'
.shufflePermutationIndices <- function(pind,seed) {

    set.seed(seed)
    lapply(pind,
           function(x) sample(x, length(x))
    )

} # .shufflePermutationIndices


#' Internal function to generate a randomized expression matrix
#'
#' @param ncounts    A matrix of normalized read counts.
#' @param pind      Permutation indices such as returned by
#'   \code{\link{.buildPermutationIndices}}.
#' @param seed      Seed to reproduce exact same sampling.
#'
#' @return \code{ncount} with shuffled row names (gene symbols). Shuffling is
#'   performed within rows of comparable average expression.
#'
#'
.buildPermutatedCountMatrix <- function(ncounts, pind,seed) {

    symbols <- rownames(ncounts)
    rind <- .shufflePermutationIndices(pind,seed)

    for (i in seq_len(length(pind)))
        symbols[pind[[i]]] <- symbols[rind[[i]]]
    rownames(ncounts) <- symbols

    ncounts

}  # .buildPermutatedCountMatrix


#' Internal function to fit a censored Gaussian distribution
#'
#' @description Maximum-likelihood estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param verbose   Provide details on computations.
#' @param file.name   The file name of a PDF file.
#'
#' @return A list with the mean (\code{mu}) and standard deviation
#'   (\code{sigma}) estimates.
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
.getGaussianParam <- function(d, title, verbose = FALSE, file.name = NULL) {

    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq=FALSE, main=paste0(title, " / censored normal"),
                       xlab = "Spearman correlation", breaks = 30)
    }

    # initial fit with a Gaussian over ]-infty;+infty[]
    mu <- mean(d)
    sigma <- stats::sd(d)
    if (verbose){
        cat("Initial estimate of the mean: ", mu, "\n", sep="")
        cat("Initial estimate of the standard deviation: ", sigma, "\n", sep="")
    }

    # ML fit of a censored Gaussian on [-1;1]
    GaussianLL <- function(par){
        q <- stats::pnorm(1, par[1], par[2]) - stats::pnorm(-1, par[1], par[2])
        -sum(stats::dnorm(d, par[1], par[2], log=TRUE)) + length(d)*log(q)
    }
    par.0 <- c(mu, sigma)
    lo.bound <- c(mu-0.5, 0.5*sigma)
    hi.bound <- c(mu+0.5, 1.5*sigma)
    res <- stats::optim(par.0, GaussianLL, method="L-BFGS-B",
                        lower=lo.bound, upper=hi.bound,
                        control=list(maxit=2000))
    if (res$convergence != 0)
        stop("optim() could not fit the normal distribution parameters")
    mu <- res$par[1]
    sigma <- res$par[2]
    if (verbose){
        cat("Censored normal mean: ", mu, "\n", sep="")
        cat("Censored normal standard deviation: ", sigma, "\n", sep="")
    }
    q <- stats::pnorm(1, mu, sigma) - stats::pnorm(-1, mu, sigma)
    start <- stats::pnorm(-1, mu, sigma)
    params <- list(mu = mu, sigma = sigma, factor = q,
                   start = start, distrib = "censored_normal")

    # KS test D statistics
    x <- seq(-1, 1, by = 0.005)
    y <- stats::dnorm(x, mu, sigma)/q
    params$D <- as.numeric(suppressWarnings(stats::ks.test(d, y)$statistic))

    # Chi2
    x <- seq(-1, 1, by = 0.05)
    h <- graphics::hist(d, breaks=x, plot=FALSE)
    hist.rf <- h$counts/length(d)
    gauss.rf <- .cdfGaussian(h$breaks[-1], params) -
        .cdfGaussian(h$breaks[-length(h$breaks)], params)
    params$Chi2 <- sum((hist.rf-gauss.rf)**2)

    # control plot
    if (!is.null(file.name)) {
        x <- seq(-1, 1, by = 0.002)
        graphics::lines(x = x, y = stats::dnorm(x, mu, sigma)/q,
                        col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Model",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
        # fn <- gsub("pdf$", "txt", file.name)
        # write.table(d, file=fn, row.names = FALSE)
    }

    params

}  # .getGaussianParam


#' Internal function to compute a censored Gaussian CDF
#'
#' @param x   A vector of observed values.
#' @param par A list containing the censored Gaussian model parameters.
#'
#' @return A vector of probabilities P(X<x|par).
.cdfGaussian <- function(x, par){

    (stats::pnorm(x, par$mu, par$sigma) - par$start) / par$factor

} # .cdfGaussian


#' Internal function to fit a censored mixed-Gaussian distribution
#'
#' @description Maximum-likelihood estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param verbose   Provide details on computations.
#' @param file.name   The file name of a PDF file.
#'
#' @return A list with the mean (\code{mu}) and standard deviation
#'   (\code{sigma}) estimates of each distribution along with the
#'   weight alpha applied to the first distribution.
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
.getMixedGaussianParam <- function(d, title, verbose = FALSE, file.name = NULL) {
    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq=FALSE,
                       main=paste0(title, " / censored mixed normal"),
                       xlab = "Spearman correlation", breaks = 30)
    }

    # ML fit of a censored mixed-Gaussian on [-1;1]
    mixedGaussianLL <- function(par){
        alpha <- par[1]
        mu1 <- par[2]
        sigma1 <- par[3]
        mu2 <- par[4]
        sigma2 <- par[5]
        q <- alpha*stats::pnorm(1, mu1, sigma1) +
            (1-alpha)*stats::pnorm(1, mu2, sigma2) -
            (alpha*stats::pnorm(-1, mu1, sigma1) +
                 (1-alpha)*stats::pnorm(-1, mu2, sigma2)
            )
        -sum(log(alpha*stats::dnorm(d, mu1, sigma1) +
                     (1-alpha)*stats::dnorm(d, mu2, sigma2))
        ) + length(d)*log(q)
    }
    mu <- mean(d)
    sigma <- stats::sd(d)
    par.0 <- c(0.7, mu-0.1, 0.75*sigma, mu+0.1, 3*sigma)
    lo.bound <- c(0.5, mu-0.5, 0.5*sigma, mu-0.2, 2*sigma)
    hi.bound <- c(1.0, mu+0.2, 1.5*sigma, mu+0.4, 10*sigma)
    res <- stats::optim(par.0, mixedGaussianLL, method="L-BFGS-B",
                        lower=lo.bound, upper=hi.bound,
                        control=list(maxit=2000))
    if (res$convergence != 0)
        stop("optim() could not fit the mixed normal distribution parameters")
    if (verbose)
        cat("Censored mixed normal parameters (alpha, mean1, sd1, mean2, sd2): ",
            paste(res$par, collapse=", "), "\n", sep="")
    alpha <- res$par[1]
    mu1 <- res$par[2]
    sigma1 <- res$par[3]
    mu2 <- res$par[4]
    sigma2 <- res$par[5]
    q <- alpha*stats::pnorm(1, mu1, sigma1) +
        (1-alpha)*stats::pnorm(1, mu2, sigma2) -
        (alpha*stats::pnorm(-1, mu1, sigma1) +
             (1-alpha)*stats::pnorm(-1, mu2, sigma2)
        )
    start <- alpha*stats::pnorm(-1, mu1, sigma1) +
        (1-alpha)*stats::pnorm(-1, mu2, sigma2)
    params <- list(alpha=alpha, mu1=mu1, sigma1=sigma1, mu2=mu2, sigma2=sigma2,
                   factor=q, start=start, distrib="censored_mixed_normal")

    # KS test D statistics
    x <- seq(-1, 1, by = 0.005)
    y <- alpha*stats::dnorm(x, mu1, sigma1) +
        (1-alpha)*stats::dnorm(x, mu2, sigma2)/q
    params$D <- as.numeric(suppressWarnings(stats::ks.test(d, y)$statistic))

    # Chi2
    x <- seq(-1, 1, by = 0.05)
    h <- graphics::hist(d, breaks=x, plot=FALSE)
    hist.rf <- h$counts/length(d)
    mixed.rf <- .cdfMixedGaussian(h$breaks[-1], params) -
        .cdfMixedGaussian(h$breaks[-length(h$breaks)], params)
    params$Chi2 <- sum((hist.rf-mixed.rf)**2)

    # control plot
    if (!is.null(file.name)) {
        x <- seq(-1, 1, by = 0.002)
        graphics::lines(x = x, y = alpha*stats::dnorm(x, mu1, sigma1) +
                            (1-alpha)*stats::dnorm(x, mu2, sigma2)/q,
                        col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Model",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
    }

    params

}  # .getMixedGaussianParam


#' Internal function to compute a censored mixed-Gaussian CDF
#'
#' @param x   A vector of observed values.
#' @param par A list containing the censored mixed-Gaussian model parameters.
#'
#' @return A vector of probabilities P(X<x|par).
.cdfMixedGaussian <- function(x, par){

    (par$alpha*stats::pnorm(x, par$mu1, par$sigma1) +
         (1-par$alpha)*stats::pnorm(x, par$mu2, par$sigma2) -
         par$start
    ) / par$factor

} # .cdfMixedGaussian


#' Internal function to fit an empirical distribution
#'
#' @description Based on stats::ecdf.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param verbose   Provide details on computations.
#' @param file.name   The file name of a PDF file.

#' @return A list with the step function implementing the CDF of
#'   the empirical distribution (\code{empirCDF}).
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
.getEmpiricalParam <- function(d, title, verbose = FALSE, file.name = NULL) {

    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq=FALSE, main=paste0(title, " / empirical"),
                       xlab = "Spearman correlation", breaks=30)
    }

    empir <- stats::ecdf(d)

    # control plot
    if (!is.null(file.name)) {
        step <- 0.005
        x <- seq(-1, 1, by = step)
        cd <- empir(x)
        n <- length(x)
        left <- cd[-c(n-1,n)]
        right <- cd[-c(1,2)]
        dens <- c(0, (right-left)/2/step, 0)
        graphics::lines(x = x, y = dens, col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Model",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
    }

    list(empirCDF = empir, distrib = "empirical")

}  # .getEmpiricalParam


#' Internal function to compute an empirical CDF
#'
#' @param x   A vector of observed values.
#' @param par A list containing the step function implementing the CDF.
#'
#' @return A vector of probabilities P(X<x|par).
.cdfEmpirical <- function(x, par){
    par$empirCDF(x)

} # .cdfEmpirical


#' Internal function to fit a Gaussian kernel-based empirical distribution
#'
#' @description Based on stats::density.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param verbose   Provide details on computations.
#' @param file.name   The file name of a PDF file.
#' @param n  The number of grid points for density FFT
#'
#' @return A list with the step function implementing the CDF of
#'   the empirical distribution (\code{kernelCDF}).
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
.getKernelEmpiricalParam <- function(d, title, verbose = FALSE,
                                     file.name = NULL, n=512) {

    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq=FALSE, main=paste0(title, " / kernel empirical"),
                       xlab = "Spearman correlation", breaks=30)
    }

    df <- stats::density(d, from=-1, to=1, n=n)
    cd <- cumsum(df$y)
    cd <- cd/cd[n]
    params <- list(kernelCDF = stats::stepfun(df$x, c(0, cd)),
                   distrib = "kernel_empirical")

    # KS test D statistics
    params$D <- as.numeric(suppressWarnings(stats::ks.test(d, df$y)$statistic))

    # Chi2
    x <- seq(-1, 1, by = 0.05)
    h <- graphics::hist(d, breaks=x, plot=FALSE)
    hist.rf <- h$counts/length(d)
    kernel.rf <- .cdfKernelEmpirical(h$breaks[-1], params) -
        .cdfKernelEmpirical(h$breaks[-length(h$breaks)], params)
    params$Chi2 <- sum((hist.rf-kernel.rf)**2)

    # control plot
    if (!is.null(file.name)) {
        # left <- cd[-c(n-1,n)]
        # right <- cd[-c(1,2)]
        # step <- 1/(n/2)
        # dens <- c(0, (right-left)/2/step, 0)
        # graphics::lines(x = df$x, y = dens, col = "blue", type = "l")
        graphics::lines(df, col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Model",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
    }

    params

}  # .getKernelEmpiricalParam


#' Internal function to compute a Gaussian kernel-based empirical CDF
#'
#' @param x   A vector of observed values.
#' @param par A list containing the step function implementing the CDF.
#'
#' @return A vector of probabilities P(X<x|par).
.cdfKernelEmpirical <- function(x, par){
    par$kernelCDF(x)

} # .cdfKernelEmpirical


#' Internal function to fit a censored (alpha-) stable distribution
#'
#' @description Maximum-likelihood estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param verbose   Provide details on computations.
#' @param file.name   The file name of a PDF file.
#'
#' @return A list with the four stable distribution parameters (S0
#'   representation).
#'
#'   If \code{file.name} is provided, a control plot is generated in a PDF with
#'   a data histogram and the fitted Gaussian. \code{title} is used to give this
#'   plot a main title.
#' @importFrom stabledist pstable dstable
.getAlphaStableParam <- function(d, title, verbose = FALSE, file.name = NULL) {
    if (!is.null(file.name)) {
        grDevices::pdf(file = file.name, width = 4, height = 4,
                       pointsize = 10, useDingbats = FALSE)
        graphics::hist(d, freq=FALSE,
                       main=paste0(title, " / censored stable"),
                       breaks = 30, xlab = "Spearman correlation")
    }

    # ML fit of a censored stable on [-1;1]
    stableLL <- function(par){
        q <- stabledist::pstable(1, alpha=par[1], beta=par[2],
                     gamma=(par[3])**2,delta=par[4]) -
            stabledist::pstable(-1, alpha=par[1], beta=par[2],
                    gamma=(par[3])**2, delta=par[4])
        -sum(stabledist::dstable(d, alpha=par[1], beta=par[2],
                                 gamma=(par[3])**2, delta=par[4], log=TRUE)
             ) + length(d)*log(q)
    }
    par.0 <- c(1.5, 0.5, sqrt(stats::sd(d)), mean(d))
    if (verbose)
        cat(paste0("Starting stable distribution parameter estimation. ",
                   "This can take a few dozens of minutes...\n"))
    res <- stats::optim(par.0, stableLL,
                        control=list(reltol=1e-5, maxit=800,
                                     trace=ifelse(verbose,1,0)))
    if (res$convergence != 0)
        stop("optim() could not fit the censored stable distribution")
    if (verbose)
        cat("Censored stable parameters (alpha, beta, gamma, delta): ",
            paste(res$par, collapse=", "), "\n", sep="")
    alpha <- res$par[1]
    beta <- res$par[2]
    gamma <- (res$par[3])**2
    delta <- res$par[4]
    q <- stabledist::pstable(1, alpha=alpha, beta=beta,
                             gamma=gamma, delta=delta) -
        stabledist::pstable(-1, alpha=alpha, beta=beta,
                            gamma=gamma, delta=delta)
    start <- stabledist::pstable(-1, alpha=alpha, beta=beta,
                                 gamma=gamma, delta=delta)

    # control plot
    if (!is.null(file.name)) {
        x <- seq(-1, 1, by = 0.002)
        graphics::lines(x = x, y = stabledist::dstable(x, alpha=alpha,
                            beta=beta, gamma=gamma, delta=delta)/q,
                        col = "blue", type = "l")
        graphics::legend(x = "topright", lty = 1, legend = "Model",
                         col = "blue", bty = "n", pt.cex = 0.5)
        grDevices::dev.off()
    }
    list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, factor=q,
         start=start, distrib="censored_stable")

}  # .getAlphaStableParam


#' Internal function to compute a censored alpha-) stable CDF
#'
#' @param x   A vector of observed values.
#' @param par A list containing the censored stable model parameters.
#'
#' @return A vector of probabilities P(X<x|par).
#' @importFrom stabledist pstable dstable
.cdfAlphaStable <- function(x, par){

    (stabledist::pstable(x, alpha=par$alpha, beta=par$beta, gamma=par$gamma,
                         delta=par$delta) - par$start) / par$factor

} # .cdfAlphaStable


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
#' @param seed      Seed to reproduce exact same sampling.
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
                             min.pw.size = 5, min.positive = 4,seed=123) {

    pindices <- .buildPermutationIndices(ncounts)
    r.ds <- prepareDataset(ncounts, normalize = FALSE, method = "ALREADY")
    if (foreach::getDoParWorkers() > 1)
        foreach::foreach(k = seq_len(n.rand), .combine = c) %dopar% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices,seed)
            r.LR <- .getCorrelatedLR(r.ds, min.cor = min.cor)
            list(.checkReceptorSignaling(r.ds, r.LR,
                        with.complex = with.complex, max.pw.size = max.pw.size,
                        min.pw.size = min.pw.size, min.positive = min.positive)
            )
        }
    else
        foreach::foreach(k = seq_len(n.rand), .combine = c) %do% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices,seed)
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
#' @param seed      Seed to reproduce exact same sampling.
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
.getEmpiricalNullCorrLR <- function(ncounts, n.rand = 5, min.cor = -1,seed=123) {

    pindices <- .buildPermutationIndices(ncounts)
    r.ds <- prepareDataset(ncounts, normalize = FALSE, method = "ALREADY")

    if (foreach::getDoParWorkers() > 1)
        foreach::foreach(k = seq_len(n.rand), .combine = 'c',
            .packages="BulkSignalR"
            ) %dopar% {
            ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices,seed)
            list(.getCorrelatedLR(r.ds, min.cor = min.cor))
        }
    else
    foreach::foreach(k = seq_len(n.rand), .combine = c) %do% {
        ncounts(r.ds) <- .buildPermutatedCountMatrix(ncounts, pindices,seed)
        list(.getCorrelatedLR(r.ds, min.cor = min.cor))
    }

}  # .getEmpiricalNullCorrLR

