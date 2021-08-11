#' Internal function to generate expression matrix permutation indices
#'
#' @param ncounts         A matrix or table of normalized read counts.
#' @param n.bins          Number of bins.
#' @return A list containing a vectors of row indices. The vector at index \code{i} in the list
#' contains the row indices of rows with mean normalized read count in bin \code{i}.
.buildPermutationIndices <- function(ncounts,n.bins=20){

  rm <- rowMeans(ncounts,na.rm=TRUE)
  breaks <- stats::quantile(rm,prob=(0:n.bins)/n.bins)
  breaks[1] <- 0
  lapply(2:length(breaks),function(i) which(rm>breaks[i-1] & rm<=breaks[i]))

} # .buildPermutationIndices


#' Internal function to shuffle permutation indices
#'
#' @param pind      Permutation indices such as returned by \code{.buildPermutationIndices}.
#' @return A list with same structure as \code{pind} but shuffled indices within each bin.
.shufflePermutationIndices <- function(pind){
  lapply(pind,function(x) sample(x,length(x)))
}


#' Internal function to generate a randomized expression matrix
#'
#' @param ncounts    A matrix or table of normalized read counts.
#' @param pind      Permutation indices such as returned by \code{.buildPermutationIndices}.
#' @return \code{ncount} with shuffled row names (gene symbols). Shuffling is performed within rows with comparable average expression.
.buildPermutatedCountMatrix <- function(ncounts,pind){

  symbols <- rownames(ncounts)
  rind <- .shufflePermutationIndices(pind)
  for (i in 1:length(pind))
    symbols[pind[[i]]] <- symbols[rind[[i]]]
  rownames(ncounts) <- symbols
  ncounts

} # .buildPermutatedMatrix


#' Internal function sampling the (empirical) null distribution downstream the receptors
#'
#' Perform receptor downstream analysis with \code{checkReceptorSignaling} based on randomized expression data
#' and ligand-receptor pairs selected from the same randomized data.
#'
#' @param ncounts         A matrix or table of normalized read counts.
#' @param n.rand          The number of repetitions.
#' @param min.cor         The minimum ligand-receptor correlation required.
#' @param max.pw.size     Maximum pathway size to consider from the pathway reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway reference.
#' @param min.positive    Minimum number of target genes to be found in a given pathway.
#' @param with.complex    A logical indicating whether receptor co-complex members should be included in the target genes.
#' @return A list of \code{n.rand} tables such as output by \code{checkReceptorSignaling}. Each table is computed from
#' a randomized expression matrix (randomized \code{ncounts}).
#'
#' A large number of correlations (ligand-receptor and receptor-downstream target genes) is reported in each table. Therefore,
#' \code{n.rand} should be given a modest value to avoid unnecessarily long computations.
#'
#' See \code{\link{checkReceptorSignaling}} for more details about the parameters.
#'
#' This function is made visible to allow advanced users to develop their own statistical models or to
#' estimate their own statistical significance.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds.RT.null <- .getEmpiricalNull(ds$ncounts)
#' }
#' @importFrom foreach %do% %dopar%
.getEmpiricalNull <- function(ncounts,n.rand=5,min.cor=-1,with.complex=TRUE,max.pw.size=200,min.pw.size=5,min.positive=4){

  pindices <- .buildPermutationIndices(ncounts)
  r.ds <- prepareDataset(ncounts,rep("void",ncol(ncounts)),normalize=FALSE)
  if (foreach::getDoParWorkers()>1)
    foreach::foreach(k=1:n.rand,.combine=c) %dopar% {
      r.ds$ncounts <- .buildPermutatedCountMatrix(ncounts,pindices)
      r.LR <- getCorrelatedLR(r.ds,min.cor=min.cor)
      r.LR <- checkReceptorSignaling(r.ds,r.LR,with.complex=with.complex,max.pw.size=max.pw.size,min.pw.size=min.pw.size,min.positive=min.positive)
      list(r.LR$merged.pairs)
    }
  else
    foreach::foreach(k=1:n.rand,.combine=c) %do% {
      r.ds$ncounts <- .buildPermutatedCountMatrix(ncounts,pindices)
      r.LR <- getCorrelatedLR(r.ds,min.cor=min.cor)
      r.LR <- checkReceptorSignaling(r.ds,r.LR,with.complex=with.complex,max.pw.size=max.pw.size,min.pw.size=min.pw.size,min.positive=min.positive)
      list(r.LR$merged.pairs)
    }
  
} # .getEmpiricalNull


#' Internal function sampling of the (empirical) null distribution of ligand-receptor correlations
#'
#' Perform a ligand-receptor correlation analysis based on randomized expression data.
#'
#' @param ncounts         A matrix or table of normalized read counts.
#' @param n.rand          The number of repetitions.
#' @param min.cor         The minimum ligand-receptor correlation required.
#' @return A list of \code{n.rand} tables such as output by \code{getCorrelatedLR}. Each table is computed from
#' a randomized expression matrix (randomized \code{ncounts}).
#'
#' A large number of correlations is reported in each table. Therefore,
#' \code{n.rand} should be given a modest value to avoid unnecessarily long computations.
#'
#' See \code{\link{getCorrelatedLR}} for more details about the parameters.
#'
#' This function is made visible to allow advanced users to develop their own statistical models or to
#' estimate their own statistical significance.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds.LR.null <- .getEmpiricalNullCorrLR(ds$ncounts,sample.set="normal vs. tumor")
#' }
#' @importFrom foreach %do% %dopar%
.getEmpiricalNullCorrLR <- function(ncounts,n.rand=5,min.cor=-1){

  pindices <- .buildPermutationIndices(ncounts)
  r.ds <- prepareDataset(ncounts,rep("void",ncol(ncounts)),normalize=FALSE)
  foreach::foreach(k=1:n.rand,.combine=c) %do% {
    r.ds$ncounts <- .buildPermutatedCountMatrix(ncounts,pindices)
    r.LR <- getCorrelatedLR(r.ds,min.cor=min.cor)
    list(r.LR$putative.pairs)
  }
  
} # .getEmpiricalNullCorrLR



# ===========================================================================================
# Alternative distribution methods -----------------------------------------------------------------


#' Internal function to score pathways
#'
#' Perform hypergeometric tests and collect intersection sizes.
#'
#' @param genes     The set of differentialls expressed genes.
#' @param sp2pw     A pathway reference table.
#' @param fdr.proc  The procedure for adjusting P-values accoring to \code{\link[multtest]{mt.rawp2adjp}}.
#' @param col.id    Column index or name in \code{pw} for the gene symbols.
#' @param col.term  Column index or name in \code{pw} for the pathway IDs.
#' @param col.descr Column index or name in \code{pw} for the pathway names.
#' @param max.pw.size     Maximum pathway size to consider from the pathway reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway reference.
#' @param min.inter       Minimum number of regulated genes in a pathway for this pathway to be reported.
#' @return A table of pathways from the reference table \code{sp2pw} with different statistics and P-values.
.pwHyperScores <- function(genes,sp2pw,fdr.proc=c("BH","Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BY","ABH","TSBH"),
                           col.id,col.term,col.descr,max.pw.size=200,min.pw.size=5,min.inter=2){

  fdr.proc <- match.arg(fdr.proc)
  terms <- unique(sp2pw[[col.term]])
  n.genes <- length(unique(sp2pw[[col.id]]))
  res <- NULL
  for (term in terms){
    in.pw <- unique(sp2pw[sp2pw[[col.term]]==term,col.id])
    if ((length(in.pw)>=min.pw.size) && (length(in.pw)<=max.pw.size)){
      prots <- intersect(genes,in.pw)
      n <- length(prots)
      if (n >= min.inter){
        pval <- 1-stats::phyper(q=n-1,m=length(in.pw),n=n.genes-length(in.pw),k=length(genes))
        res <- rbind(res,data.frame(term=term,descr=sp2pw[sp2pw[[col.term]]==term,col.descr][1],intersect.size=n,
                                    n.pw=length(in.pw),cov=n/length(in.pw),pval=pval,prots=paste0(prots,collapse=","),stringsAsFactors=FALSE))
      }
    }
  }
  if (!is.null(res))
    if (dim(res)[1] > 1){
      rawp <- res[["pval"]]
      adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
      qval <- adj$adjp[order(adj$index),fdr.proc]
      res <- cbind(res,data.frame(qval=qval))
    }
  res

} # .pwHyperScore


#' Internal function to search for differentially regulated pathways
#'
#' Perform differential gene analysis with edgeR and score pathways.
#'
#' @param ncounts   A matrix or table of normalized read counts.
#' @param types     A character vector defining the sample types.
#' @param induced   A character vector defining the sample types that correspond to the "induced" condition
#' @param normal   A character vector defining the sample types that correspond to the "normal" condition
#' @param max.pw.size     Maximum pathway size to consider from the pathway reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway reference.
#' @param max.qval        Maximum Q-value imposed to genes.
#' @param min.FC          Minimum log2 fold-change (in absolute value) imposed to genes.
#' @param min.count       Minimum median expression value across all the samples imposed to genes.
#' @return A table of pathways as returned by \code{\link{.pwHyperScores}}.
#'
#' The 2-condition induced versus normal is performed using edgeR \code{\link[edgeR]{glmFit}}.
#'
#' No normalization is performed with edgeR since the matrix is supposed to be normalized already. Note that it is possible
#' to provide \code{\link{prepareDataset}} with a matrix that would have been pre-normalized using
#' \code{\link[edgeR]{calcNormFactors}} if so desired. (Remember that \code{\link[edgeR]{calcNormFactors}} does not normalize the matrix, it only
#' returns the normalization factors.)
.getDifferentialPathways <- function(ncounts,types,normal,induced,max.pw.size=200,min.pw.size=5,max.qval=0.01,min.FC=log2(1.5),min.count=100){

  if (!requireNamespace("edgeR",quietly=TRUE))
    stop("Package \"edgeR\" needed for this function to work. Please install it.")

  dge <- edgeR::DGEList(ncounts[,types%in%normal | types%in%induced],genes=rownames(ncounts))
  rownames(dge$counts) <- rownames(ncounts)
  conditions <- types[types%in%normal | types%in%induced]
  conditions[conditions%in%normal] <- "normal"
  conditions[conditions%in%induced] <- "induced"
  cl <- factor(conditions)
  design <- stats::model.matrix(~0+cl)
  colnames(design) <- gsub("^cl","",colnames(design))
  cm <- limma::makeContrasts(induced-normal,levels=design)
  y <- edgeR::estimateDisp(dge,design,robust=TRUE)
  fit.y <- edgeR::glmFit(y,design)
  comparison <- "induced - normal"
  lrt <- edgeR::glmLRT(fit.y,contrast=cm[,comparison])
  sel.r <- edgeR::topTags(lrt,adjust.method="BH",n=nrow(y$counts))
  sel.r <- sel.r[sel.r$table$FDR<max.qval & abs(sel.r$table$logFC)>min.FC & apply(ncounts[rownames(sel.r),],1,stats::median)>min.count,]
  sel.genes <- rownames(sel.r)

  .pwHyperScores(sel.genes,reactome,col.id="Gene name",col.term="Reactome ID",col.descr="Reactome name",max.pw.size=max.pw.size,min.pw.size=min.pw.size)

} # .getDifferentialPathways


# ===============================================================================
# Parameter learning ------------------------------------------------------------


#' Internal function to fit a Gaussian distribution
#'
#' Maximum-likelihood estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param file.name   The file name of a PDF file.
#' @return A list with the mean (\code{mu}) and standard deviation (\code{sigma}) estimates.
#'
#' If \code{file.name} is provided, a control plot is generated in a PDF with a data histogram and
#' the fitted Gaussian. \code{title} is used to give this plot a main title.
.getGaussianParam <- function(d,title,file.name=NULL){
  if (!is.null(file.name)){
    grDevices::pdf(file=file.name,width=4,height=4,pointsize=10,useDingbats=FALSE)
    graphics::hist(d,freq=FALSE,main=title,xlab="Spearman correlation")
  }
  mu <- mean(d)
  sigma <- stats::sd(d)
  if (!is.null(file.name)){
    x <- seq(-1,1,by=0.01)
    graphics::lines(x=x,y=stats::dnorm(x,mu,sigma),col="blue",type="l")
    graphics::legend(x="topright",lty=1,legend="Normal",col="blue",bty="n",pt.cex=0.5)
    grDevices::dev.off()
  }
  list(mu=mu,sigma=sigma)

} # .getGaussianParam


#' Internal function to fit a Gamma distribution
#'
#' Method of moments estimators are used.
#'
#' @param d   A vector of values to fit.
#' @param title     A plot title.
#' @param file.name   The file name of a PDF file.
#' @return A list with the shape (\code{k}) and scale (\code{theta}) estimates.
#'
#' If \code{file.name} is provided, a control plot is generated in a PDF with a data histogram and
#' the fitted Gamma (plus a Gaussian for comparison). \code{title} is used to give this plot a main title.
.getGammaParam <- function(d,title,file.name=NULL){
  z <- 1+d
  if (!is.null(file.name)){
    grDevices::pdf(file=file.name,width=4,height=4,pointsize=10,useDingbats=FALSE)
    graphics::hist(z,freq=FALSE,main=title,xlab="Spearman correlation")
    x <- seq(0,2,by=0.01)
    mu <- mean(d)
    sigma <- stats::sd(d)
    graphics::lines(x=x,y=stats::dnorm(x-1,mu,sigma),col="blue",type="l")
  }
  mu <- mean(z)
  s2 <- stats::var(z)
  k <- mu**2 / s2
  theta <- s2/mu
  if (!is.null(file.name)){
    graphics::lines(x=x,y=stats::dgamma(x,shape=k,scale=theta),col="red",type="l")
    graphics::legend(x="topright",lty=1,legend=c("Normal","Gamma"),col=c("blue","red"),bty="n",pt.cex=0.5)
    grDevices::dev.off()
  }
  list(k=k,theta=theta)

} # .getGammaParam


#' Training of BulkSignalR model parameters
#'
#' Unique entry point for training the parameters behind BulkSignalR statistical models.
#'
#' @param ds        A BulkSignalR data set, i.e., a list with the normalized read counts and sample types.
#' @param plot.folder   A folder name for generating control plots.
#' @param verbose       A logical activating progress messages for the user.
#' @param n.rand.LR     The number of random expression matrices to use for learning the ligand-receptor correlation distribution.
#' @param n.rand.RT     The number of random expression matrices to use for learning the receptor-target genes correlation distribution.
#' @param with.complex    A logical indicating whether receptor co-complex members should be included in the target genes.
#' @param max.pw.size     Maximum pathway size to consider from the pathway reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway reference.
#' @param min.positive    Minimum number of target genes to be found in a given pathway.
#' @param quick           A logical indicating whether approximate parameters for the receptor-target correlations should be used.
#' @return The BulkSignalR data set \code{ds} extended with a slot \code{ds$param} that contains all the null model parameters
#' as well as parameters that were used to specify how they were learnt. Provided \code{induced} and \code{normal} are set, the
#' alternative model parameters are also added into \code{ds$param}. \code{ds$param} is a multi-level named list.
#'
#' In a reference pathway, i.e., a Reactome pathway or the genes of a GOBP term, the target genes are the
#' genes coding for proteins forming a complex with the receptor and the genes in the pathway downstream the receptor,
#' which are given as regulated by the pathway. If \code{with.complex} is set to \code{FALSE}, then only the
#' regulated genes are considered. Participation to a complex and being regulated as well as the pathway directed topologies
#' are defined by Reactome and KEGG pathways as provided by PathwayCommons. Those parameters should be identical to the
#' values intended when searching for ligand-receptor pairs with \code{\link{getCorrelatedLR}} and
#' \code{\link{checkReceptorSignaling}}. Although the statistical distributions are rather robust it is not advisable
#' to use different parameters that could nonetheless introduce biases, but for saving compute time and exploring.
#'
#' The maximum pathway size is used to limit the redundancy inherent to GOBP and Reactome. The minimum pathway size is
#' used to avoid overspecific, noninformative results.
#'
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#' }
learnParameters <- function(ds,plot.folder=NULL,verbose=FALSE,n.rand.LR=5,n.rand.RT=2,
                            with.complex=TRUE,max.pw.size=200,min.pw.size=5,min.positive=4,quick=TRUE){

  n.rand.LR <- trunc(n.rand.LR)
  if (n.rand.LR < 1)
    stop("Parameter n.rand.LR muste be an integer > 0")
  n.rand.RT <- trunc(n.rand.RT)
  if (n.rand.RT < 1)
    stop("Parameter n.rand.RT muste be an integer > 0")

  param <- list()
  param$plot.folder <- plot.folder
  param$with.complex <- with.complex
  param$max.pw.size <- max.pw.size
  param$min.pw.size <- min.pw.size
  param$min.positive <- min.positive
  param$quick <- quick

  # LR correlation null ----------------

  induced.samples <- 1:ncol(ds$ncounts)
  if (verbose)
    cat("Learning ligand-receptor correlation null distribution...\n")
  ds.LR.null <- .getEmpiricalNullCorrLR(ds$ncounts[,induced.samples],n.rand=n.rand.LR)
  rc <- ds.LR.null[[1]]$corr
  if (length(ds.LR.null) > 1)
    for (i in 2:length(ds.LR.null))
      rc <- c(rc,ds.LR.null[[2]]$corr)
  param$LR.0$n.rand <- n.rand.LR
  param$LR.0$n <- length(rc)

  # Gaussian model
  if (!is.null(plot.folder))
    file.name <- paste0(plot.folder,"LR-null-hist-gaussian.pdf")
  else
    file.name <- NULL
  gp <- .getGaussianParam(rc,"LR correlation (null)",file.name)
  param$LR.0$norm$mu <- gp$mu
  param$LR.0$norm$sigma <- gp$sigma

  # Gamma model (and comparison with Gaussian)
  if (!is.null(plot.folder))
    file.name <- paste0(plot.folder,"LR-null-hist-gamma.pdf")
  else
    file.name <- NULL
  gap <- .getGammaParam(rc,"LR correlation (null)",file.name)
  param$LR.0$gamma$k <- gap$k
  param$LR.0$gamma$theta <- gap$theta

  # RT correlation null ------------------------------------

  if (quick){
    # RT correlations are assumed to be equal to LR correlations
    if (verbose)
      cat("Quick learning, receptor-target correlation null distribution assumed to be equal to ligand-receptor...\n")
    param$RT.0$n <- param$LR.0$n
    param$RT.0$norm$mu <- param$LR.0$norm$mu
    param$RT.0$norm$sigma <- param$LR.0$norm$sigma
    param$RT.0$gamma$k <- param$LR.0$gamma$k
    param$RT.0$gamma$theta <- param$LR.0$gamma$theta
  }
  else{
    # RT correlations are actually learnt
    if (verbose)
      cat("Learning receptor-target correlation null distribution...\n")
    ds.RT.null <- .getEmpiricalNull(ds$ncounts[,induced.samples],n.rand=n.rand.RT,with.complex=with.complex,max.pw.size=max.pw.size,min.pw.size=min.pw.size,min.positive=min.positive)
    param$RT.0$n.rand <- n.rand.RT
    t <- ds.RT.null[[1]]
    if (length(ds.RT.null) > 1)
      for (i in 2:length(ds.RT.null))
        t <- rbind(t,ds.RT.null[[2]])
    above <- unlist(strsplit(t$all.corr,split="\\|"))
    r.corrs <- NULL
    for (i in 1:length(above)){
      corr <- as.numeric(strsplit(above[i],split=";")[[1]])
      r.corrs <- c(r.corrs,corr)
    }
    param$RT.0$n <- length(r.corrs)
    
    # Gaussian model
    if (!is.null(plot.folder))
      file.name <- paste0(plot.folder,"RT-null-hist-gaussian.pdf")
    else
      file.name <- NULL
    gp <- .getGaussianParam(r.corrs,"RT correlation (null)",file.name)
    param$RT.0$norm$mu <- gp$mu
    param$RT.0$norm$sigma <- gp$sigma
    
    # Gamma model
    if (!is.null(plot.folder))
      file.name <- paste0(plot.folder,"RT-null-hist-gamma.pdf")
    else
      file.name <- NULL
    gap <- .getGammaParam(r.corrs,"RT correlation (null)",file.name)
    param$RT.0$gamma$k <- gap$k
    param$RT.0$gamma$theta <- gap$theta
  }

  # Results ----------------------------------
  if (verbose)
    cat("Learning of statistical model parameters completed\n")
  ds$param <- param
  ds

} # learnParameters

