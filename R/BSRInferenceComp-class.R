library(methods)

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
#' print('BSRInferenceComp class')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1,"random.example")
#' 
setClass("BSRInferenceComp",
         contains="BSRInference",
         slots=c(cmp.name="character",
                 src.cmp.name="character",
                 tg.pval="list",
                 tg.logFC="list",
                 tg.expr="list"),
         prototype=list(
           cmp.name="happy",
           src.cmp.name="",
           LRinter=data.frame(L="A", R="B", pw.id="123",
                              pw.name="one pw", pval=1.0, qval=1.0, L.logFC=2,
                              R.logFC=1.5, LR.pval=0.6, LR.corr=0.5,
                              rank=2, len=50, rank.pval=0.6, rank.corr=0.34,
                              LR.score=0.6, L.expr=2, R.expr=3,
                              stringsAsFactors=FALSE),
           tg.pval=list(c(0.05,0.1,0.008)),
           tg.logFC=list(c(-1,0,2)),
           tg.expr=list(c(1,2,3))
         ))

setValidity("BSRInferenceComp",
            function(object) {
              if (!is.character(object@cmp.name))
                return("cmp.name is not of character type")
              if (!is.character(object@src.cmp.name))
                return("src.cmp.name is not of character type")
              if (length(object@cmp.name) == 0)
                return("cmp.name must have a length > 0")
              if (!is.list(object@tg.pval))
                return("tg.pval is not a list")
              if (!is.list(object@tg.logFC))
                return("tg.logFC is not a list")
              if (!is.list(object@tg.expr))
                return("tg.expr is not a list")
              
              TRUE
            }
)

setMethod("show", "BSRInferenceComp",
          function(object) {
            callNextMethod()
            cat("Cluster comparison name:", object@cmp.name, "\n")
            cat("Source cluster comparison name:", object@src.cmp.name, "\n")
          }
)


# Accessors & setters ========================================================

if (!isGeneric("cmpName")) {
  if (is.function("cmpName"))
    fun <- cmpName
  else
    fun <- function(x) standardGeneric("cmpName")
  setGeneric("cmpName", fun)
}
#' Comparison name accessor
#'
#' @name cmpName
#' @aliases cmpName,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("cmpName", "BSRInferenceComp", function(x) x@cmp.name)

if (!isGeneric("cmpName<-")) {
  if (is.function("cmpName<-"))
    fun <- `cmpName<-`
  else
    fun <- function(x,value) standardGeneric("cmpName<-")
  setGeneric("cmpName<-", fun)
}
#' Comparison name setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("cmpName<-", "BSRInferenceComp", function(x, value){
  x@cmp.name <- value
  methods::validObject(x)
  x
})


if (!isGeneric("srcCmpName")) {
  if (is.function("srcCmpName"))
    fun <- srcCmpName
  else
    fun <- function(x) standardGeneric("srcCmpName")
  setGeneric("srcCmpName", fun)
}
#' Source comparison name accessor
#'
#' @name srcCmpName
#' @aliases srcCmpName,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("srcCmpName", "BSRInferenceComp", function(x) x@src.cmp.name)

if (!isGeneric("srcCmpName<-")) {
  if (is.function("srcCmpName<-"))
    fun <- `srcCmpName<-`
  else
    fun <- function(x,value) standardGeneric("srcCmpName<-")
  setGeneric("srcCmpName<-", fun)
}
#' Source comparison name setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("srcCmpName<-", "BSRInferenceComp", function(x, value){
  x@src.cmp.name <- value
  methods::validObject(x)
  x
})


if (!isGeneric("tgPval")) {
  if (is.function("tgPval"))
    fun <- tgPval
  else
    fun <- function(x) standardGeneric("tgPval")
  setGeneric("tgPval", fun)
}
#' Target gene P-values accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("tgPval", "BSRInferenceComp", function(x) x@tg.pval)

if (!isGeneric("tgPval<-")) {
  if (is.function("tgPval<-"))
    fun <- `tgPval<-`
  else
    fun <- function(x,value) standardGeneric("tgPval<-")
  setGeneric("tgPval<-", fun)
}
#' Target gene P-values setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("tgPval<-", "BSRInferenceComp", function(x, value){
  x@tg.pval <- value
  methods::validObject(x)
  x
})


if (!isGeneric("tgLogFC")) {
  if (is.function("tgLogFC"))
    fun <- tgLogFC
  else
    fun <- function(x) standardGeneric("tgLogFC")
  setGeneric("tgLogFC", fun)
}
#' Target gene logFC accessor
#'
#' @name tgLogFC
#' @aliases tgLogFC,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("tgLogFC", "BSRInferenceComp", function(x) x@tg.logFC)

if (!isGeneric("tgLogFC<-")) {
  if (is.function("tgLogFC<-"))
    fun <- `tgLogFC<-`
  else
    fun <- function(x,value) standardGeneric("tgLogFC<-")
  setGeneric("tgLogFC<-", fun)
}
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("tgLogFC<-", "BSRInferenceComp", function(x, value){
  x@tg.logFC <- value
  methods::validObject(x)
  x
})


if (!isGeneric("tgExpr")) {
  if (is.function("tgExpr"))
    fun <- tgExpr
  else
    fun <- function(x) standardGeneric("tgExpr")
  setGeneric("tgExpr", fun)
}
#' Target gene expression accessor
#'
#' @name tgExpr
#' @aliases tgExpr,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("tgExpr", "BSRInferenceComp", function(x) x@tg.expr)

if (!isGeneric("tgExpr<-")) {
  if (is.function("tgExpr<-"))
    fun <- `tgExpr<-`
  else
    fun <- function(x,value) standardGeneric("tgExpr<-")
  setGeneric("tgExpr<-", fun)
}
#' Target gene expression setter (internal use only)
#' @param x BSRInferenceComp object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("tgExpr<-", "BSRInferenceComp", function(x, value){
  x@tg.expr <- value
  methods::validObject(x)
  x
})


# simplified table views =======================================================

if (!isGeneric("LRinterShort")) {
  if (is.function("LRinterShort"))
    fun <- LRinterShort
  else
    fun <- function(x) standardGeneric("LRinterShort")
  setGeneric("LRinterShort", fun)
}
#' Simplified LRinter accessor reporting the essential columns
#'
#' @name LRinterShort
#' @aliases LRinterShort,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("LRinterShort", "BSRInferenceComp",
          function(x) x@LRinter[,c("L", "R", "pw.id", "pw.name",
                                   "qval", "L.logFC", "R.logFC", "len")]
)


if (!isGeneric("LRinterScore")) {
  if (is.function("LRinterScore"))
    fun <- LRinterScore
  else
    fun <- function(x) standardGeneric("LRinterScore")
  setGeneric("LRinterScore", fun)
}
#' Simplified LRinter accessor with focus on the LR-score
#'
#' @name LRinterScore
#' @aliases LRinterScore,BSRInferenceComp-method
#' @param x BSRInferenceComp object
#' @export
setMethod("LRinterScore", "BSRInferenceComp",
          function(x) x@LRinter[,c("L", "R", "pw.id", "pw.name",
                                   "LR.score", "L.expr", "R.expr", "len")]
)


# Rescoring & updating ================================================================

if (!isGeneric("rescoreInference")) {
  if (is.function("rescoreInference"))
    fun <- rescoreInference
  else
    fun <- function(obj, ...) standardGeneric("rescoreInference")
  setGeneric("rescoreInference", fun)
}
#' Inference re-scoring
#'
#' A method to re-score an existing BSRInferenceComp object
#' (P- and Q-value estimations).
#'
#' @name rescoreInference
#' @aliases rescoreInference,BSRInferenceComp-method
#'
#' @param obj BSRInferecenceComp object.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @details A BSRInferenceComp object should be created by calling 
#' \code{"\link[=BSRClusterComp-class]{initialInference}"}
#'
#' @return A BSRInferenceComp object.
#'
#' @export
#' @examples
#' print('rescoreInference')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # rescore
#' bsrinf.less <- rescoreInference(bsrinf, param=param(bsrdm.comp), rank.p=0.75)
#'
setMethod("rescoreInference", "BSRInferenceComp", function(obj, param, rank.p=0.55,
            fdr.proc=c("BH", "Bonferroni", "Holm",
                       "Hochberg", "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {
  
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  
  # extract the necessary data from the BSRInferenceComp object
  pairs <- LRinter(obj)
  t.genes <- tGenes(obj)
  tg.pval <- tgPval(obj)
  tg.corr <- tgCorr(obj)

  # recompute P-values, LR- and LRT-scores
  for (i in 1:nrow(pairs)){
    tg <- t.genes[[i]]
    spvals <- tg.pval[[i]]
    
    # get the LR correlation P-value
    p.lr <- pairs$LR.pval[i]
    
    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    len <- pairs$len[i]
    r <- min(max(1, trunc(rank.p*len)), len)
    pairs$rank[i] <- r
    rank.pval <- spvals[r]
    # r-1 P-values are > rank.pval, prob to have r-1 or less
    # P-values > rank.pval is given by a binomial with success rate
    # equal to the probability to get a P-value > rank.pval, i.e.,
    # 1-rank.pval. If rank.pval is low (i.e., highly significant),
    # it becomes difficult to get as little as r-1 P-values > rank.pval by chance!
    p.rt <- stats::pbinom(r-1, len, 1-rank.pval) # cdf is punif here!
    pairs$pval[i] <- p.lr*p.rt
    pairs$rank.pval[i] <- rank.pval
    pairs$rank.corr[i] <- tg.corr[[i]][r]
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


if (!isGeneric("updateInference")) {
  if (is.function("updateInference"))
    fun <- updateInference
  else
    fun <- function(obj, ...) standardGeneric("updateInference")
  setGeneric("updateInference", fun)
}
#' Inference updating
#'
#' A method to update the data underlying statistical significance estimations
#' prior to rescoring for an existing BSRInferenceComp object
#' (P- and Q-value estimations as well as LR-score).
#'
#' @name updateInference
#' @aliases updateInference,BSRInferenceComp-method
#'
#' @param obj BSRInferecenceComp object.
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
#' \code{"\link[=BSRClusterComp-class]{initialInference}"}
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
#' Note that correlations are set to 1 to avoid
#' lengthy computations with scRNA-seq data and multiple cell
#' populations.
#' 
#' The main function of this method is to support our SingleCellSignalR v2
#' package.
#'
#' @export
#' @examples
#' print('updateInference')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # update
#' stats$pval <- stats$pval/100
#' stats$logFC <- stats$logFC + 0.5
#' bsrcc.2 <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrinf.updated <- updateInference(bsrinf, bsrcc.2)
#'
setMethod("updateInference", "BSRInferenceComp", function(obj, bsrcc, ncounts,
                                src.bsrcc=NULL, rank.p=0.55,
                                max.pval=0.01, min.logFC=1, min.LR.score=0,
                                neg.receptors=FALSE,
                                pos.targets=FALSE, neg.targets=FALSE,
                                min.t.logFC=0.5, min.positive=2,
                                fdr.proc=c("BH","Bonferroni","Holm","Hochberg",
                                           "SidakSS","SidakSD","BY","ABH","TSBH")){

  if (!is(obj, "BSRInferenceComp"))
    stop("obj must be a BSRInferenceComp object")
  if (!is(bsrcc, "BSRClusterComp"))
    stop("bsrcc must be a BSRClusterComp object")
  if (!is.null(src.bsrcc) && !is(src.bsrcc, "BSRClusterComp"))
    stop("src.bsrcc must be a BSRClusterComp object")
    
  # get BSRInferenceComp object details
  inter <- LRinter(obj)
  L <- ligands(obj)
  R <- receptors(obj)
  t <- tGenes(obj)
  c <- tgCorr(obj)
  lfc <- tgLogFC(obj)
  p <- tgPval(obj)
  e <- tgExpr(obj)
  par <- infParam(obj)
  mu <- par$mu
  logTransf <- par$log.transformed.data
  R.stats <- stats(bsrcc)
  if (!is.null(src.bsrcc))
    L.stats <- stats(src.bsrcc) # paracrine
  else
    L.stats <- R.stats # autocrine
  
  # update parameters
  par$colA <- colA(bsrcc)
  par$colB <- colB(bsrcc)
  if (is.null(src.bsrcc))
    par$inference.type <- "autocrine"
  else{
    par$inference.type <- "paracrine"
    par$src.colA <- colA(src.bsrcc)
    par$src.colB <- colB(src.bsrcc)
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

  # if (is.null(src.bsrcc))
  #   corlr <- stats::cor(t(ncounts[, c(colA(bsrcc),colB(bsrcc))]), method = "spearman")
  # else
  #   corlr <- stats::cor(t(ncounts[, c(colA(bsrcc),colA(src.bsrcc))]), method = "spearman")
  # for (i in seq_len(nrow(inter)))
  #   inter$LR.corr[i] <- corlr[inter$L[i], inter$R[i]]
  inter$LR.corr <- 1
  
  # LR-score
  if (logTransf)
    sq <- sqrt(inter$L.expr * inter$R.expr)
  else
    sq <- sqrt(log1p(inter$L.expr)/log(2) * log1p(inter$R.expr)/log(2))
  inter$LR.score <- sq/(mu+sq)

  # select on L & R as well as LR-scores
  good <- L.stats[inter$L,"pval"] <= max.pval & inter$L.logFC >= min.logFC &
    R.stats[inter$R,"pval"] <= max.pval
  if (neg.receptors)
    good <- good & abs(inter$R.logFC) >= min.logFC
  else
    good <- good & inter$R.logFC >= min.logFC
  good <- good & inter$LR.score >= min.LR.score
  if (sum(good) == 0)
    stop("No selection")    
  inter <- inter[good,]
  L <- L[good]
  R <- R[good]
  t <- t[good]
  c <- c[good]
  lfc <- lfc[good]
  p <- p[good]
  e <- e[good]

  # assign correct logFC, P-values, correlations, and expression to the targets
  # and select the targets
  keep <- NULL
  for (i in seq_len(nrow(inter))){
    genes <- t[[i]]
    logfc <- R.stats[genes, "logFC"]
    if (pos.targets || neg.targets){
      if (pos.targets)
        genes <- genes[logfc >= min.t.logFC]
      else
        genes <- genes[logfc <= -min.t.logFC]
    }
    if (length(genes) < min.positive)
      keep <- c(keep, FALSE)
    else{
      # if all conditions are met, list all target genes with
      # their regulation P-values in a data frame
      # row. Target genes are sorted wrt P-values in decreasing
      # order to keep the compatibility with correlation analysis,
      # where the most significant values are at the end.
      pv <- R.stats[genes, "pval"]
      o <- order(pv, decreasing=TRUE)
      pv <- pv[o]
      logfc <- R.stats[genes, "logFC"]
      logfc <- logfc[o]
      expr <- R.stats[genes, "expr"]
      expr <- expr[o]
      genes <- genes[o]
      co <- rep(1, length(genes)) #corlr[inter[i, "R"], genes]
      t[[i]] <- genes
      c[[i]] <- co
      lfc[[i]] <- logfc
      p[[i]] <- pv
      e[[i]] <- expr
      len <- length(genes)
      inter[i, "len"] <- len
      rank <- min(max(1, trunc(rank.p*len)), len)
      inter[i, "rank"] <- rank
      # rank.corr and rank.pval are updated by calling rescoreInference() below
      keep <- c(keep, TRUE)
    }
  }
  if (sum(keep) == 0)
    stop("No selection")    
  inter <- inter[keep,]
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
  tGenes(obj) <- t
  tgCorr(obj) <- c
  tgPval(obj) <- p
  tgLogFC(obj) <- lfc
  tgExpr(obj) <- e
  infParam(obj) <- par

  # rescore and return
  rescoreInference(obj, par, rank.p=rank.p, fdr.proc=fdr.proc)

}) # updateInferenceComp




# Reduction and pathway stat methods ===========================================


if (!isGeneric("reduceToBestPathway")) {
  if (is.function("reduceToBestPathway"))
    fun <- reduceToBestPathway
  else
    fun <- function(obj, ...) standardGeneric("reduceToBestPathway")
  setGeneric("reduceToBestPathway", fun)
}
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
#' print('reduceToBestPathway')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redBP  <- reduceToBestPathway(bsrinf)
#'
#' @importFrom rlang .data
setMethod("reduceToBestPathway", "BSRInferenceComp", function(obj) {
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual data representation
  
  # get best p-value pathway per LR pair
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  tg.expr <- list()
  
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
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    tg.expr <- c(tg.expr, obj@tg.expr[j])
    LRinter <- rbind(LRinter, pairs[j,])
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@tg.expr <- tg.expr
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
#' print('reduceToReceptor')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redR  <- reduceToReceptor(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToReceptor", "BSRInferenceComp", function(obj){
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual data representation
  
  if (infParam(obj)$ligand.reduced)
    stop("Already reduced to receptor") # because ligands were reduced
  
  # pool the ligands
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  tg.expr <- list()
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
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    tg.expr <- c(tg.expr, obj@tg.expr[j])
    to.add <- pairs[j,]
    to.add[1, "L"] <- paste0("{", paste(unique(lig$L), collapse=";"), "}")
    to.add[1, "LR.score"] <- max(lig$LR.score)
    to.add[1, "L.expr"] <- max(lig$L.expr)
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@tg.expr <- tg.expr
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
#' print('reduceToLigand')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redL  <- reduceToLigand(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToLigand", "BSRInferenceComp", function(obj){
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual representation
  
  if (infParam(obj)$receptor.reduced)
    stop("Already reduced to ligand") # because receptors were reduced
  
  # pool the receptors
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  tg.expr <- list()
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
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    tg.expr <- c(tg.expr, obj@tg.expr[j])
    to.add <- pairs[j,]
    to.add[1, "R"] <- paste0("{", paste(unique(rec$R), collapse=";"), "}")
    to.add[1, "LR.score"] <- max(lig$LR.score)
    to.add[1, "R.expr"] <- max(lig$R.expr)
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@tg.expr <- tg.expr
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
#' print('reduceToPathway')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' @importFrom rlang .data
setMethod("reduceToPathway", "BSRInferenceComp", function(obj){
  
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
  tg.pval <- list()
  tg.logFC <- list()
  tg.expr <- list()
  LRinter <- NULL
  pairs <- obj@LRinter
  for (p in unique(pairs$pw.id)){
    j <- which(pairs$pw.id==p)[1]
    ligands <- c(ligands, list(unique(pairs$L[pairs$pw.id==p])))
    receptors <- c(receptors, list(unique(pairs$R[pairs$pw.id==p])))
    t.genes <- c(t.genes, obj@t.genes[j])
    tg.corr <- c(tg.corr, obj@tg.corr[j])
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    tg.expr <- c(tg.expr, obj@tg.expr[j])
    to.add <- pairs[j,]
    to.add[1, "L"] <- paste0("{", paste(unique(pairs$L[pairs$pw.id==p]),
                                        collapse=";"), "}")
    to.add[1, "R"] <- paste0("{", paste(unique(pairs$R[pairs$pw.id==p]),
                                        collapse=";"), "}")
    to.add[1, "LR.score"] <- max(pairs$LR.score[pairs$pw.id==p])
    to.add[1, "R.expr"] <- max(pairs$R.expr[pairs$pw.id==p])
    to.add[1, "L.expr"] <- max(pairs$L.expr[pairs$pw.id==p])
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@tg.expr <- tg.expr
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
#' \code{"\link[=BSRInferenceComp-class]{scoreLRGeneSignatures}"}
#'
#' @name getLRGeneSignatures
#' @aliases getLRGeneSignatures,BSRInferenceComp-method
#'
#' @param obj    BSRInferenceComp object.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param with.pw.id    A logical indicating whether the ID of a pathway
#' should be concatenated to its name.
#' @return A BSRSignatureComp object containing a gene signature for each triple
#' ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is
#' comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}.
#' @export
#' @examples
#' print('getLRGeneSignatures')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
#'                     expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reductions
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' bsrsig.redP <- getLRGeneSignatures(bsrinf.redP,qval.thres=0.001)
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom methods new
setMethod("getLRGeneSignatures", "BSRInferenceComp", function(obj,
                                                          pval.thres=NULL, qval.thres=NULL, with.pw.id=FALSE){
  
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
  pairs <- pairs[selected,]
  ligands <- ligands(obj)[selected]
  receptors <- receptors(obj)[selected]
  if (with.pw.id)
    pathways <- paste0(pairs$pw.id, ": ", pairs$pw.name)
  else
    pathways <- pairs$pw.name
  t.genes <- tGenes(obj)[selected]
  t.corrs <- tgCorr(obj)[selected]
  t.pvals <- tgPval(obj)[selected]
  t.logFCs <- tgLogFC(obj)[selected]
  
  for (i in seq_len(nrow(pairs))){
    tg <- t.genes[[i]]
    t.genes[[i]] <- tg[pairs$rank[i]:length(tg)]
    tc <- t.corrs[[i]]
    t.corrs[[i]] <- tc[pairs$rank[i]:length(tc)]
    tp <- t.pvals[[i]]
    t.pvals[[i]] <- tp[pairs$rank[i]:length(tp)]
    tl <- t.logFCs[[i]]
    t.logFCs[[i]] <- tl[pairs$rank[i]:length(tl)]
  }
  
  new("BSRSignatureComp", pathways=pathways, ligands=ligands,
      receptors=receptors, t.genes=t.genes, tg.corr=t.corrs,
      tg.pval=t.pvals, tg.logFC=t.logFCs,
      cmp.name=cmpName(obj))
  
}) # getLRGeneSignatures


