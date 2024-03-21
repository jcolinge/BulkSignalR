library(methods)

#' BulkSignalR cluster comparison-based inference object
#'
#' An S4 class to represent ligand-receptor interactions inferred from
#' a comparison of two clusters of samples. This class inherits from
#' BSRInference.
#'
#' @slot cmp.name  The name of the BSRClusterComp object in a BSRDataModelComp
#' object comp list.
#' @slot tg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot tg.logFC  A list of target gene logFC, one
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
                 tg.pval="list",
                 tg.logFC="list"),
         prototype=list(
           cmp.name="happy",
           LRinter=data.frame(L="A", R="B", pw.id="123",
                              pw.name="one pw", pval=1.0, qval=1.0, L.logFC=2,
                              R.logFC=1.5, LR.pval=0.6, LR.corr=0.5,
                              rank=2, len=50, rank.pval=0.6, rank.corr=0.34,
                              stringsAsFactors=FALSE),
           tg.pval=list(c(0.05,0.1,0.008)),
           tg.logFC=list(c(-1,0,2))
         ))

setValidity("BSRInferenceComp",
            function(object) {
              if (!is.character(object@cmp.name))
                return("cmp.name is not of character type")
              if (length(object@cmp.name) == 0)
                return("cmp.name must have a length > 0")
              if (!is.list(object@tg.pval))
                return("tg.pval is not a list")
              if (!is.list(object@tg.logFC))
                return("tg.logFC is not a list")
              
              TRUE
            }
)

setMethod("show", "BSRInferenceComp",
          function(object) {
            callNextMethod()
            cat("Cluster comparison name:", object@cmp.name, "\n")
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
  
  # recompute P-values
  for (i in 1:nrow(pairs)){
    tg <- t.genes[[i]]
    spvals <- tg.pval[[i]]
    
    # get the LR correlation P-value
    p.lr <- pairs$LR.pval[i]
    
    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    len <- pairs$len[i]
    r <- min(max(1, trunc(rank.p*len)), len)
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
  lg.logFC <- list()
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
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
#' the best ligand-receptor pair that
#' was in this pathway.
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
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
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
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


