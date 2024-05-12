library(methods)

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
#'                    expr=runif(n, 0, 10))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp, max.pval=1,"random.example")
#' 
setClass("BSRDataModelComp",
         contains=c("BSRDataModel"),
         slots=c(comp="list",
                 mu="numeric"),
         prototype=list(
           initial.organism="hsapiens",
           initial.orthologs=list("A","B","C"),
           ncounts=matrix(1.0,nrow=2,ncol=1,dimnames=list(c("A","B"),"C")),
           log.transformed=FALSE,
           normalization="UQ",
           param=list(spatial.smooth=FALSE),
           comp=list(),
           mu=1.0
         ))

setValidity("BSRDataModelComp",
            function(object) {
              if (length(object@comp) > 0){
                if (is.null(names(object@comp)))
                  return("comp names must be set")
                for (c in names(object@comp))
                  if (!is(object@comp[[c]], "BSRClusterComp"))
                    return("comp must contain objects of class BSRClusterComp")
              }
              if (!is.numeric(object@mu) || object@mu<=0)
                  return("mu must be numeric and >0")

              TRUE
            }
)

setMethod("show", "BSRDataModelComp",
          function(object) {
            callNextMethod()
            cat("mu: ", object@mu, "\n", sep="")
            cat("Defined comparisons:\n")
            utils::str(object@comp)
          }
)


# Accessors & setters ========================================================

if (!isGeneric("comp")) {
  if (is.function("comp"))
    fun <- comp
  else
    fun <- function(x) standardGeneric("comp")
  setGeneric("comp", fun)
}
#' Comparisons list accessor
#'
#' @name comp
#' @aliases comp,BSRDataModelComp-method
#' @param x object BSRDataModelComp 
#' @export
setMethod("comp", "BSRDataModelComp", function(x) x@comp)

if (!isGeneric("comp<-")) {
  if (is.function("comp<-"))
    fun <- `comp<-`
  else
    fun <- function(x, value) standardGeneric("comp<-")
  setGeneric("comp<-", fun)
}
#' Comparisons list setter (internal use only, use addComparison() otherwise)
#'
#' @param x object BSRDataModelComp 
#' @param value value to be set for BSRDataModelComp
#' @keywords internal 
setMethod("comp<-", "BSRDataModelComp", function(x,value){
  x@comp <- value
  methods::validObject(x)
  x
})


if (!isGeneric("mu")) {
  if (is.function("mu"))
    fun <- mu
  else
    fun <- function(x) standardGeneric("mu")
  setGeneric("mu", fun)
}
#' Mu accessor
#'
#' @name mu
#' @aliases mu,BSRDataModelComp-method
#' @param x object BSRDataModelComp 
#' @export
setMethod("mu", "BSRDataModelComp", function(x) x@mu)

if (!isGeneric("mu<-")) {
  if (is.function("mu<-"))
    fun <- `mu<-`
  else
    fun <- function(x, value) standardGeneric("mu<-")
  setGeneric("mu<-", fun)
}
#' Mu setter (internal use only)
#'
#' @param x object BSRDataModelComp 
#' @param value value to be set for BSRDataModelComp
#' @keywords internal 
setMethod("mu<-", "BSRDataModelComp", function(x,value){
  x@mu <- value
  methods::validObject(x)
  x
})


# generation of BSRDataModelComp from BSRDataModel =====================

#' @title conversion of BSRDataModel into BSRDataModelComp 
#'
#' @description In case ligand-receptor inferences should be obtained
#' based on gene/protein regulation P-values comparing two clusters of
#' samples, it is necessary to first promote the BSRDataModel object that
#' contains the count matrix into a BSRDataModelComp object able to contain
#' a list of such cluster pairs comparisons. This function performs this
#' promotion adding an empty list of comparisons.
#' @param bsrdm    A BSRDataModel object.
#'
#' @export
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' 
as.BSRDataModelComp <- function(bsrdm){
  
  if (!is(bsrdm, "BSRDataModel"))
    stop("bsrdm must be of class BSRDataModel")
  m <- mean(ncounts(bsrdm))
  if (!logTransformed(bsrdm))
    m <- log1p(m) # approximate of mu on the log1p-transformed matrix
  new("BSRDataModelComp", bsrdm, comp=list(), mu=m)
  
} # as.BSRDataModelComp


# defining/adding a cluster comparison ========================================

if (!isGeneric("defineClusterComp")) {
  if (is.function("defineClusterComp"))
    fun <- defineClusterComp
  else
    fun <- function(obj, ...) standardGeneric("defineClusterComp")
  setGeneric("defineClusterComp", fun)
}
#' Definition of the comparison between two clusters of samples
#'
#' Define the columns of the expression matrix that belong to each cluster,
#' and store the result of the cluster differences statistical analysis
#' obtained by an external tool such as edgeR, DESeq2, etc.
#'
#' @name defineClusterComp
#' @aliases defineClusterComp,BSRDataModelComp-method
#'
#' @param obj    A BSRDataModelComp object output by
#'   \code{\link{as.BSRDataModelComp}}.
#' @param colA   Cluster A column indices.
#' @param colB   Cluster B column indices.
#' @param stats  A data.frame containing statistics about the differential
#' analysis cluster A versus B. \code{stats} must contain at least the
#' columns 'pval' (for P-values), 'logFC' for log-fold-changes A/B, and
#' 'expr' for the expression of the genes in cluster A.
#'
#' @details Create a BSRClusterComp object describing a comparison
#' of two clusters of columns taken from the expression matrix
#' in the BSRDataModelComp object \code{obj}. Such a cluster comparison
#' description is the basis for inferring LRIs from differential
#' expression P-values instead of correlation analysis.
#' 
#' The rows of \code{stats} must be in the same order as those of the count
#' matrix in \code{obj}. Alternatively, \code{stats}
#' rows can be named and a 1-1 correspondence must exist between these names
#' and those of the count matrix.
#'
#' @return A BSRClusterComp object.
#'
#' @export
#'
#' @examples
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
#' bsrinf <- initialInference(bsrdm.comp, max.pval=1,"random.example")
#' 
#' @importFrom methods new
setMethod("defineClusterComp", "BSRDataModelComp", function(obj, colA,
                                                        colB, stats){
  
  if (!is.integer(colA))
    stop("colA must contain integer indices")
  if (!is.integer(colB))
    stop("colB must contain integer indices")
  if (length(intersect(colA, colB)) > 0)
    stop("colA and colB must be disjoint")
  if (any(colA < 1 | colA > ncol(ncounts(obj))))
    stop("colA indices must fall in [1; ncol(ncounts)]")
  if (any(colB < 1 | colB > ncol(ncounts(obj))))
    stop("colB indices must fall in [1; ncol(ncounts)]")
  if (!is.data.frame(stats))
    stop("stats must be a data.frame")
  if (!all(c("pval","logFC","expr") %in% names(stats)))
    stop("stats data.frame must contain columns named 'pval', 'logFC', and 'expr'")
  if (nrow(stats) != nrow(ncounts(obj)))
    stop("stats and ncounts(obj) number of rows differ")
  if (!is.null(rownames(stats)) &&
      (sum(rownames(stats) %in% rownames(ncounts(obj))) != nrow(stats)))
    stop("stats rownames defined but do not all match ncounts(obs)")
  if (is.null(rownames(stats)))
    rownames(stats) <- rownames(ncounts(obj))
  
  new("BSRClusterComp", colA=colA, colB=colB, stats=stats)
  
}) # defineClusterComp


if (!isGeneric("addClusterComp")) {
  if (is.function("addClusterComp"))
    fun <- addClusterComp
  else
    fun <- function(obj, ...) standardGeneric("addClusterComp")
  setGeneric("addClusterComp", fun)
}
#' Add a comparison between two clusters of samples
#'
#' Add a comparison to a BSRDataModelComp object.
#'
#' @name addClusterComp
#' @aliases addClusterComp,BSRDataModelComp-method
#'
#' @param obj    A BSRDataModelComp object output by
#'   \code{\link{as.BSRDataModelComp}}.
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
#' @importFrom methods new
setMethod("addClusterComp", "BSRDataModelComp", function(obj, cmp,
                                                         cmp.name){
  
  if (!is(cmp, "BSRClusterComp"))
    stop("cmp must be of class BSRClusterComp")
  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (cmp.name %in% names(comp(obj)))
    stop("cmp.name is already in the list of comparisons")
  
  tmp <- c(comp(obj), list(cmp))
  names(tmp)[length(tmp)] <- cmp.name
  comp(obj) <- tmp
  obj
  
}) # addClusterComp


if (!isGeneric("removeClusterComp")) {
  if (is.function("removeClusterComp"))
    fun <- removeClusterComp
  else
    fun <- function(obj, ...) standardGeneric("removeClusterComp")
  setGeneric("removeClusterComp", fun)
}
#' Remove a comparison from a BSRDataModelComp object.
#'
#' @name removeClusterComp
#' @aliases removeClusterComp,BSRDataModelComp-method
#'
#' @param obj    A BSRDataModelComp object output by
#'   \code{\link{as.BSRDataModelComp}}.
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
#' 
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' bsrdm.comp
#' bsrdm.comp <- removeClusterComp(bsrdm.comp, "random.example")
#' bsrdm.comp
#' 
setMethod("removeClusterComp", "BSRDataModelComp", function(obj, cmp.name){

  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (!(cmp.name %in% names(comp(obj))))
    stop("cmp.name must be in the list of comparisons")
  
  if (length(comp(obj)) == 1)
    comp(obj) <- list()
  else
    comp(obj) <- comp(obj)[names(comp(obj)) != cmp.name]
  
  obj
  
}) # removeClusterComp


# performing initial inference ===========================================


if (!isGeneric("initialInference")) {
  if (is.function("initialInference"))
    fun <- initialInference
  else
    fun <- function(obj, ...) standardGeneric("initialInference")
  setGeneric("initialInference", fun)
}
#' Inference of ligand-receptor interactions based on regulation
#'
#' This method supports two configurations that we refer to
#' as paracrine and autocrine.
#' 
#' In the autocrine case, a single cluster comparison name is provided.
#' In the corresponding cluster comparison, a group of samples A was
#' compared to a group of samples B to determine fold-changes and associated
#' P-values. The inferred ligand-receptor interactions take place in the
#' samples of group A. They are paracrine interactions in the case of
#' single-cell data or they take place in the same tissue represented by
#' cluster A. A typical single-cell example would be a population of
#' macrophages (group A) compared to all the other populations (group B) to
#' represent specific increased or decreased expression in macrophages. The
#' resulting ligand-receptor interactions will be autocrine interactions
#' that are exacerbated (or reduced depending on the chosen parameters) in
#' macrophages.
#' 
#' In the paracrine case, two cluster comparison names must be provided.
#' For instance, a first comparison coul involved macrophages versus all
#' the other cell populations as above. The second comparison could be
#' B-cells against all the other populations. Now, calling initialInference()
#' with comparison macrophages vs. the rest and, as source comparison, B-cells
#' vs. the rest, will result in inferring interactions between B-cells
#' (ligands) and macrophages (receptors and downstream pathways). To obtain
#' macrophages to B-cells paracrine interactions, it is necessary to call the
#' method a second time with permuted cluster comparison names. Another example
#' in spatial transcriptomics could be two thin bands at the boundary of two
#' tissue regions, one emitting the ligand and the other one expressing the
#' receptor.
#'
#' In this initial inference, all the receptor-containing pathways are reported,
#' see reduction functions to reduce this list.
#' 
#' @name initialInference
#' @aliases initialInference,BSRDataModelComp-method
#'
#' @param obj       A BSRDataModelComp object.
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
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
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
#' @param min.t.logFC     The minimum log2 fold-change allowed for
#'   targets in case pos.targets or neg.targets are used.
#' @param pos.targets   A logical imposing that all the network targets must
#'   display positive logFC, i.e. logFC >= min.t.logFC.
#' @param neg.targets   A logical imposing that all the network targets must
#'   display negative logFC, i.e. logFC <= - min.t.logFC.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#' @param use.full.network  A logical to avoid limiting the reference network
#' to the detected genes and use the whole reference network.
#'
#' @details Perform the initial ligand-receptor inference. Initial means that
#' no reduction is applied. All the (ligand, receptor, downstream pathway)
#' triples are reported, i.e., a given LR pair may appear multiple times
#' with different pathways downstream the receptor. Specific reduction
#' functions are available from the package to operate subsequent
#' simplifications based on the BSRInferenceComp object created by this method.
#' 
#' Here, ligand-receptor interactions are inferred based on gene or protein
#' regulation-associated P-values when comparing two clusters of samples. Since
#' a BSRDataModelComp object can contain several such comparisons, the name
#' of the comparison to use must be specified (parameter \code{cmp.name}).
#' 
#' Note that since the introduction of the use.full.network parameter
#' (April 29, 2024), the pathway sizes are always computed before potential
#' intersection with the observed data (use.full.network set to FALSE) for
#' consistency. Accordingly, the minimum and maximum pathway default values
#' have been raised from 5 & 200 to 10 & 600 respectively. By default,
#' use.full.network is set to FALSE.
#' 
#' In addition to statistical significance estimated according to BulkSignalR
#' statistical model, we compute SingleCellSignalR LR-score, L and R
#' cluster average expression. In the paracrine case, L average expression
#' is taken from the source cluster.
#'
#' @return A BSRInferenceComp object with initial inferences set.
#'
#' @export
#'
#' @examples
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
#' @importFrom methods new
setMethod("initialInference", "BSRDataModelComp", function(obj, cmp.name, src.cmp.name=NULL, rank.p=0.55,
                                                         max.pval=0.01, min.logFC=1, neg.receptors=FALSE,
                                                         pos.targets=FALSE, neg.targets=FALSE,
                                                         min.t.logFC=0.5, restrict.genes=NULL,
                                                         use.full.network=FALSE,
                                                         reference=c("REACTOME-GOBP","REACTOME","GOBP"),
                                                         max.pw.size=600, min.pw.size=10, min.positive=4,
                                                         restrict.pw=NULL, with.complex=TRUE,
                                                         fdr.proc=c("BH","Bonferroni","Holm","Hochberg",
                                                                    "SidakSS","SidakSD","BY","ABH","TSBH")){
  
  if (!(cmp.name %in% names(comp(obj))))
    stop("cmp.name must exist in the names of comparisons contained in obj")
  if (!is.null(src.cmp.name) && !(src.cmp.name %in% names(comp(obj))))
    stop("src.cmp.name must exist in the names of comparisons contained in obj")
  reference <- match.arg(reference)
  fdr.proc <- match.arg(fdr.proc)
  if (min.logFC <= 0)
    stop("min.logFC must be >0")
  if (min.t.logFC <= 0)
    stop("min.t.logFC must be >0")
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  if (neg.targets && pos.targets)
    stop("neg.targets and pos.targets cannot be TRUE simultaneously")
  
  # retrieve the BSRClusterComp object(s)
  cc <- comp(obj)[[cmp.name]]
  if (!is.null(src.cmp.name))
    scc <- comp(obj)[[src.cmp.name]]
  else
    scc <- NULL

  # store inference parameters and retrieve relevant L & R
  inf.param <- list()
  inf.param$log.transformed.data <- logTransformed(obj)
  inf.param$mu <- mu(obj)
  inf.param$colA <- colA(cc)
  inf.param$colB <- colB(cc)
  if (is.null(src.cmp.name))
    inf.param$inference.type <- "autocrine"
  else{
    inf.param$inference.type <- "paracrine"
    inf.param$src.colA <- colA(scc)
    inf.param$src.colB <- colB(scc)
  }
  inf.param$max.pval <- max.pval
  inf.param$min.logFC <- min.logFC
  inf.param$neg.receptors <- neg.receptors
  inf.param$pos.targets <- pos.targets
  inf.param$neg.targets <- neg.targets
  inf.param$min.t.logFC <- min.t.logFC
  inf.param$restrict.genes <- restrict.genes
  inf.param$use.full.network <- use.full.network
  lr <- .getRegulatedLR(obj, cc, scc, max.pval=max.pval, min.logFC=min.logFC,
                        neg.receptors=neg.receptors, restrict.genes=restrict.genes)
  
  # apply BSR model on the targets
  inf.param$reference <- reference
  inf.param$min.pw.size <- min.pw.size
  inf.param$max.pw.size <- max.pw.size
  inf.param$with.complex <- with.complex
  inf.param$min.positive <- min.positive
  inf.param$restrict.pw <- restrict.pw
  pairs <- .checkRegulatedReceptorSignaling(obj, cc, lr, reference=reference,
                                            pos.targets=pos.targets, neg.targets=neg.targets,
                                            min.t.logFC=min.logFC, use.full.network=use.full.network,
                                            min.pw.size=min.pw.size, max.pw.size=max.pw.size,
                                            min.positive=min.positive, with.complex=with.complex,
                                            restrict.pw=restrict.pw)
  inf.param$fdr.proc <- fdr.proc
  inf.param$rank.p <- rank.p
  
  # compute P-values
  inter <- .pValuesRegulatedLR(pairs, param(obj), rank.p=rank.p, fdr.proc=fdr.proc)
  
  # compute LR-score for compatibility with SingleCellSignalR version 1
  if (is.null(scc))
    inter$L.expr <- stats(cc)[inter$L, "expr"]
  else
    inter$L.expr <- stats(scc)[inter$L, "expr"]
  inter$R.expr <- stats(cc)[inter$R, "expr"]
  if (inf.param$log.transformed.data)
    sq <- sqrt(inter$L.expr * inter$R.expr)
  else
    sq <- sqrt(log1p(inter$L.expr) * log1p(inter$R.expr))
  inter$LR.score <- sq/(inf.param$mu+sq)

  # prepare the accompanying lists  
  ligands <- strsplit(inter$L, ";")
  receptors <- strsplit(inter$R, ";")
  tg <- strsplit(inter$target.genes, ";")
  tgpval <- lapply(strsplit(inter$target.pval, ";"),
                   function(x) as.numeric(x))
  tglogFC <- lapply(strsplit(inter$target.logFC, ";"),
                   function(x) as.numeric(x))
  tgcorr <- lapply(strsplit(inter$target.corr, ";"),
                   function(x) as.numeric(x))
  tgexpr <- lapply(strsplit(inter$target.expr, ";"),
                   function(x) as.numeric(x))
  inf.param$ligand.reduced <- FALSE
  inf.param$receptor.reduced <- FALSE
  inf.param$pathway.reduced <- FALSE
  
  # instantiate the object
  if (is.null(src.cmp.name))
    src.cmp.name.char <- ""
  else
    src.cmp.name.char <- src.cmp.name
  new("BSRInferenceComp", LRinter=inter[,c("L","R","pw.id","pw.name","pval","qval",
                                           "L.logFC","R.logFC","LR.pval","LR.corr",
                                           "rank","len","rank.pval","rank.corr",
                                           "LR.score","L.expr","R.expr")],
      ligands=ligands, receptors=receptors, t.genes=tg, tg.corr=tgcorr,
      tg.pval=tgpval, tg.logFC=tglogFC, tg.expr=tgexpr, inf.param=inf.param,
      cmp.name=cmp.name, src.cmp.name=src.cmp.name.char)
  
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
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1,  "random.example")
#' 
#' # reduction
#' bsrinf.red <- reduceToBestPathway(bsrinf)
#' 
#' # signature extraction and scoring
#' bsrsig.red <- getLRGeneSignatures(bsrinf.red, qval.thres=1e-6)
#' scores.red <- scoreLRGeneSignatures(bsrdm.comp, bsrsig.red,
#'                          name.by.pathway=TRUE, rownames.LRP=TRUE)
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom methods is
setMethod("scoreLRGeneSignatures", "BSRDataModelComp", function(obj,
                                                            sig, LR.weight=0.5, robust=FALSE,
                                                            name.by.pathway=FALSE, abs.z.score=FALSE,rownames.LRP=FALSE){
  
  if (!is(sig, "BSRSignatureComp"))
    stop("sig must be a BSRSignature object")
  if (LR.weight<=0 || LR.weight>=1)
    stop("LRweight must reside in (0;1)")

  # retrieve the BSRClusterComp object
  cmp.name <- cmpName(sig)
  if (!(cmp.name %in% names(comp(obj))))
    stop("The comparison name in sig is not present in obj")
  cc <- comp(obj)[[cmp.name]]

  # species management  
  if( initialOrganism(obj)!="hsapiens" )
    all.genes <- unlist(initialOrthologs(obj))
  else
    all.genes <- rownames(ncounts(obj))
  
  # get the ncount matrix with proper columns
  ncounts <- ncounts(obj)[, c(colA(cc), colB(cc))]
  
  # intersect signature gene names with RNA-seq data
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
