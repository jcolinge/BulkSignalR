#' Get correlated ligand-receptor pairs
#'
#' Compute the Spearman correlations of all the ligand-receptor pairs in LRdb and return those above a minimum value.
#'
#' @param ds              A BulkSignalR data set.
#' @param min.cor         The minimum correlation required.
#' @param restrict.genes  A list of gene symbols that restricts ligands and receptors.
#' @return A list containing a table of putative ligand-receptor pairs along with
#' their correlations above \code{min.cor}. This list is the first step of a ligand-receptor analysis.
#'
#' The \code{restrict.genes}
#' parameter is used for special cases where LRdb must be further restricted to a subset.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds.LR <- getCorrelatedLR(ds)
#' }
getCorrelatedLR <- function(ds,min.cor=0.3,restrict.genes=NULL){

  if ((min.cor < -1) || (min.cor > 1))
    stop("min.cor must lie in [-1;+1]")

  lrgenes <- intersect(c(SingleCellSignalR::LRdb$ligand,SingleCellSignalR::LRdb$receptor),rownames(ds$ncounts))
  if (!is.null(restrict.genes))
    lrgenes <- intersect(lrgenes,restrict.genes)
  corlr <- stats::cor(t(ds$ncounts[lrgenes,]),method="spearman")
  pairs <- foreach::foreach(i=1:nrow(SingleCellSignalR::LRdb),.combine=rbind) %do% {
    if (SingleCellSignalR::LRdb$ligand[i]%in%rownames(corlr) && SingleCellSignalR::LRdb$receptor[i]%in%rownames(corlr))
      data.frame(L=SingleCellSignalR::LRdb$ligand[i],R=SingleCellSignalR::LRdb$receptor[i],
                 corr=corlr[SingleCellSignalR::LRdb$ligand[i],SingleCellSignalR::LRdb$receptor[i]],stringsAsFactors=FALSE)
    else
      NULL
  }
  good <- pairs$corr>=min.cor
  sum(good)
  list(putative.pairs=pairs[good,])

} # getCorrelatedLR


#' Internal function to check receptor signaling downstream
#'
#' @param lr              A ligand-receptor analysis such as returned by \code{getCorrelatedLR()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway IDs).
#' @param rncounts        A matrix or table of normalized read counts with at least all the ligands, receptors, and genes in the reference pathways.
#' @param id.col          Column index or name in \code{pw} for the pathway IDs.
#' @param gene.col        Column index or name in \code{pw} for the gene symbols.
#' @param pw.col          Column index or name in \code{pw} for the pathway names.
#' @param min.positive    Minimum number of target genes to be found in a given pathway.
#' @param with.complex    A logical indicating whether receptor co-complex members should be included in the target genes.
#' @return A table reporting all the ligand-receptor pairs provided in \code{lr}
#' along with the pathways found and data about target gene correlations with the receptor.
#'
#' This BulkSignalR internal function is made visible for advanced users who might want to use it to implement
#' alternative analyses.
#' @export
.downstreamSignaling <- function(lr,pw,pw.size,rncounts,id.col,gene.col,pw.col,min.positive,with.complex=TRUE){

  # define interaction types
  control.int <- "controls-expression-of"
  incomplex.int <- "in-complex-with"
  directed.int <- c("controls-state-change-of","catalysis-precedes","controls-expression-of","controls-transport-of","controls-phosphorylation-of")
  if (with.complex)
    correlated.int <- union(control.int,incomplex.int)
  else
    correlated.int <- control.int

  # compute downstream correlations
  corrg <- stats::cor(t(rncounts),method="spearman")
  # global computation above is faster than for the receptors only
  corrg <- corrg[unique(lr$putative.pairs$R),]
  reg.proc <- foreach::foreach(r=unique(lr$putative.pairs$R),.combine=rbind) %dopar% {
    # reg.proc <- NULL
    # for (r in unique(lr$putative.pairs$R)){
    pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
    if (length(pa)>0){
      best.2nd <- foreach::foreach(p=pa,.combine=rbind) %do% {
        # best.2nd <- NULL
        # for (p in pa){
        int <- SingleCellSignalR::PwC_ReactomeKEGG[SingleCellSignalR::PwC_ReactomeKEGG$a.gn%in%pw[pw[[id.col]]==p,gene.col] & SingleCellSignalR::PwC_ReactomeKEGG$b.gn%in%pw[pw[[id.col]]==p,gene.col],]
        directed <- int$type%in%directed.int
        d.int <- unique(rbind(int[directed,c("a.gn","b.gn")],int[!directed,c("a.gn","b.gn")],int[!directed,c("b.gn","a.gn")]))
        g <- igraph::graph_from_data_frame(d.int,directed=TRUE)
        if (r%in%d.int$a.gn || r%in%d.int$b.gn){
          target.genes <- setdiff(c(int[int$type%in%incomplex.int & int$a.gn==r,"b.gn"],int[int$type%in%incomplex.int & int$b.gn==r,"a.gn"],int[int$type%in%directed.int,"b.gn"]),r)
          sp <- igraph::shortest.paths(g,r,target.genes)
          target.genes <- colnames(sp)[!is.infinite(sp[r,])]
          if (length(target.genes)>=min.positive){
            c <- corrg[r,target.genes]
            o <- order(c)
            c <- c[o]
            target.genes <- target.genes[o]
            data.frame(pathway=p,corr.2nd=c[length(c)-1],all.corr=paste(c,collapse=";"),target.genes=paste(target.genes,collapse=";"),len=length(c),stringsAsFactors=FALSE)
          }
          else
            NULL
        }
        else
          NULL
      }
      if (!is.null(best.2nd))
        data.frame(R=r,corr=max(best.2nd$corr.2nd),pathways=paste(best.2nd$pathway,collapse="|"),
                   all.corr=paste(best.2nd$all.corr,collapse='|'),target.genes=paste(best.2nd$target.genes,collapse='|'),
                   len=paste(best.2nd$len,collapse='|'),stringsAsFactors=FALSE)
      else
        NULL
    }
    else
      NULL
  }
  rownames(reg.proc) <- reg.proc$R
  conf.pairs <- lr$putative.pairs[lr$putative.pairs$R%in%reg.proc$R,]
  conf.pairs$pwid <- reg.proc[conf.pairs$R,"pathways"]
  conf.pairs$all.corr <- reg.proc[conf.pairs$R,"all.corr"]
  conf.pairs$len <- reg.proc[conf.pairs$R,"len"]
  conf.pairs$corr.pw <- reg.proc[conf.pairs$R,"corr"]
  conf.pairs$target.genes <- reg.proc[conf.pairs$R,"target.genes"]
  pw.name <- unique(pw[,c(id.col,pw.col)])
  pw2name <- stats::setNames(pw.name[[2]],pw.name[[1]])
  conf.pairs$pwname <- foreach::foreach(pl=conf.pairs$pwid,.combine=c) %do% {
    paste(foreach::foreach (id=unlist(strsplit(pl,"\\|")),.combine=c)%do%{pw2name[id]},collapse="|")
  }

  conf.pairs

} # .downstreamSignaling


#' Check receptor signaling downstream
#'
#' Assess the existence of correlations between a receptor, part of a ligand-receptor pair, and
#' genes coding for proteins forming a complex with the receptor or genes regulated by the
#' receptor downstream signaling.
#'
#' @param ds              A BulkSignalR data set.
#' @param lr              A ligand-receptor analysis such as returned by \code{getCorrelatedLR()}.
#' @param method          Which pathway reference should be used ("reactome" for Reactome,"GOBP" for GO Biological Process,
#'                        or "Reactome-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway reference.
#' @param min.positive    Minimum number of target genes to be found in a given pathway.
#' @param restrict.pw     A list of pathway IDs to restrict the application of the function.
#' @param with.complex    A logical indicating whether receptor co-complex members should be included in the target genes.
#' @return A ligand-receptor analysis extending \code{lr} content with a table reporting all the ligand-receptor pairs provided in \code{lr}
#' along with the pathways found and data about target gene correlations with the receptor. Those data are the lists
#' (in semi-colon-separated format) of target genes and their Spearman correlations with the receptor. The lists are sorted
#' according to the correlation coefficient.
#'
#' In a pathway of the reference, i.e., a Reactome pathway or the genes of a GOBP term, the target genes are the
#' genes coding for proteins forming a complex with the receptor and the genes in the pathway downstream the receptor,
#' which are given as regulated by the pathway. If \code{with.complex} is set to \code{FALSE}, then only the
#' regulated genes are considered. Participation to a complex and being regulated as well as the pathway directed topologies
#' are defined by Reactome and KEGG pathways as provided by PathwayCommons.
#'
#' The maximum pathway size is used to limit the redundancy inherent to GOBP and Reactome. The minimum pathway size is
#' used to avoid overspecific, noninformative results.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' }
checkReceptorSignaling <- function(ds,lr,method=c("reactome-GOBP","reactome","GOBP")[1],max.pw.size=200,min.pw.size=5,min.positive=4,restrict.pw=NULL,with.complex=TRUE){

  results <- lr
  method <- toupper(method)
  if (!(method %in% c("REACTOME-GOBP","REACTOME","GOBP")))
    stop("invalid method parameter")
  results$method <- method

  # reactome pathways
  if (method %in% c("REACTOME-GOBP","REACTOME")){
    react <- reactome[reactome$`Gene name`%in%rownames(ds$ncounts),]
    if (!is.null(restrict.pw))
      react <- react[react$`Reactome ID`%in%restrict.pw,]
    pw.size <- table(react$`Reactome ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- react[react$`Gene name`%in%lr$putative.pairs$R,"Reactome ID"]
    pw.size <- pw.size[names(pw.size)%in%contains.receptors]
    corgenes <- unique(c(lr$putative.pairs$R,react[react$`Reactome ID`%in%names(pw.size),"Gene name"]))
    results$reactome.pairs <- .downstreamSignaling(lr,react,pw.size,ds$ncounts[corgenes,],"Reactome ID","Gene name","Reactome name",min.positive,with.complex=with.complex)
  }

  # GOBP
  if (method %in% c("REACTOME-GOBP","GOBP")){
    go <- gobp[gobp$`Gene name`%in%rownames(ds$ncounts),]
    if (!is.null(restrict.pw))
      go <- go[go$`GO ID` %in%restrict.pw,]
    pw.size <- table(go$`GO ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- go[go$`Gene name`%in%lr$putative.pairs$R,"GO ID"]
    pw.size <- pw.size[names(pw.size)%in%contains.receptors]
    corgenes <- unique(c(lr$putative.pairs$R,go[go$`GO ID`%in%names(pw.size),"Gene name"]))
    results$gobp.pairs <- .downstreamSignaling(lr,go,pw.size,ds$ncounts[corgenes,],"GO ID","Gene name","GO name",min.positive,with.complex=with.complex)
  }

  # merge
  if (method == "REACTOME-GOBP"){
    pairs <- unique(rbind(results$reactome.pairs[,1:2],results$gobp.pairs[,1:2]))
    react.keys <- paste(results$reactome.pairs[[1]],results$reactome.pairs[[2]],sep="|")
    gobp.keys <- paste(results$gobp.pairs[[1]],results$gobp.pairs[[2]],sep="|")
    results$merged.pairs <- rbind(results$reactome.pairs,results$gobp.pairs[!(gobp.keys%in%react.keys),])
  }
  else if (method == "REACTOME")
    results$merged.pairs <- results$reactome.pairs
  else
    results$merged.pairs <- results$gobp.pairs

  results

} # checkReceptorSignaling


