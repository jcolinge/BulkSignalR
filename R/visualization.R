#' Basic statistics about hit pathways
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @param LLR.thres     Log-likelihood threshold.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @return A table with the pathways selected after the chosen threshold was applied to rows in \code{pairs}.
#' Each pathway is reported along with various statistics:
#' the number of selected receptors in this pathway, the total number of receptors described this pathway,
#' the number of selected ligand-receptor pairs hitting this pathway, and the total number of
#' ligand-receptor pairs described taht could hit this pathway.
#'
#' Obviously, one could imagine computing enrichments in receptors or ligand-receptor pairs
#' based on such statistics, but the actual meaning of such an analysis would be ambiguous since
#' the pathways were already selected as significantly regulated by the receptor. We thus did not implement
#' this (elementary) computation.
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
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' pw.stat <- getPathwayStats(pairs.p,qval.thres=0.01)
#' }
getPathwayStats <- function(pairs,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")

  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>=LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<=pval.thres,]
    else
      pairs <- pairs[pairs$qval<=qval.thres,]

    # pairs reduced to the receptor & pathway names
    pairs.R <- unique(pairs[,c("R","pw.id","pw.name")])
    upw <- unique(pairs[,c("pw.id","pw.name")])
    pw.names <- stats::setNames(upw$pw.name,upw$pw.id)

    # number of "hits" on a pathway with or without the ligand combinatorial contribution
    pw.ids <- table(pairs$pw.id)
    pw.ids.R <- table(pairs.R$pw.id)

    # number of ligands for each receptor
    R.n.comb <- table(SingleCellSignalR::LRdb$receptor)

    foreach::foreach(id=names(pw.ids),.combine=rbind) %do% {

      # number of receptors that are in the current pathway, depending on whether it is a GOBP or Reactome pathway
      if (regexpr("^R-",id) != -1)
        Rs <- intersect(SingleCellSignalR::LRdb$receptor,reactome[reactome$`Reactome ID`==id,"Gene name"])
      else
        Rs <- intersect(SingleCellSignalR::LRdb$receptor,gobp[gobp$`GO ID`==id,"Gene name"])

      # non-combinatorial version (ignore ligands)
      tot.R <- length(Rs)

      # ligand combinatorics included
      tot.LR <- sum(R.n.comb[Rs])

      data.frame(pw.id=id,pw.name=pw.names[id],n.R=pw.ids.R[id],tot.R=tot.R,n.LR=pw.ids[id],tot.LR=tot.LR,stringsAsFactors=FALSE)
    }

} # getPathwayStats


#' Extract gene signatures of LR pair activity
#'
#' Obtains gene signatures reflecting ligand-receptor as well as receptor downstream activity to
#' score ligand-receptor pairs across samples subsequently with \code{scoreLRGeneSignatures()}.
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @param LLR.thres     Log-likelihood threshold.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param signed        A logical indicating whether correlations were be considered with their signs upon P-value or LLR estimations.
#' @return A table with a gene signature for each triple ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}. In case \code{signed==TRUE},
#' the rank is defined for correlation absolute values.
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
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' signatures <- getPathwayStats(pairs.p,qval.thres=0.01)
#' }
getLRGeneSignatures <- function(pairs,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL,signed=TRUE){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")

  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>=LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<=pval.thres,]
    else
      pairs <- pairs[pairs$qval<=qval.thres,]
    LR <- reduceToBestPathway(pairs)

    foreach::foreach(i=1:nrow(LR),.combine=rbind) %do% {
      tg <- unlist(strsplit(LR$target.genes[i],split=";"))
      if (signed)
        signature <- tg[LR$rank[i]:length(tg)]
      else{
        corr <- as.numeric(unlist(strsplit(LRl$all.corr[i],split=";")))
        o <- order(corr**2)
        signature <- tg[o][LR$rank[i]:length(tg)]
      }
      data.frame(L=LR$L[i],R=LR$R[i],best.pw.id=LR$pw.id[i],best.pw.name=LR$pw.name[i],signature=paste(signature,collapse=";"),stringsAsFactors=FALSE)
    }

} # getLRGeneSignatures


#' Score ligand-receptor gene signatures
#'
#' Compute ligand-receptor gene signature scores over a BulkSIgnalR data set.
#'
#' @param ds         A BulkSignalR data set.
#' @param sig           Gene signatures computed by \code{getLRGeneSignatures()}.
#' @param sample.types  A character vector defining which samples in \code{ds} should be used.
#' @param LR.weight    A number between 0 and 1 defining the relative weight of the ligand and the receptor in the signature.
#' @param robust       A logical indicating that z-scores should be computed with median and MAD instead of mean and standard deviation.
#' @param rename.by.pathway     A logical indicating whether row names of the result score matrix should be pathway names instead of ligand-receptor gene symbols.
#' @return A matrix containing the scores of each ligand-receptor gene signature in each sample.
#'
#' A selection of samples can be obtained by using the \code{sample.types} parameter.
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
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' p.red.P <- reduceToPathway(pairs.p)
#' signatures <- getPathwayStats(p.red.P,qval.thres=0.01)
#' scores <- scoreLRGeneSignatures(ds,signatures,rename.by.pathway=TRUE)
#' }
scoreLRGeneSignatures <- function(ds,sig,sample.types=NULL,LR.weight=0.5,robust=FALSE,rename.by.pathway=FALSE){

  if (LR.weight<=0 || LR.weight>=1)
    stop("LRweight must reside in (0;1)")

  # check whether ligands or receptors were reduced
  red.mode.L <- (sum(regexpr("^{",sig$L,perl=TRUE)==1)==nrow(sig) && sum(regexpr("}$",sig$L,perl=TRUE)==nchar(sig$L))==nrow(sig))
  red.mode.R <- (sum(regexpr("^{",sig$R,perl=TRUE)==1)==nrow(sig) && sum(regexpr("}$",sig$R,perl=TRUE)==nchar(sig$R))==nrow(sig))

  # intersect signature gene names with RNA-seq data
  tg <- foreach::foreach(i=1:nrow(sig),.combine=c) %do% {
    list(intersect(rownames(ds$ncounts),unlist(strsplit(sig$signature[i],split=";"))))
  }
  if (rename.by.pathway)
    names(tg) <- sig$best.pw.name
  else
    names(tg) <- paste(sig$L,sig$R,sep=" / ")
  good <- sapply(tg,length) > 0
  tg <- tg[good]
  sig <- sig[good,]
  if (red.mode.L)
    ligands <- unlist(strsplit(gsub("}$","",gsub("^{","",sig$L,perl=TRUE)),split=";"))
  else
    ligands <- sig$L
  if (red.mode.R)
    receptors <- unlist(strsplit(gsub("}$","",gsub("^{","",sig$R,perl=TRUE)),split=";"))
  else
    receptors <- sig$R
  pool <- unique(c(unlist(tg),receptors,ligands))

  # convert normalized counts into z-scores
  if (is.null(sample.types))
    sample.types <- unique(ds$types)
  if (ds$log.transformed)
    ncounts <- ds$ncounts
  else
    ncounts <- data.matrix(log10(1+ds$ncounts[pool,ds$types%in%sample.types]))
  if (robust)
    z <- (ncounts-apply(ncounts,1,stats::median))/apply(ncounts,1,stats::mad)
  else
    z <- (ncounts-rowMeans(ncounts))/apply(ncounts,1,stats::sd)

  # compute the LR gene signatures
  res <- matrix(0,nrow=nrow(sig),ncol=ncol(ncounts),dimnames=list(names(tg),colnames(ncounts)))
  if (red.mode.R || red.mode.L)
    for (i in 1:nrow(sig)){
      ligands <- unlist(strsplit(gsub("}$","",gsub("^{","",sig$L[i],perl=TRUE)),split=";"))
      if (length(ligands)==1)
        mL <- z[ligands,]
      else
        mL <- as.numeric(colSums(z[ligands,])/length(ligands))
      receptors <- unlist(strsplit(gsub("}$","",gsub("^{","",sig$R[i],perl=TRUE)),split=";"))
      if (length(receptors)==1)
        mR <- z[receptors,]
      else
        mR <- as.numeric(colSums(z[receptors,])/length(receptors))
      res[i,] <- as.numeric(LR.weight*0.5*(mL+mR)+(1-LR.weight)*colSums(z[tg[[i]],])/length(tg[[i]]))
    }
  else
    for (i in 1:nrow(sig))
      res[i,] <- as.numeric(LR.weight*0.5*(z[sig$L[i],]+z[sig$R[i],])+(1-LR.weight)*colSums(z[tg[[i]],])/length(tg[[i]]))

  res

} # scoreLRGeneSignatures


#' Gene signature scoring
#'
#' Scores generic gene signatures over the samples of a BulkSignalR data set.
#'
#' @param ds         A BulkSignalR data set.
#' @param ref.signatures           Gene signatures.
#' @param sample.types  A character vector defining which samples in \code{ds} should be used.
#' @param robust       A logical indicating that z-scores should be computed with median and MAD instead of mean and standard deviation.
#' @return A matrix containing the scores of each gene signature in each sample. Note that ligand-receptor gene
#' signature scores should be computed with \code{scoreLRGeneSignatures()} instead.
#'
#' A selection of samples can be obtained by using the \code{sample.types} parameter.
#' @export
#' @examples
#' \dontrun{
#' tme.signatures <- fread("TME-population-signatures.txt",data.table=FALSE)
#' tme.scores <- scoreSignatures(ds,tme.signatures)
#' }
scoreSignatures <- function(ds,ref.signatures,sample.types=NULL,robust=FALSE){

  ref.signatures <- ref.signatures[ref.signatures$gene%in%rownames(ds$ncounts),]
  pool <- unique(ref.signatures$gene)

  # convert normalized counts into z-scores
  if (is.null(sample.types))
    sample.types <- unique(ds$types)
  if (ds$log.transformed)
    ncounts <- ds$ncounts
  else
    ncounts <- data.matrix(log10(1+ds$ncounts[pool,ds$types%in%sample.types]))
  if (robust)
    z <- (ncounts-apply(ncounts,1,stats::median))/apply(ncounts,1,stats::mad)
  else
    z <- (ncounts-rowMeans(ncounts))/apply(ncounts,1,stats::sd)

  # compute the gene signature scores
  pop <- unique(ref.signatures$signature)
  sig <- matrix(0,nrow=length(pop),ncol=ncol(ncounts),dimnames=list(pop,colnames(ncounts)))
  for (p in pop)
    sig[p,] <- as.numeric(colSums(z[ref.signatures$gene[ref.signatures$signature==p],])/sum(ref.signatures$signature==p))

  sig

} # scoreSignatures


#' Internal function to cut extreme values from a matrix
#'
#' @param m         A matrix.
#' @param p          Proportion of top and bottom values for thresholding.
#' @return A matrix with values beyond top and bottom thresholds repaced by the latter thresholds.
.cutExtremeValues <- function(m,p){
  thres.lo <- stats::quantile(m,prob=p)
  thres.hi <- stats::quantile(m,prob=1-p)
  m[m>thres.hi] <- thres.hi
  m[m<thres.lo] <- thres.lo
  m
}


#' Heatmap function for LR scores
#'
#' Generate a PDF file with a heatmap representing ligand-receptor gene signature scores.
#'
#' @param mat.c         A matrix with the signature scores such as output by \code{scoreLRGeneSignatures()}.
#' @param file.name     A PDF file name.
#' @param dend.row       A precomputed row dendrogram.
#' @param dend.spl       A precompute sample (column) dendrogram.
#' @param cols           A vector of colors to use for the heatmap.
#' @param width          PDF width.
#' @param height         PDF height.
#' @param pointsize      PDF pointsize.
#' @param bottom.annotation  \code{ComplexHeatmap} package bottom annotations.
#' @param n.col.clust    Number of column clusters.
#' @param n.row.clust    Number of row clusters.
#' @param gap.size       Gap size between clusters.
#' @param cut.p          Proportion of top and bottom values for thresholding.
#' @param row.names      A logical to turn on/off the display of row (gene) names.
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen, it is advisable to provide a
#' PDF file name and plot in a file (much faster).
#'
#' Extreme values (top and bottom) can be replaced by global quantiles at \code{cut.p} and \code{1-cut.p}
#' to avoid color scales shrunk by a few outliers.
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap} package to propose a simple way
#' of representing signature scores. If more advanced features are needed or more graphic parameters
#' should be controled, users should implement their own function, e.g., strarting from the code
#' of \code{simpleHeatmap()}.
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
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' p.red.P <- reduceToPathway(pairs.p)
#' signatures <- getPathwayStats(p.red.P,qval.thres=0.01)
#' scores <- scoreLRGeneSignatures(ds,signatures,rename.by.pathway=TRUE)
#' simpleHeatmap(scores,"example.pdf",width=9,height=5)
#' }
simpleHeatmap <- function(mat.c,file.name=NULL,dend.row=NULL,dend.spl=NULL,cols=NULL,width,height,pointsize=4,bottom.annotation=NULL,
                          n.col.clust=0,n.row.clust=0,gap.size=0.5,cut.p=0.01,row.names=TRUE){

  if (!requireNamespace("ComplexHeatmap",quietly=TRUE))
    stop("Package \"ComplexHeatmap\" needed for this function to work. Please install it.")
  if (cut.p<0 || cut.p>0.1)
    stop("cut.p must lie in [0;0.1]")

  if (is.null(cols))
    cols <- c(grDevices::colorRampPalette(c("royalblue3","white"),space="Lab")(trunc(-100*min(mat.c))),grDevices::colorRampPalette(c("white","orange"),space="Lab")(trunc(100*max(mat.c)))[-1])

  if (cut.p!=0)
    mat.c <- .cutExtremeValues(mat.c,cut.p)
  if (is.null(dend.spl)){
    di.spl <- stats::dist(t(mat.c))
    hc.spl <- stats::hclust(di.spl,method="ward.D")
    dend.spl <- stats::as.dendrogram(hc.spl)
  }
  if (is.null(dend.row)){
    di.gene <- stats::dist(mat.c)
    hc.gene <- stats::hclust(di.gene,method="ward.D")
    dend.row <- stats::as.dendrogram(hc.gene)
  }

  if (!is.null(file.name))
    grDevices::pdf(file.name,width=width,height=height,pointsize=pointsize,useDingbats=FALSE)
  if (n.row.clust>0)
    if (n.col.clust>0)
      print(ComplexHeatmap::Heatmap(mat.c,cluster_rows=dend.row,cluster_columns=dend.spl,col=cols,show_row_names=row.names,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                    row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,bottom_annotation=bottom.annotation,
                    split=n.row.clust,gap=grid::unit(gap.size,"mm"),column_split=n.col.clust,column_gap=grid::unit(gap.size,"mm")))
  else
    print(ComplexHeatmap::Heatmap(mat.c,cluster_rows=dend.row,cluster_columns=dend.spl,col=cols,show_row_names=row.names,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                  row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,bottom_annotation=bottom.annotation,
                  split=n.row.clust,gap=grid::unit(gap.size,"mm")))
  else
    if (n.col.clust)
      print(ComplexHeatmap::Heatmap(mat.c,cluster_rows=dend.row,cluster_columns=dend.spl,col=cols,show_row_names=row.names,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                    row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,bottom_annotation=bottom.annotation,
                    column_split=n.col.clust,column_gap=grid::unit(gap.size,"mm")))
  else
    print(ComplexHeatmap::Heatmap(mat.c,cluster_rows=dend.row,cluster_columns=dend.spl,col=cols,show_row_names=row.names,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                  row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,bottom_annotation=bottom.annotation))
  if (!is.null(file.name))
    grDevices::dev.off()

} # simpleHeatmap


#' Heatmap function for LR scores and additional data
#'
#' Generate a stack of heatmaps. The top heatmap represents ligand-receptor gene signature scores,
#' while the bottom heatmap represents the second, user-chosen data set to feature along with gene signatures.
#'
#' @param mat.c         A matrix with the signature scores such as output by \code{scoreLRGeneSignatures()}.
#' @param mat.e         A second matrix to feature below mat.c.
#' @param file.name     A PDF file name.
#' @param dend.row       A precomputed row dendrogram for \code{mat.c}.
#' @param dend.e         A precomputed row dendrogram for \code{mat.e}.
#' @param dend.spl       A precompute sample (column) dendrogram.
#' @param cols           A vector of colors to use for the heatmap of \code{mat.c}.
#' @param cols.e           A vector of colors to use for the heatmap of \code{mat.e}.
#' @param width          PDF width.
#' @param height         PDF height.
#' @param pointsize      PDF pointsize.
#' @param cut.p          Proportion of top and bottom values for thresholding.
#' @param row.names.1      A logical to turn on/off the display of row (gene) names in the top heatmap.
#' @param row.names.2      A logical to turn on/off the display of row (gene) names in the bottom heatmap.
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen, it is advisable to provide a
#' PDF file name and plot in a file (much faster).
#'
#' Extreme values (top and bottom) can be replaced by global quantiles at \code{cut.p} and \code{1-cut.p}
#' to avoid color scales shrunk by a few outliers.
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap} package to propose a simple way
#' of representing signature scores. If more advanced features are needed or more graphic parameters
#' should be controled, users should implement their own function, e.g., strarting from the code
#' of \code{dualHeatmap()}.
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
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' p.red.P <- reduceToPathway(pairs.p)
#' signatures <- getPathwayStats(p.red.P,qval.thres=0.01)
#' scores <- scoreLRGeneSignatures(ds,signatures,rename.by.pathway=TRUE)
#' data(tme.signatures,package="BulkSignalR")
#' tme.scores <- scoreSignatures(ds,tme.signatures)
#' dualHeatmap(mat.c,mat.e,"example-with-TME.pdf",width=9,height=7,pointsize=4)
#' }
dualHeatmap <- function(mat.c,mat.e,file.name=NULL,dend.row=NULL,dend.spl=NULL,dend.e=NULL,cols=NULL,cols.e=NULL,width,height=6,pointsize=4,cut.p=0.01,vert.p=0.9,ros.names.1=TRUE,row.names.2=TRUE){

  if (!requireNamespace("ComplexHeatmap",quietly=TRUE))
    stop("Package \"ComplexHeatmap\" needed for this function to work. Please install it.")
  if (cut.p<0 || cut.p>0.1)
    stop("cut.p must lie in [0;0.1]")
  if (vert.p<0.05 || vert.p>0.95)
    stop("vert.p must lie in [0.05;0.95]")
  
  if (is.null(cols))
    cols <- c(grDevices::colorRampPalette(c("royalblue3","white"),space="Lab")(trunc(-100*min(mat.c))),grDevices::colorRampPalette(c("white","orange"),space="Lab")(trunc(100*max(mat.c)))[-1])
  if (is.null(cols.e))
    cols.e <- c(grDevices::colorRampPalette(c("deepskyblue","white"),space="Lab")(trunc(-100*min(mat.e))),grDevices::colorRampPalette(c("white","tomato"),space="Lab")(trunc(100*max(mat.e)))[-1])

  if (cut.p!=0){
    mat.c <- .cutExtremeValues(mat.c,cut.p)
    mat.e <- .cutExtremeValues(mat.e,cut.p)
  }
  if (is.null(dend.spl)){
    di.spl <- stats::dist(t(mat.c))
    hc.spl <- stats::hclust(di.spl,method="ward.D")
    dend.spl <- stats::as.dendrogram(hc.spl)
  }
  if (is.null(dend.row)){
    di.gene <- stats::dist(mat.c)
    hc.gene <- stats::hclust(di.gene,method="ward.D")
    dend.row <- stats::as.dendrogram(hc.gene)
  }
  if (is.null(dend.e)){
    di.e <- stats::dist(mat.e)
    hc.e <- stats::hclust(di.e,method="ward.D")
    dend.e <- stats::as.dendrogram(hc.e)
  }

  hm.LR <- ComplexHeatmap::Heatmap(mat.c,cluster_rows=dend.row,cluster_columns=dend.spl,col=cols,show_row_names=row.names.1,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                   row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,height=vert.p*height)
  hm.e <- ComplexHeatmap::Heatmap(mat.e,cluster_rows=dend.e,cluster_columns=dend.spl,col=cols.e,show_row_names=row.names.2,show_column_names=FALSE,use_raster=TRUE,raster_device="png",raster_quality=8,
                  row_names_gp=grid::gpar(fontsize=pointsize),show_row_dend=TRUE,height=(1-vert.p)*height)

  if (!is.null(file.name))
    grDevices::pdf(file.name,width=width,height=height,pointsize=pointsize,useDingbats=FALSE)
  import::from(ComplexHeatmap,"%v%")
  ComplexHeatmap::draw(hm.LR %v% hm.e,gap=grid::unit(1,"mm"))
  if (!is.null(file.name))
    grDevices::dev.off()

} # dualHeatmap