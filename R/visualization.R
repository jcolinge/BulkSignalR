# Pure display functions (for convenience) =====================================

#' Internal function to cut extreme values from a matrix
#'
#' @param m         A matrix.
#' @param p          Proportion of top and bottom values for thresholding.
#' @return A matrix with values beyond top and bottom thresholds repaced by
#' the latter thresholds.
.cutExtremeValues <- function(m, p){
    thres.lo <- stats::quantile(m, prob=p)
    thres.hi <- stats::quantile(m, prob=1-p)
    m[m>thres.hi] <- thres.hi
    m[m<thres.lo] <- thres.lo
    m
}


#' Heatmap function for LR scores
#'
#' Generate a PDF file with a heatmap representing ligand-receptor gene
#' signature scores.
#'
#' @param mat.c         A matrix with the signature scores such as output by
#' \code{scoreLRGeneSignatures()}.
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
#' @param row.names      A logical to turn on/off the display of row (gene)
#' names.
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen,
#' it is advisable to provide a
#' PDF file name and plot in a file (much faster).
#'
#' Extreme values (top and bottom) can be replaced by global quantiles at
#' \code{cut.p} and \code{1-cut.p}
#' to avoid color scales shrunk by a few outliers.
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap}
#' package to propose a simple way
#' of representing signature scores. If more advanced features are needed or
#' more graphic parameters
#' should be controlled, users should implement their own function.
#' @export
#' @examples
#' \dontrun{
#' }
#' @import ComplexHeatmap
#' @import circlize
simpleHeatmap <- function(mat.c, width, height, file.name=NULL, dend.row=NULL,
                          dend.spl=NULL, cols=NULL, pointsize=4,
                          bottom.annotation=NULL, n.col.clust=0, n.row.clust=0,
                          gap.size=0.5, cut.p=0.01, row.names=TRUE){

    if (!requireNamespace("ComplexHeatmap",quietly=TRUE))
        stop(paste0("Package \"ComplexHeatmap\" needed for this function ",
                    "to work. Please install it."))
    if (cut.p<0 || cut.p>0.1)
        stop("cut.p must lie in [0;0.1]")

    if (is.null(cols)){
        if (!requireNamespace("circlize",quietly=TRUE))
            stop(paste0("Package \"circlize\" needed for this function to ",
                        "work (generation of color scale). Please install it."))
        cols <- colorRamp2(breaks=c(min(mat.c), 0, max(mat.c)),
                           colors=c("royalblue3","white","orange"))
    }

    if (cut.p!=0)
        mat.c <- .cutExtremeValues(mat.c, cut.p)
    if (is.null(dend.spl)){
        di.spl <- stats::dist(t(mat.c))
        hc.spl <- stats::hclust(di.spl, method="ward.D")
        dend.spl <- stats::as.dendrogram(hc.spl)
    }
    if (is.null(dend.row)){
        di.gene <- stats::dist(mat.c)
        hc.gene <- stats::hclust(di.gene, method="ward.D")
        dend.row <- stats::as.dendrogram(hc.gene)
    }

    if (!is.null(file.name))
        grDevices::pdf(file.name, width=width, height=height,
                       pointsize=pointsize, useDingbats=FALSE)
    if (n.row.clust>0)
        if (n.col.clust>0)
            print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl, col=cols, show_row_names=row.names,
                show_column_names=FALSE, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                split=n.row.clust, gap=grid::unit(gap.size,"mm"),
                column_split=n.col.clust, column_gap=grid::unit(gap.size,"mm")))
    else
        print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl,col=cols,show_row_names=row.names,
                show_column_names=FALSE, use_raster=TRUE,raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                split=n.row.clust,gap=grid::unit(gap.size,"mm")))
    else
        if (n.col.clust)
            print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl, col=cols, show_row_names=row.names,
                show_column_names=FALSE, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                column_split=n.col.clust,column_gap=grid::unit(gap.size,"mm")))
    else
        print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl,col=cols, show_row_names=row.names,
                show_column_names=FALSE, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE,bottom_annotation=bottom.annotation))
    if (!is.null(file.name))
        grDevices::dev.off()

} # simpleHeatmap


#' Heatmap function for LR scores and additional data
#'
#' Generate a stack of heatmaps. The top heatmap represents ligand-receptor
#' gene signature scores,
#' while the bottom heatmap represents the second, user-chosen data set to
#' feature along with gene signatures.
#'
#' @param mat.c         A matrix with the signature scores such as output by
#' \code{\link{scoreLRGeneSignatures}}.
#' @param mat.e         A second matrix to feature below mat.c.
#' @param file.name     A PDF file name.
#' @param dend.row       A precomputed row dendrogram for \code{mat.c}.
#' @param dend.e         A precomputed row dendrogram for \code{mat.e}.
#' @param dend.spl       A precompute sample (column) dendrogram.
#' @param cols           A vector of colors to use for the heatmap of
#' \code{mat.c}.
#' @param cols.e           A vector of colors to use for the heatmap of
#' \code{mat.e}.
#' @param width          PDF width.
#' @param height         PDF height.
#' @param pointsize      PDF pointsize.
#' @param cut.p          Proportion of top and bottom values for thresholding.
#' @param vert.p         Proportion of total height dedicated to the top heatmap.
#' @param row.names.1      A logical to turn on/off the display of row (gene)
#' names in the top heatmap.
#' @param row.names.2      A logical to turn on/off the display of row (gene)
#' names in the bottom heatmap.
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen, it
#' is advisable to provide a
#' PDF file name and plot in a file (much faster).
#'
#' Extreme values (top and bottom) can be replaced by global quantiles at
#' \code{cut.p} and \code{1-cut.p}
#' to avoid color scales shrunk by a few outliers.
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap}
#' package to propose a simple way
#' of representing signature scores. If more advanced features are needed or
#' more graphic parameters
#' should be controlled, users should implement their own function.
#' @export
#' @examples
#' \dontrun{
#' }
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom ComplexHeatmap %v%
dualHeatmap <- function(mat.c, mat.e, width, height, file.name=NULL,
            dend.row=NULL, dend.spl=NULL, dend.e=NULL, cols=NULL, cols.e=NULL,
            pointsize=4, cut.p=0.01, vert.p=0.9, row.names.1=TRUE,
            row.names.2=TRUE){

    if (!requireNamespace("ComplexHeatmap",quietly=TRUE))
        stop(paste0("Package \"ComplexHeatmap\" needed for this function ",
                    "to work. Please install it."))
    if (cut.p<0 || cut.p>0.1)
        stop("cut.p must lie in [0;0.1]")
    if (vert.p<0.05 || vert.p>0.95)
        stop("vert.p must lie in [0.05;0.95]")

    if (is.null(cols)){
        if (!requireNamespace("circlize",quietly=TRUE))
            stop(paste0("Package \"circlize\" needed for this function to ",
                        "work (generation of color scale). Please install it."))
        cols <- colorRamp2(breaks=c(min(mat.c), 0, max(mat.c)),
                           colors=c("royalblue3","white","orange"))
    }
    if (is.null(cols.e)){
        if (!requireNamespace("circlize",quietly=TRUE))
            stop(paste0("Package \"circlize\" needed for this function to ",
                        "work (generation of color scale). Please install it."))
        cols.e <- colorRamp2(breaks=c(min(mat.e), 0, max(mat.e)),
                           colors=c("deepskyblue","white","tomato"))
    }
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

    hm.LR <- ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl, col=cols, show_row_names=row.names.1,
                show_column_names=FALSE, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, height=vert.p*height)
    hm.e <- ComplexHeatmap::Heatmap(mat.e, cluster_rows=dend.e,
                cluster_columns=dend.spl, col=cols.e,
                show_row_names=row.names.2, show_column_names=FALSE,
                use_raster=TRUE, raster_device="png", raster_quality=8,
                raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, height=(1-vert.p)*height)

    if (!is.null(file.name))
        grDevices::pdf(file.name, width=width, height=height,
                       pointsize=pointsize, useDingbats=FALSE)
    # import::from(ComplexHeatmap, "%v%")
    ComplexHeatmap::draw(hm.LR %v% hm.e, gap=grid::unit(1,"mm"))
    if (!is.null(file.name))
        grDevices::dev.off()

} # dualHeatmap


#' Generic gene signature scoring
#'
#' Scores generic gene signatures over the samples of a BSRDataModel object.
#'
#' @param ds         A BSRDataModel object.
#' @param ref.signatures           Gene signatures.
#' @param robust       A logical indicating that z-scores should be computed
#' with median and MAD instead of mean and standard deviation.
#' @details This function relies on a simple average of gene z-scores
#' over each signature. It is no replacement for mode advanced methods
#' such as CYBERSORT or BisqueRNA. It is provided for convenience.
#' @return A matrix containing the scores of each gene signature in each sample.
#' Note that ligand-receptor gene
#' signature scores should be computed with \code{\link{scoreLRGeneSignatures}}
#' instead.
#' @export
#' @examples
#' \dontrun{
#' }
scoreSignatures <- function(ds, ref.signatures, robust=FALSE){

    if (!is(ds, "BSRDataModel"))
        stop("ds must be BSRDataModel object")

    ref.signatures <- ref.signatures[ref.signatures$gene %in%
                                         rownames(ncounts(ds)),]
    pool <- unique(ref.signatures$gene)

    # convert normalized counts into z-scores
    if (logTransformed(ds))
        ncounts <- 2**ncounts(ds)[pool,]
    else
        ncounts <- ncounts(ds)
    if (robust)
        z <- (ncounts-apply(ncounts,1,stats::median))/apply(ncounts,1,stats::mad)
    else
        z <- (ncounts-rowMeans(ncounts))/apply(ncounts,1,stats::sd)

    # compute the gene signature scores
    pop <- unique(ref.signatures$signature)
    sig <- matrix(0, nrow=length(pop), ncol=ncol(ncounts),
                  dimnames=list(pop, colnames(ncounts)))
    for (p in pop){
        n <- sum(ref.signatures$signature==p)
        if (n > 1)
            sig[p,] <- as.numeric(colSums(
                z[ref.signatures$gene[ref.signatures$signature==p], ])/n)
        else
            sig[p,] <- as.numeric(
                z[ref.signatures$gene[ref.signatures$signature==p],])
    }

    sig

} # scoreSignatures


#' Chord Diagram of LR interactions with correlations
#'
#' @description By default, chord diagrams will be plot on disk.
#'
#' @param interactions Dataframe with Ligand, Receptor and LR.corr
#  This is reformated from a bsrinf object in order to fit circlize
#' format for plotting and let the user specify interactions he
#' wants to highlight.
#' @param path Path where to plot
#' @param filename Filename for the plot
#' @param pw.id.filter One Pathway ID accepted only to 
#  retrieve respective LR interactions.
#' @param ligands Vector of respective ligands 
#' for the LR pair that you want to 
#' highlight in the chord diagram. 
#' @param receptors Vector of respective receptors
#' for the LR pair that you want to highlight
#' in the chord diagram. 
#' @param limit Number of interactions you can visualize.
#  Maximum set to 30.
#' @param format png / svg / pdf. By default means it will 
#' plot in a pdf format.
#' @return NULL
#' @import ComplexHeatmap
#' @import circlize  
#'
#' @export
#' @examples
#' print('chord.diagram.LR')
chord.diagram.LR  <- function(interactions,path="./",filename="chord",
    pw.id.filter="R-RSA-17821",ligands=c("L1"),receptors=c("R1"),
    limit=20,format="pdf") {

    print("chord.diagram.LR")

    if (limit >= 30){
       cat("Number of selected interactions is too large",limit,".\n")
       stop("Number of visualised interactions sould be less than 30.\n")
    }

    pair.to.highlight <- paste(ligands,receptors,sep='-')

    dataframe.bsrinf<- data.frame(
        ligands=unlist(ligands(bsrinf)),
        receptors=unlist(receptors(bsrinf)),
        corr = LRinter(bsrinf)$LR.corr,
        pw.id=LRinter(bsrinf)$pw.id,
        pathways=LRinter(bsrinf)$pw.name,
        pval=LRinter(bsrinf)$pval,
        qval=LRinter(bsrinf)$qval
        )
    dataframe.bsrinf$pair <- paste(dataframe.bsrinf$ligands,dataframe.bsrinf$receptors,sep="-")
    dataframe.bsrinf <- dataframe.bsrinf[dataframe.bsrinf$pw.id==pw.id.filter,]
    if(dim(dataframe.bsrinf)[1]==0)
        stop("ID was not found")

    if (dim(dataframe.bsrinf)[1] < limit) {
        limit <- dim(dataframe.bsrinf)[1] 
        cat("Limit is too high. Only ", limit, " interactions are found.\n")
    }

    dataframe.bsrinf <- dataframe.bsrinf[order(dataframe.bsrinf$qval),]
    dataframe.bsrinf <- dataframe.bsrinf[1:limit,]

    # FROM -> LIGANDS -> GREEN
    # TO -> RECEPTORS -> RED

    #cr <- colorRampPalette(c("#FF9900","#CC0000"))(max(interactions$value)-min(interactions$value)+1)
    #cr <- colorRamp2(c(min(dataframe.bsrinf$corr), max(dataframe.bsrinf$corr)), c("lightblue","#FF9900"))
    cr <- colorRamp2(c(min(dataframe.bsrinf$corr), max(dataframe.bsrinf$corr)), c("white","#febd17"))

    #grid.col <- c(receptors = "red", ligands = "green")

    myList.ligands <- rep("gray25",times=length(dataframe.bsrinf$ligands))
    names(myList.ligands)  <- as.list(dataframe.bsrinf$ligands)
    
    myList.receptors <- rep("#7fbb00",times=length(dataframe.bsrinf$receptors))
    names(myList.receptors)  <- as.list(dataframe.bsrinf$receptors)

    myList <- c(myList.receptors,myList.ligands)

    link.col <- rep("dodgerblue3",nrow(dataframe.bsrinf))
    
    link.lwd <- rep(1,nrow(dataframe.bsrinf))
    
    link.width <- rep(0.12,nrow(dataframe.bsrinf))

    dataframe.bsrinf <- dataframe.bsrinf[order(1:limit),]

    index.filter <- which(pair.to.highlight %in% dataframe.bsrinf$pair)
    print(index.filter)

    link.col[index.filter] <- "#d40000"
    link.lwd[index.filter] <- 3
    link.width[index.filter] <- 0.15

    if (format=="svg")#3.6
        svg(paste0(path,"/",filename,".svg"),width=4, height=4) #inch
    if (format=="png")
        png(paste0(path,"/",filename,".png"))
    if (format=="pdf")
        pdf(paste0(path,"/",filename,".pdf"))

    interactions <- data.frame(from=dataframe.bsrinf$ligands,to=dataframe.bsrinf$receptors,value=dataframe.bsrinf$corr)
    chordDiagramFromDataFrame(interactions,
          #grid.col = grid.col, 
          col=cr,
          annotationTrack = "grid", 
          grid.col = myList, 
          transparency = 0.7,
          preAllocateTracks = 1,   
          directional = 1,
          direction.type = "arrows",
          link.arr.length = link.width,
          link.arr.width  = link.width,
          link.arr.type  = "triangle", 
          link.arr.lty = par("lty"),
          link.arr.lwd = link.lwd, 
          link.arr.col = link.col,
          big.gap = 2, 
          small.gap = 1)
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  circos.text(mean(xlim), ylim[1] + 1.9, sector.name, 
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.7)

  circos.axis(h="top",labels=FALSE,minor.ticks=FALSE,
                    major.tick.length=1,
                    major.at=c(xlim), sector.index=sector.name,
                    track.index=2)
  })

  # LEGEND 

    lgd_points = Legend(labels = c("Ligands", "Receptors"), type = "points", pch = 16,
        legend_gp = gpar(col = c("gray25","#7fbb00")), title_position = "topleft", 
            labels_gp = gpar( font = 6),
        title = "LR") 
   
    lgd_links = Legend(at = c(round(min(interactions$value),digits = 2) , round(max(interactions$value),digits = 2)), col_fun = cr, 
         title = "Correlation",direction ="horizontal"   ,   
         grid_width = unit(0.9, "mm") ,     grid_height = unit(1.3, "mm") , 

        labels_gp = gpar( font = 6),title_position = "topcenter",
       )# border = "black"

    lgd_list_vertical = packLegend(lgd_points)

     draw(lgd_list_vertical, x = unit(2, "mm"), y = unit(2, "mm"), just = c("left", "bottom"))
    draw(lgd_links, x = unit(2.7, "inch"), y = unit(2, "mm"),just = c("left", "bottom"))

     circos.clear()

     dev.off()

 NULL   
}

