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

#' Bubble Plot to explore LR & Pathways 
#'
#' Quick check to observe LR - Pathways association
#' with their respective correlation
#' and Q-values.
#'
#' @param bsrinf     BulkSignalR inference object.
#' @param pathways  Vector or pathway names
#' @param threshold     BulkSignalR data model object.
#' @param filter.L     Filter a list of ligands.
#' @param filter.R     Filter a list of receptors.
#' @param path      Path directory to plot.
#' @param filename     An output filename.
#' @param color     Main color used for the gradient.
#' @param image_width     Global image width size.
#' @param image_height    Global image height size.
#' @param format   File format svg (default)/png/pdf

#' @return  A plot is created with the given filename. 
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap}
#' package to propose a simple way
#' of representing LR - Pathways association
#' with their respective correlation
#' and Q-values.
#' @export
#' @examples
#' \dontrun{
#' }
#' @import ggplot2
bubblePlot.PathwaysLR <- function(bsrinf,
    pathways=c("Cell surface interactions at the vascular wall"),
    threshold = 0.001,
    filter.L = NULL, 
    filter.R = NULL,
    path=".",
    filename="DotPlot",
    color= "#16a647",
    image_width  = 16, 
    image_height = 7,
    format="svg" ) {

    if (!dir.exists(path)) {
        stop("Directory is invalid.")
    }

    filtered.brinf <- as.data.frame(LRinter(bsrinf))
    filtered.brinf <- filtered.brinf[filtered.brinf$qval < threshold,]

    if(!is.null(filter.R) | !is.null(filter.L))
        filtered.brinf <- filtered.brinf[filtered.brinf$L %in% filter.L | filtered.brinf$R %in% filter.R ,]

    filtered.brinf$LR <- paste (filtered.brinf$L,filtered.brinf$R,sep= "->")

    filtered.brinf <- filtered.brinf[,c("LR","pw.name","LR.corr","qval")] 
    filtered.brinf$log10 <- -log10(filtered.brinf$qval)

    filtered.brinf <- filtered.brinf[filtered.brinf$pw.name %in% pathways,]
   
    if(dim(filtered.brinf)[1]==0)
        stop("The pathways you have selected do no exist.")

    limit.P <- 5
    if(length(pathways)>=limit.P){
        cat("We recommend less than ",limit.P," pathways.","\n")
        stop("Too many pathways were given in input.")
    }
    limit.LR <- 50
    if(length(unique(filtered.brinf$LR))>=limit.LR){
        cat("Too many LR interactions detected (",length(unique(filtered.brinf$LR)) ,").","\n")
        cat("We recommend less than ", limit.LR, " LR interactions to visualize.","\n")
        stop("Try to reduce (Qval-Threshold, number of pathways...).\n")
    }

    cat(length(unique(filtered.brinf$LR))," LR interactions detected.\n")

    # Adjust image
    width <- ( length(unique(filtered.brinf$LR)) * 0.5 ) / 2.54
    if(width > image_width ) image_width <- width
    height <- (length(pathways) * 0.5 ) / 2.54
    if(height  > image_height ) image_height  <- height

    print(image_height)
    print(image_width)

    format <- match.arg(format)

    if (format=="svg")
        svg(file=paste0(path,filename,".svg")
         ,width=(image_width+3)/2.54, height=image_height/2.54) 

    if (format=="png")
        png(file=paste0(path,filename,".png")
         ,width=(image_width+3)/2.54, height=image_height/2.54) 

    if (format=="pdf")
       pdf(file=paste0(path,filename,".pdf")
            ,width=(image_width+3)/2.54, height=image_height/2.54) 


    print(ggplot2::ggplot(filtered.brinf, aes(x = LR, y = pw.name)) + 
      ggplot2::geom_point(aes(size = -log10(LR.corr), fill = log10 ), alpha = 0.75, shape = 21) + 
       labs(x= "", y = "", size = "-log10 (LR.corr)", fill = "-log10 (Qval)")  + 
       ggplot2::theme(legend.key=ggplot2::element_blank(),  
       legend.key.size = unit(0.2, "cm"),
       legend.position = "right",  legend.box = "horizontal",
       axis.text.x = ggplot2::element_text(colour = "black", size = 6, angle = 90, vjust = 0.3, hjust = 1), 
       axis.text.y = ggplot2::element_text(colour = "black", face = "bold", size = 6), 
       legend.text = ggplot2::element_text(size = 5, face ="bold", colour ="black"), 
       legend.title = ggplot2::element_text(size = 6, face = "bold"), 
       panel.background = ggplot2::element_blank(), 
       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1.2), 
       panel.grid.major.y = ggplot2::element_line(colour = "grey95")) +    
       ggplot2::scale_fill_gradient(low = "white",high = color,space = "Lab",na.value = "grey50",
       guide = "colourbar",aesthetics = "fill") +
       ggplot2::scale_y_discrete(limits = rev(levels(filtered.brinf$pw.name))) )
   
   grDevices::dev.off()

}


#' Heatmap function for gene expression of signature
#'
#' Generate one heatmap re-used by 
#' by  \code{\link{scoreLRGeneSignatures}}
#'
#' @param counts  Matrice of counts exported from a
#  BulksignalR data model object.
#' @param width   Heatmap  width size.
#' @param height      Heatmap  height size.
#' @param image_width     Global image width size.
#' @param image_height    Global image height size.
#' @param scoring Vector of scored sample for a 
#' a previously choosen pathway.
#' @param cols.scoring   Fixed colorRamp2 object.
#' @param palette   Palette from 
#' HCL colormaps supported by ComplexHeatmap.
#' @param show_column_names   Add column names on heatmap.
#'
#' @return ComplexHeatmap object. 
#'
#' This is a convenience function that relies 
#' on the \code{ComplexHeatmap}
#' package to propose a simple way
#' of representing expression of genes involved in a specific
#' pathway.
#' @export
#' @examples
#' \dontrun{
#' }
#' @import ComplexHeatmap
#' @import circlize   
customheatmap <- function(counts, 
    width=5, 
    height=10 ,
    scoring = c(-1.5,0,4,5,6.1,0.3) ,
    cols.scoring,
    palette= "Blues 3",
    show_column_names = FALSE
    ) {

    print("customheatmap")

    counts <- data.matrix(counts)
    counts.scaled = t(scale(t(counts)))

    cols <- circlize::colorRamp2(breaks=c(-1, 0, 3), hcl_palette = palette,reverse=TRUE) 

    top.annotation <- HeatmapAnnotation(    
     border = c(scoring = TRUE),
     show_legend = FALSE,
     simple_anno_size = unit(2.5, "mm"),
     show_annotation_name = FALSE,
          scoring  = as.vector(scoring),
          col = list( scoring =   cols.scoring ) 
    )

    di.gene <- stats::dist(counts.scaled)
    hc.gene <- stats::hclust(di.gene, method = "ward.D")
    dend.row <- stats::as.dendrogram(hc.gene)

    di.spl <- stats::dist(t(counts.scaled))
    hc.spl <- stats::hclust(di.spl, method="ward.D")
    dend.spl <- stats::as.dendrogram(hc.spl)

     ComplexHeatmap::ht_opt(heatmap_border = TRUE,
          annotation_border = FALSE
     ) 

     p <- ComplexHeatmap::Heatmap(counts.scaled, 
       cluster_rows = dend.row,   cluster_columns = dend.spl,
       show_row_dend = FALSE,  show_column_dend = TRUE,
       col = cols, show_row_names = TRUE, show_column_names = show_column_names, 
       use_raster = TRUE,raster_device = "png", raster_quality = 8,  raster_by_magick = FALSE,
       rect_gp = gpar(col= "white"),  row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 4),
       top_annotation = top.annotation,
       show_heatmap_legend = FALSE, 
       width =unit(width, "cm"),
       height =unit(height, "cm"),# all heatmap components (excluding the legends)
       column_gap = grid::unit(0.5, "mm")
     )

    return(p)
} # customheatmap

#' Heatmap function for gene expression of signature
#'
#' Generate a list of heatmaps for ligand,
#' receptor and target genes
#' for a specific pathway
#'
#' @param pathway        Pathway name
#' @param bsrdm     BulkSignalR data model object.
#' @param bsrsig     BulkSignalR signature object.
#' @param path     Absolute path to write your ouput.
#' @param filename     An output filename.
#' @param image_width     Global image width size.
#' @param image_height    Global image height size.
#' @param format   File format svg (default)/png/pdf
#' @param show_column_names   Add column names on heatmap.

#' @return  A plot is created with the given filename. 
#'
#' This is a convenience function that relies on the \code{ComplexHeatmap}
#' package to propose a simple way
#' of representing expression of genes involved in a specific
#' pathway.
#' @export
#' @examples
#' \dontrun{
#' }
#' @import ComplexHeatmap
#' @import circlize
signatureHeatmaps <- function(
        pathway="Cell surface interactions at the vascular wall",
        bsrdm=bsrdm,
        bsrsig=bsrsig.redPP,
        path=".",
        filename="ExpresssionSignature",
        image_width  = 15, 
        image_height = 8,
        format=c("svg","png","pdf"),
        show_column_names= FALSE
        ){

    print("signatureHeatmaps")

    #check path
    if (!dir.exists(path)) {
        stop("Directory is invalid.")
    }

    idx.path.sig    <- which(pathways(bsrsig) == pathway)

    if(rlang::is_empty(idx.path.sig)){
        stop("Pathway is not defined in signature.")
    }

    scoresPathway <- scoreLRGeneSignatures(bsrdm,bsrsig,name.by.pathway=TRUE,rownames.LRP=FALSE)

    counts <- as.data.frame(bsrdm@ncounts)

    filter.L      <- ligands(bsrsig)[[idx.path.sig]]

    counts.L <- counts[filter.L,] 
    palette.L     <- "RdPu"
    cols.L <- circlize::colorRamp2(breaks=c(-1, 0, 1), hcl_palette = palette.L,reverse=TRUE) 

    filter.R      <- receptors(bsrsig)[[idx.path.sig]]
    filter.T      <- tGenes(bsrsig)[[idx.path.sig]]

    # Remove in receptors, genes that are potential targets.
    filter.R <- filter.R[! filter.R %in% filter.T]

    counts.R <- counts[filter.R,] 
    palette.R     <- "YlGn"
    cols.R <- circlize::colorRamp2(breaks=c(-1, 0, 1), hcl_palette = palette.R,reverse=TRUE) 

    counts.T <- counts[filter.T,] 
    palette.T     <- "Blues 3"
    cols.T <- circlize::colorRamp2(breaks=c(-1, 0, 1), hcl_palette = palette.T,reverse=TRUE) 

    cols.scoring  <- circlize::colorRamp2(breaks=c(-1,0,1),colors=c("blue","white","red"))

    abundance.samples <- dim(counts.T)[2]
    abundance.genes <- c(dim(counts.T)[1],dim(counts.R)[1],dim(counts.L)[1])
  
    # Adjust size image / heatmaps
    heatmap_ind_height <- (max(abundance.genes)*0.5)/2.54

    heatmap_ind_witdth <- (abundance.samples*0.5)/2.54

    width <-  (heatmap_ind_witdth+3) * 3
    height <- heatmap_ind_height + 6

    if(width > image_width ) image_width <- width
    if(height  > image_height ) image_height  <- height

    # width and height in cm
    p.T <- customheatmap(counts=counts.T,
     width=heatmap_ind_witdth,height=heatmap_ind_height, 
         scoring=as.vector(scoresPathway[idx.path.sig,]),
      palette= palette.T,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    p.R <- customheatmap(counts=counts.R, 
        width=heatmap_ind_witdth,height=heatmap_ind_height ,
         scoring=as.vector(scoresPathway[idx.path.sig,]),
          palette= palette.R,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    p.L <- customheatmap(counts=counts.L, 
        width= heatmap_ind_witdth,height=heatmap_ind_height,
         scoring=as.vector(scoresPathway[idx.path.sig,]),
          palette= palette.L,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    format <- match.arg(format)

    if (format=="svg")
        svg(file=paste0(path,filename,".svg")
         ,width=(image_width+6)/2.54, height=image_height/2.54) 

    if (format=="png")
        png(file=paste0(path,filename,".png")
         ,width=(image_width+6)/2.54, height=image_height/2.54) 

    if (format=="pdf")
       pdf(file=paste0(path,filename,".pdf")
            ,width=(image_width+6)/2.54, height=image_height/2.54) 

     grid.newpage()#grid
     pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))
     
     pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
     ComplexHeatmap::draw(p.L, newpage = FALSE)
     upViewport()

     pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
     ComplexHeatmap::draw(p.R, newpage = FALSE)
     upViewport()

     pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
     ComplexHeatmap::draw(p.T, newpage = FALSE)
     upViewport()

    lgd_score = ComplexHeatmap::Legend( col_fun = cols.scoring, direction = "horizontal",  
    title = "Gene signature scores",title_position = "topleft",
    grid_height = unit(0.2, "mm"),
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6,fontface="bold")
    )

    lgd_heatmap.T = ComplexHeatmap::Legend( col_fun = cols.T, direction = "horizontal",  
    title = "Expression target",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = gpar( fontsize = 6),
    title_gp = gpar( fontsize = 6,fontface="bold")
    )

    lgd_heatmap.R = ComplexHeatmap::Legend( col_fun = cols.R, direction = "horizontal",  
    title = "Expression receptor",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = gpar( fontsize = 6),
    title_gp = gpar( fontsize = 6,fontface="bold")
    )

    lgd_heatmap.L = ComplexHeatmap::Legend( col_fun = cols.L, direction = "horizontal",  
    title = "Expression ligand",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = gpar( fontsize = 6),
    title_gp = gpar( fontsize = 6,fontface="bold")
    )

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
    ComplexHeatmap::draw(lgd_score, y = unit( (image_height+2)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.L,y = unit( (image_height+5)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.R, y = unit( (image_height+7)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.T,y = unit( (image_height+9)/2, "cm"))

    upViewport()

    grDevices::dev.off()

} # signatureHeatmaps


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
#' @param row.names      A logical to turn on/off the display 
#' of row (gene) names.
#' @param column.names      A logical to turn on/off the display
#' of column (sample) names.
#' 
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
                          gap.size=0.5, cut.p=0.01, row.names=TRUE,
                          column.names = TRUE 
                          ){

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
                show_column_names=column.names, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                split=n.row.clust, gap=grid::unit(gap.size,"mm"),
                column_split=n.col.clust, column_gap=grid::unit(gap.size,"mm")))
    else
        print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl,col=cols,show_row_names=row.names,
                show_column_names=column.names, use_raster=TRUE,raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                split=n.row.clust,gap=grid::unit(gap.size,"mm")))
    else
        if (n.col.clust)
            print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl, col=cols, show_row_names=row.names,
                show_column_names=column.names, use_raster=TRUE, raster_device="png",
                raster_quality=8, raster_by_magick=FALSE,
                row_names_gp=grid::gpar(fontsize=pointsize),
                show_row_dend=TRUE, bottom_annotation=bottom.annotation,
                column_split=n.col.clust,column_gap=grid::unit(gap.size,"mm")))
    else
        print(ComplexHeatmap::Heatmap(mat.c, cluster_rows=dend.row,
                cluster_columns=dend.spl,col=cols, show_row_names=row.names,
                show_column_names=column.names, use_raster=TRUE, raster_device="png",
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

#' Alluvial plot (sankey like)
#'
#' @description Representation of the links
#' between Ligands,Receptors and Pathways.
#'
#' @param bsrsig object bsrinf inference.
#' @param keywords vector of pathways.
#' @param qval threshold over Q-value.
#' @param type filter on Ligand, Recptor or pathway id.
#' @return ggplot object to print
#' This is a convenience function that relies on the \code{ggalluvial}
#' package to propose a simple way
#' of representing Ligands, Receptors 
#  and underlying Pathways associated. 
#' @import ggplot2
#' @import ggalluvial
#' @export
#' @examples
#' print('ggalluvial')
alluvial.plot <- function(bsrinf,keywords=c("COL4A1"),type=c("L","R","pw.id"),qval=1) {

     interactions <- data.frame(
           L =  unlist(ligands(bsrinf)) 
          ,R =  unlist(receptors(bsrinf))
          ,pw.name = LRinter(bsrinf)$pw.name
          ,pw.id = LRinter(bsrinf)$pw.id
          ,qval = LRinter(bsrinf)$qval
          ,targets = sapply(tGenes(bsrinf), paste, collapse = "  ")

          )

    type <- match.arg(type)

     if(type=="L")
          subset.interactions                <- interactions[interactions$L %in% keywords,]
     if(type=="R")
          subset.interactions                <- interactions[interactions$R %in% keywords,]
     if(type=="pw.id")
          subset.interactions                <- interactions[interactions$pw.id %in% keywords,]  

       if (dim(subset.interactions)[1]==0){
        cat(paste(keywords, collapse = ' ')," for ", type, " not found.","\n", sep="")
        stop("Try another value for filtering.")
    }
     subset.interactions                <- subset.interactions[subset.interactions$qval < qval,]  
     subset.interactions$count <- 1
     subset.interactions       <- subset.interactions[,c("L" ,   "R"  ,  "pw.name" ,"count")]

     stratum <- ggalluvial::StatStratum

     pl <- ggplot2::ggplot(as.data.frame(subset.interactions),
        aes(y = count, axis1 = L, axis2 = R,axis3 = pw.name)) +
        ggalluvial::geom_alluvium(aes(fill = R), width = 1/12) +
        ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
        ggplot2::geom_label(stat = stratum,  aes(label = ggplot2::after_stat(stratum)))+
        ggplot2::scale_x_discrete(limits = c("L", "R","pw.name"), expand = c(0.5, 0.5)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::ggtitle("Ligand-Receptor Interactions & Underlying Pathways")

     pl <- pl + ggplot2::theme_bw()
     pl <- pl + ggplot2::theme(legend.position = "none")
     pl <- pl + ggplot2::theme(axis.title = ggplot2::element_blank()
                       ,text = ggplot2::element_text(size = 12)
                       , plot.title = ggplot2::element_text(hjust = .5)
                       , axis.text.y = ggplot2::element_blank()
                       , axis.text.x = ggplot2::element_blank()
                       , axis.ticks = ggplot2::element_blank()  
                       , panel.grid = ggplot2::element_blank())

 return (pl)
} #alluvial.plot


#' Chord Diagram of LR interactions with correlations
#'
#' @description By default, chord diagrams will be plot on disk.
#'
#' @param bsrinf bsrinf object 
#' @param path Path where to plot
#' @param filename Filename for the plot
#' @param pw.id.filter One Pathway ID accepted only to 
#  retrieve respective LR interactions.
#' @param ligand Ligand
#' for the LR pair that you want to 
#' highlight in the chord diagram. 
#' @param receptor Receptor
#' for the LR pair that you want to highlight
#' in the chord diagram. 
#' @param limit Number of interactions you can visualize.
#  Maximum set to 30.
#' @param format png / svg / pdf. By default means it will 
#' plot in a pdf format.
#' @param width width of image 
#' @param height height of image 
#' @return NULL
#' @import ComplexHeatmap
#' @import circlize  
#'
#' @export
#' @examples
#' print('chord.diagram.LR')
chord.diagram.LR  <- function(bsrinf,path="./",filename="chord",
    pw.id.filter="R-RSA-17821",ligand="L1",receptor="R1",
    limit=20,format=c("svg","png","pdf"),width=4,height=4) {

    format <- match.arg(format)

    print("chord.diagram.LR")

    if (limit >= 30){
       cat("Number of selected interactions is too large",limit,".\n")
       stop("Number of visualised interactions sould be less than 30.\n")
    }

    pair.to.highlight <- paste(ligand,receptor,sep='-')

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
    if (!pair.to.highlight %in% dataframe.bsrinf$pair) {
        cat("Highlighted LR pair ", pair.to.highlight, " was not found for",pw.id.filter,".\n")
    }

    dataframe.bsrinf <- dataframe.bsrinf[dataframe.bsrinf$pw.id==pw.id.filter,]

    if(dim(dataframe.bsrinf)[1]==0)
        stop("Pathway ID was not found")

    if (dim(dataframe.bsrinf)[1] < limit) {
        limit <- dim(dataframe.bsrinf)[1] 
        cat("Only ", limit, " interactions were found.\n")
    }

    dataframe.bsrinf <- dataframe.bsrinf[order(dataframe.bsrinf$qval),]
    dataframe.bsrinf <- dataframe.bsrinf[1:limit,]

    cr <- colorRamp2(c(min(dataframe.bsrinf$corr), max(dataframe.bsrinf$corr)), c("white","#febd17"))

    myList.ligands <- rep("gray25",times=length(dataframe.bsrinf$ligands))
    names(myList.ligands)  <- as.list(dataframe.bsrinf$ligands)
    
    myList.receptors <- rep("#7fbb00",times=length(dataframe.bsrinf$receptors))
    names(myList.receptors)  <- as.list(dataframe.bsrinf$receptors)

    myList <- c(myList.receptors,myList.ligands)

    link.col <- rep("dodgerblue3",nrow(dataframe.bsrinf))
    
    link.lwd <- rep(1,nrow(dataframe.bsrinf))
    
    link.width <- rep(0.12,nrow(dataframe.bsrinf))

    index.filter <- which(pair.to.highlight == dataframe.bsrinf$pair)

    link.col[index.filter] <- "#d40000"
    link.lwd[index.filter] <- 3
    link.width[index.filter] <- 0.15

    if (format=="svg")#3.6
        svg(paste0(path,"/",filename,".svg"),width=width, height=height) #inch
    if (format=="png")
        png(paste0(path,"/",filename,".png"),width=width, height=height)
    if (format=="pdf")
        pdf(paste0(path,"/",filename,".pdf"),width=width, height=height)

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

