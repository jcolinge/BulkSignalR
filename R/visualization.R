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
#' @param width     Global image width size.
#' @param height    Global image height size.
#' @param pointsize Global pointsize.
#' @param format   File format svg (default),png or pdf.
#'
#' @return  A plot is created with the given filename. 
#'
#' This is a convenience function to propose a simple way
#' of representing LR - Pathways association
#' with their respective correlation
#' and Q-values.
#'
#' @export
#' @examples
#' print('bubblePlotPathwaysLR')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' pathways <- c("PD-1 signaling","Interferon gamma signaling")
#' bubblePlotPathwaysLR(bsrinf,
#'    pathways = pathways, 
#'    threshold = 1,
#'    path = "./",
#'    color = "red",
#'    filename  = "sdc_bubble", 
#'    width  = 16, 
#'    height = 7,
#'    pointsize = 8
#'    )  
#' @import ggplot2
bubblePlotPathwaysLR <- function(bsrinf,
    pathways=c("Cell surface interactions at the vascular wall"),
    threshold = 1,
    filter.L=NULL, 
    filter.R=NULL,
    path="./",
    filename="bubblePlot",
    color="#16a647",
    width=16, 
    height=7,
    pointsize=6, 
    format=c("svg","png","pdf") ) {

    if (!dir.exists(path)) {
        stop("Directory is invalid.")
    }

    filtered.brinf <- as.data.frame(LRinter(bsrinf))
    filtered.brinf <- filtered.brinf[filtered.brinf$qval < threshold,]

    if(!is.null(filter.R) | !is.null(filter.L))
        filtered.brinf <- filtered.brinf[filtered.brinf$L %in% filter.L | filtered.brinf$R %in% filter.R ,]

    filtered.brinf$LR <- paste (filtered.brinf$L,filtered.brinf$R,sep= "->")

    filtered.brinf <- filtered.brinf[,c("LR","pw.name","LR.corr","qval")] 
    filtered.brinf$log10.qval <- -log10(filtered.brinf$qval)
    filtered.brinf$log10.LR.corr <- -log10(filtered.brinf$LR.corr)

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
    width.fit <- ( length(unique(filtered.brinf$LR)) * 0.5 ) / 2.54
    if(width.fit > width ) width <- width.fit
    height.fit <- (length(pathways) * 0.5 ) / 2.54
    if(height.fit  > height ) height  <- height.fit


    format <- match.arg(format)

    #if (!is.null(filename))
        if (format=="svg")
            grDevices::svg(file=paste0(path,filename,".svg")
             ,width=(width+3)/2.54, height=height/2.54) 

        if (format=="png")
            grDevices::png(file=paste0(path,filename,".png")
             ,width=(width+3)/2.54, height=height/2.54) 

        if (format=="pdf")
           grDevices::pdf(file=paste0(path,filename,".pdf")
                ,width=(width+3)/2.54, height=height/2.54) 

    #else grDevices::dev.new(width=(width+3)/2.54, height=height/2.54)
    #ggsave
    print(ggplot2::ggplot(filtered.brinf, 
       ggplot2::aes_string(x = "LR", y = "pw.name")) + 
       ggplot2::geom_point(ggplot2::aes_string(size = "log10.LR.corr",
       fill = "log10.qval" ), alpha = 0.75, shape = 21) + 
       labs(x= "", y = "", size = "-log10 (LR.corr)",
        fill = "-log10 (Qval)")  + 
       ggplot2::theme(legend.key=ggplot2::element_blank(),  
       legend.key.size = unit(0.2, "cm"),
       legend.position = "right",  legend.box = "horizontal",
       axis.text.x = ggplot2::element_text(colour = "black", 
       size = pointsize, angle = 90, vjust = 0.3, hjust = 1), 
       axis.text.y = ggplot2::element_text(colour = "black", 
       face = "bold", size = pointsize), 
       legend.text = ggplot2::element_text(size = pointsize-1, 
       face ="bold", colour ="black"), 
       legend.title = ggplot2::element_text(size = pointsize, 
       face = "bold"), 
       panel.background = ggplot2::element_blank(), 
       panel.border = ggplot2::element_rect(colour = "black", 
       fill = NA, size = 1.2), 
       panel.grid.major.y = ggplot2::element_line(colour = "grey95")) +    
       ggplot2::scale_fill_gradient(low = "white",
       high = color,space = "Lab",na.value = "grey50",
       guide = "colourbar",aesthetics = "fill") +
       ggplot2::scale_y_discrete(limits = rev(levels(filtered.brinf$pw.name))) )
    
    #if (!is.null(filename))
        grDevices::dev.off()

}


#' Heatmap function for gene expression of signature
#'
#' Generate one heatmap re-used by 
#' by  
#' \code{"\link[=BSRDataModel-class]{scoreLRGeneSignatures}"}
#'
#' @param counts  Matrice of counts exported from a
#  BulksignalR data model object.
#' @param width   Heatmap  width size.
#' @param height      Heatmap  height size.
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
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
.customheatmap <- function(counts, 
    width=5, 
    height=10 ,
    scoring = c(-1.5,0,4,5,6.1,0.3) ,
    cols.scoring,
    palette= "Blues 3",
    show_column_names = FALSE
    ) {

    print(".customheatmap")

    counts <- data.matrix(counts)
    counts.scaled = t(scale(t(counts)))

    cols <- circlize::colorRamp2(breaks=c(-1, 0, 3), 
        hcl_palette = palette,reverse=TRUE) 

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
       col = cols, show_row_names = TRUE, 
       show_column_names = show_column_names, 
       use_raster = TRUE,raster_device = "png", 
       raster_quality = 8,  raster_by_magick = FALSE,
       rect_gp = grid::gpar(col= "white"),  
       row_names_gp = grid::gpar(fontsize = 6),
       column_names_gp = grid::gpar(fontsize = 4),
       top_annotation = top.annotation,
       show_heatmap_legend = FALSE, 
       width =unit(width, "cm"),
       height =unit(height, "cm"),
       column_gap = grid::unit(0.5, "mm")
     )

    return(p)
} # .customheatmap

#' Heatmap function for gene expression of signature
#'
#' Generate a list of heatmaps for ligand,
#' receptor and target genes
#' for a specific pathway
#'
#' @param pathway        Pathway name
#' @param bsrdm     BulkSignalR data model object.
#' @param bsrsig     BulkSignalR signature object.
#' @param path     Path to write your ouput.
#' @param filename     An output filename.
#' @param width     Global image width size.
#' @param height    Global image height size.
#' @param format   File format svg (default)/png/pdf
#' @param show_column_names   Add column names on heatmap.

#' @return  A plot is created with the given filename. 
#'
#' This is a convenience function 
#' to propose a simple way
#' of representing expression of genes 
#' involved in a specific pathway.
#' 
#' @export
#' @examples
#' print('signatureHeatmaps')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 
#' bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP,qval.thres=0.001)
#' pathway1 <- "Elastic fibre formation"
#' signatureHeatmaps(
#'        pathway = pathway1,
#'        bsrdm = bsrdm,
#'        bsrsig = bsrsig.redPBP,
#'        path = "./",
#'        filename = "sdc_signatureheatmap",
#'        width  = 15,
#'        height = 10 ,
#'        show_column_names = TRUE)
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @import grid
signatureHeatmaps <- function(
        pathway="Cell surface interactions at the vascular wall",
        bsrdm,
        bsrsig,
        path="./",
        filename="ExpresssionSignature",
        width=15, 
        height=8,
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

    width.fit <-  (heatmap_ind_witdth+3) * 3
    height.fit <- heatmap_ind_height + 6

    if(width.fit > width ) width <- width.fit
    if(height.fit  > height ) height  <- height.fit

    # width and height in cm
    p.T <- .customheatmap(counts=counts.T,
     width=heatmap_ind_witdth,height=heatmap_ind_height, 
         scoring=as.vector(scoresPathway[idx.path.sig,]),
      palette= palette.T,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    p.R <- .customheatmap(counts=counts.R, 
        width=heatmap_ind_witdth,height=heatmap_ind_height ,
         scoring=as.vector(scoresPathway[idx.path.sig,]),
          palette= palette.R,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    p.L <- .customheatmap(counts=counts.L, 
        width= heatmap_ind_witdth,height=heatmap_ind_height,
         scoring=as.vector(scoresPathway[idx.path.sig,]),
          palette= palette.L,cols.scoring=cols.scoring,
          show_column_names = show_column_names)

    format <- match.arg(format)

    if (format=="svg")
        grDevices::svg(file=paste0(path,filename,".svg")
         ,width=(width+6)/2.54, height=height/2.54) 

    if (format=="png")
        grDevices::png(file=paste0(path,filename,".png")
         ,width=(width+6)/2.54, height=height/2.54) 

    if (format=="pdf")
       grDevices::pdf(file=paste0(path,filename,".pdf")
            ,width=(width+6)/2.54, height=height/2.54) 

     grid::grid.newpage()#grid
     grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr = 1, nc = 4)))
     
     grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
     ComplexHeatmap::draw(p.L, newpage = FALSE)
     grid::upViewport()

     grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
     ComplexHeatmap::draw(p.R, newpage = FALSE)
     grid::upViewport()

     grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
     ComplexHeatmap::draw(p.T, newpage = FALSE)
     grid::upViewport()

    lgd_score = ComplexHeatmap::Legend( col_fun = cols.scoring, direction = "horizontal",  
    title = "Gene signature scores",title_position = "topleft",
    grid_height = unit(0.2, "mm"),
    labels_gp = grid::gpar(fontsize = 6),
    title_gp = grid::gpar(fontsize = 6,fontface="bold")
    )

    lgd_heatmap.T = ComplexHeatmap::Legend( col_fun = cols.T, direction = "horizontal",  
    title = "Expression target",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = grid::gpar( fontsize = 6),
    title_gp = grid::gpar( fontsize = 6,fontface="bold")
    )

    lgd_heatmap.R = ComplexHeatmap::Legend( col_fun = cols.R, direction = "horizontal",  
    title = "Expression receptor",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = grid::gpar( fontsize = 6),
    title_gp = grid::gpar( fontsize = 6,fontface="bold")
    )

    lgd_heatmap.L = ComplexHeatmap::Legend( col_fun = cols.L, direction = "horizontal",  
    title = "Expression ligand",title_position = "topleft",
    grid_width = unit(0.7, "mm"),
    grid_height = unit(0.2, "mm") , 
    labels_gp = grid::gpar( fontsize = 6),
    title_gp = grid::gpar( fontsize = 6,fontface="bold")
    )

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 4))
    ComplexHeatmap::draw(lgd_score, y = unit( (height+2)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.L,y = unit( (height+5)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.R, y = unit( (height+7)/2, "cm"))
    ComplexHeatmap::draw(lgd_heatmap.T,y = unit( (height+9)/2, "cm"))

    grid::upViewport()

    grDevices::dev.off()

} # signatureHeatmaps


#' Heatmap function for LR scores
#'
#' Generate a PDF file with a heatmap representing ligand-receptor gene
#' signature scores.
#'
#' @param mat.c         A matrix with the signature scores such as output by
#' \code{scoreLRGeneSignatures()}.
#' @param path directory where to plot file.
#' @param filename file name.
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
#' @param column.names      A logical to turn on/off the display
#' of column (sample) names.
#' @param hcl_palette    support for HCL colormaps in ComplexHeatmap
#' using color mapping function with circlize::colorRamp2().
#' palettes are listed in grDevides::hcl.pals().
#' of row (gene) names.
#' @param reverse    A logicial to reverse or not colors in hclpalette.
#' @param format   File format svg (default)/png/pdf
#' 
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen,
#' it is advisable to provide a
#' PDF file name and plot in a file (much faster).
#'
#' If hcl_palette is set. colors parameter won't be used. 
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
#' print('simpleHeatmap')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redBP <- reduceToBestPathway(bsrinf) 
#' bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP,
#' qval.thres=0.001)
#'
#' scoresLR <- scoreLRGeneSignatures(bsrdm,bsrsig.redBP,
#'                        name.by.pathway=FALSE)
#' simpleHeatmap(scoresLR[1:20,], 
#'                   path = "./",
#'                   filename = "sdc_scoresLR",
#'                   column.names = TRUE, 
#'                   height = 5, width = 9,
#'                   pointsize = 10,
#'                   hcl_palette = "Cividis"                   
#'                   )
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
simpleHeatmap <- function(mat.c, width, height, 
                          path="./", filename="simpleHeatmap",
                          dend.row=NULL,
                          dend.spl=NULL, cols=NULL, pointsize=4,
                          bottom.annotation=NULL, n.col.clust=0, n.row.clust=0,
                          gap.size=0.5, cut.p=0.01, row.names=TRUE,
                          column.names = TRUE ,hcl_palette = NULL,
                          reverse = FALSE, format=c("svg","png","pdf")
                          ){

    if (!requireNamespace("ComplexHeatmap",quietly=TRUE))
        stop(paste0("Package \"ComplexHeatmap\" needed for this function ",
                    "to work. Please install it."))
    if (cut.p<0 || cut.p>0.1)
        stop("cut.p must lie in [0;0.1]")

    if (cut.p!=0)
        mat.c.cut <- .cutExtremeValues(mat.c, cut.p)

    if (is.null(cols)){
        if (!requireNamespace("circlize",quietly=TRUE))
            stop(paste0("Package \"circlize\" needed for this function to ",
                        "work (generation of color scale). Please install it."))

        cols <- circlize::colorRamp2(breaks=c(min(mat.c.cut), 0, max(mat.c.cut)),
                           colors=c("royalblue3","white","orange"))

        if (!is.null(hcl_palette)){
                cols <- circlize::colorRamp2(breaks=c(min(mat.c.cut), 0, max(mat.c.cut)),
                           , hcl_palette = hcl_palette,reverse = reverse)
        }    
    }

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

    format <- match.arg(format)

        if (format=="svg")
            grDevices::svg(file=paste0(path,filename,".svg")
             ,width=width, height=height) 

        if (format=="png")
            grDevices::png(file=paste0(path,filename,".png")
             ,width=width, height=height) 

        if (format=="pdf")
           grDevices::pdf(file=paste0(path,filename,".pdf")
                    , width=width, height=height,
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
    #if (!is.null(filename))
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
#' \code{"\link[=BSRInference-class]{scoreLRGeneSignatures}"}
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
#' print('simpleHeatmap')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' bsrinf.redP <- reduceToPathway(bsrinf) 
#' bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 
#' bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP,
#' qval.thres=0.001)
#'
#' scoresPathway <- scoreLRGeneSignatures(bsrdm,
#'                bsrsig.redPBP,name.by.pathway=TRUE)
#' 
#  # correlate with the immune microenvironment
#' data(immune.signatures, package="BulkSignalR")
#' imm.scores <- scoreSignatures(bsrdm, immune.signatures)
#' dualHeatmap(scoresPathway, imm.scores, width=6, height=9,
#'            file="SDC-LR-dualheatmap.pdf",
#' pointsize=4, vert.p=0.85)
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
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
        cols <- circlize::colorRamp2(breaks=c(min(mat.c), 0, max(mat.c)),
                           colors=c("royalblue3","white","orange"))
    }
    if (is.null(cols.e)){
        if (!requireNamespace("circlize",quietly=TRUE))
            stop(paste0("Package \"circlize\" needed for this function to ",
                        "work (generation of color scale). Please install it."))
        cols.e <- circlize::colorRamp2(breaks=c(min(mat.e), 0, max(mat.e)),
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
#' Scores generic gene signatures over the 
#' samples of a BSRDataModel object.
#'
#' @param ds         A BSRDataModel object.
#' @param ref.signatures           Gene signatures.
#' @param robust       A logical indicating that z-scores 
#' should be computed with median and MAD 
#' instead of mean and standard deviation.
#' @details This function relies on 
#' a simple average of gene z-scores over each signature.
#' It is no replacement for mode advanced methods
#' such as CYBERSORT or BisqueRNA. 
#' It is provided for convenience.
#' @return A matrix containing the scores of 
#' each gene signature in each sample.
#' Note that ligand-receptor gene
#' signature scores should be computed with 
#' \code{"\link[=BSRDataModel-class]{scoreLRGeneSignatures}"}
#' instead.
#' @export
#' @examples
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' data(immune.signatures, package="BulkSignalR")
#' imm.scores <- scoreSignatures(bsrdm, immune.signatures)
#' @importFrom methods is
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

#' Alluvial plot 
#'
#' @description Representation of the links
#' between Ligands,Receptors and Pathways.
#'
#' @param bsrinf object bsrinf inference.
#' @param keywords vector of pathways.
#' @param type filter on Ligand, Recptor or pathway id.
#' @param qval threshold over Q-value.
#' @param format svg / png / pdf. By default means it will 
#' plot in a svg format.
#' @param path directory where to plot file.
#' @param filename file name.
#' @param width width of image.
#' @param height height of image. 
#' @return NULL (plot printed somewhere on disk)
#' This is a convenience function that relies on the \code{ggalluvial}
#' package to propose a simple way
#' of representing Ligands, Receptors 
#  and underlying Pathways associated. 
#' @import ggplot2
#' @import ggalluvial
#' @export
#' @examples
#' print('alluvialPlot')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' alluvialPlot(bsrinf,
#'              keywords = c("COL4A1"),
#'              type = "L",
#'              qval = 0.001,
#'              path = "./",
#'              filename = "sdc_alluvial", 
#'              width  = 16, 
#'              height = 12
#'              )
alluvialPlot <- function(bsrinf,
                 keywords=c("COL4A1"),
                 type=c("L","R","pw.id"),
                 qval=1,
                 format=c("svg","png","pdf") ,
                 path="./",
                 filename = "alluvial.plot",
                 width = 10 ,
                 height = 10
                 ) {

     interactions <- data.frame(
           L =  unlist(ligands(bsrinf)) 
          ,R =  unlist(receptors(bsrinf))
          ,pw.name = LRinter(bsrinf)$pw.name
          ,pw.id = LRinter(bsrinf)$pw.id
          ,qval = LRinter(bsrinf)$qval
          ,targets = sapply(tGenes(bsrinf), paste, collapse = "  ")
        )


    if (!dir.exists(path)) {
        stop("Directory is invalid.")
    }

    type <- match.arg(type)
    format <- match.arg(format)

     if(type=="L")
          subset.interactions <- interactions[interactions$L %in% keywords,]
     if(type=="R")
          subset.interactions <- interactions[interactions$R %in% keywords,]
     if(type=="pw.id")
          subset.interactions <- interactions[interactions$pw.id %in% keywords,]  

       if (dim(subset.interactions)[1]==0){
        cat(paste(keywords, collapse = ' ')," for ", type, " not found.","\n", sep="")
        stop("Try another value for filtering.")
    }
     subset.interactions                <- subset.interactions[subset.interactions$qval < qval,]  
     subset.interactions$count <- 1
     subset.interactions       <- subset.interactions[,c("L" ,   "R"  ,  "pw.name" ,"count")]

     stratum <- ggalluvial::StatStratum

     pl <- ggplot2::ggplot(subset.interactions,
        ggplot2::aes_string(y = "count", axis1 = "L", axis2 = "R",axis3 = "pw.name")) +
        ggalluvial::geom_alluvium(ggplot2::aes_string(fill = "R"), width = 1/12) +
        ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
        ggplot2::geom_label(stat = stratum,  aes(label = ggplot2::after_stat(stratum))) +
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

    if (format=="svg")
        grDevices::svg(file=paste0(path,filename,".svg")
         ,width=width/2.54, height=height/2.54) 

    if (format=="png")
        grDevices::png(file=paste0(path,filename,".png")
         ,width=width/2.54, height=height/2.54) 

    if (format=="pdf")
       grDevices::pdf(file=paste0(path,filename,".pdf")
            ,width=width/2.54, height=height/2.54) 

    print(pl)

    grDevices::dev.off()
} #alluvialPlot


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
#' @param format svg / png / pdf. By default means it will 
#' plot in a svg format.
#' @param width width of image 
#' @param height height of image 
#' @return Circos Plot on Disk
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2 circos.par chordDiagramFromDataFrame
#' @importFrom circlize circos.trackPlotRegion circos.text circos.axis 
#' @importFrom circlize get.cell.meta.data circos.clear
#'
#' @export
#' @examples
#' print('chordDiagramLR')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' chordDiagramLR (bsrinf,
#'                  path="./",
#'                  filename="sdc_chord",
#'                  pw.id.filter="R-HSA-202733",
#'                  ligand="COL18A1",
#'                  receptor="ITGA3",
#'                  limit=20,
#'                  width=5, 
#'                  height=4.5
#'    )
chordDiagramLR  <- function(bsrinf,path="./",
    filename="chord",
    pw.id.filter="R-RSA-17821",
    ligand="L1",receptor="R1",
    limit=20,
    format=c("svg","png","pdf"),
    width=4,height=4) {

    format <- match.arg(format)

    print("chordDiagramLR")

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
    dataframe.bsrinf$pair <- paste(dataframe.bsrinf$ligands,
        dataframe.bsrinf$receptors,sep="-")

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

    cr <- circlize::colorRamp2(c(min(dataframe.bsrinf$corr),
                     max(dataframe.bsrinf$corr)), c("white","#febd17"))

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
        grDevices::svg(paste0(path,"/",filename,".svg"),width=width, height=height) #inch
    if (format=="png")
        grDevices::png(paste0(path,"/",filename,".png"),width=width, height=height)
    if (format=="pdf")
        grDevices::pdf(paste0(path,"/",filename,".pdf"),width=width, height=height)

    interactions <- data.frame(from=dataframe.bsrinf$ligands,
        to=dataframe.bsrinf$receptors,
        value=dataframe.bsrinf$corr)

    circlize::circos.par(points.overflow.warning=FALSE)

    circlize::chordDiagramFromDataFrame(interactions,
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
          link.arr.lty = "solid",
          link.arr.lwd = link.lwd, 
          link.arr.col = link.col,
          big.gap = 2, 
          small.gap = 1)
  
  circlize::circos.trackPlotRegion(track.index = 2, 
  panel.fun = function(x, y) {
  xlim = circlize::get.cell.meta.data("xlim")
  ylim = circlize::get.cell.meta.data("ylim")
  sector.name = circlize::get.cell.meta.data("sector.index")

  circlize::circos.text(mean(xlim), ylim[1] + 1.9, sector.name, 
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.7)

  circlize::circos.axis(h="top",labels=FALSE,minor.ticks=FALSE,
                    major.tick.length=1,
                    major.at=c(xlim), 
                    sector.index=sector.name,
                    track.index=2)
  })

  # LEGEND 

    lgd_points = Legend(labels = c("Ligands", "Receptors"),
        type = "points", pch = 16,
        legend_gp = grid::gpar(col = c("gray25","#7fbb00")),
        title_position = "topleft", 
        labels_gp = grid::gpar( font = 6),
        title = "LR") 
   
    lgd_links = Legend(at = c(round(min(interactions$value),
        digits = 2), round(max(interactions$value),
        digits = 2)), col_fun = cr, 
        title = "Correlation",direction ="horizontal"   ,   
        grid_width = unit(0.9, "mm") ,
        grid_height = unit(1.3, "mm") , 
        labels_gp = grid::gpar( font = 6),
        title_position = "topcenter",
       )

    lgd_list_vertical = packLegend(lgd_points)

    draw(lgd_list_vertical, x = unit(2, "mm"),
            y = unit(2, "mm"),
            just = c("left", "bottom"))
    draw(lgd_links, x = unit(2.7, "inch"),
            y = unit(2, "mm"),
            just = c("left", "bottom"))

     circlize::circos.clear()

     grDevices::dev.off()

}



#' Rasterize raw image from filepath 
#'
#' @param path.to.file path.to.file
#'
#' @return A raster image object
#'   
#' @importFrom png readPNG
#'
#' @export
#' @examples
#' print('rasterizeFromFile')
#' if (FALSE){
#'
#' my.image.as.raster <- rasterizeFromFile(path.to.file)  
#'
#' }
rasterizeFromFile <- function(path.to.file){

    img <- png::readPNG(path.to.file)

    my.image.as.raster <- grDevices::as.raster(img)  

    return(my.image.as.raster)
} # rasterizeFromFile


#' (counter-)clockwise 90 degree rotation
#'
#' @param rasterImage raster image
#' @param degrees degrees of rotation 
#'
#' @return A rotated raster image
#'   
#' @export
#' @examples
#' print('coordsFlip')
#' if(FALSE){
#' img <- png::readPNG(path.to.file)
#' my.image.as.raster <- grDevices::as.raster(path.to.file)  
#' my.image.rotated <- coordsFlip(my.image.as.raster,degrees=90)
#' 
#'  }
coordsFlip <- function(rasterImage, degrees=90) {

    if(! is.numeric(degrees))
         stop("Degrees must be numeric.", call. = FALSE)

    if(! degrees %% 90 == 0 )
         stop("Degrees must be divisible by 90", call. = FALSE)

    s <- sign(degrees)
    flip <- ifelse(s == 1, 
        \(x) t(apply(x, 2, rev)), # clockwise
        \(x) apply(x, 1, rev))    # counter-clockwise
    n <- abs(degrees / 90)
    for (i in seq_len(n)) 
        rasterImage <- flip(rasterImage)
     grDevices::as.raster(rasterImage)
}
# coordsFlip


#' Flip over horizontal/vertical axis
#'
#' @param rasterImage raster image
#' @param scale Reverse horizontal(x)/vertical(y) axis 
#'
#' @return A flipped raster image
#'   
#' @export
#' @examples
#' print('scalesReverse')
#' if(FALSE){
#' img <- png::readPNG(path.to.file)
#' my.image.as.raster <- grDevices::as.raster(path.to.file)  
#' my.image.flipped <- scalesReverse(my.image.as.raster,scale="x")
#' 
#'  }
scalesReverse <- function(rasterImage, scale = "x") {
    rasterImage <- if (scale == "x") {
        apply(rasterImage, 2, rev)
    } else {
        t(apply(rasterImage, 1, rev))
    }
     grDevices::as.raster(rasterImage)
}
# scalesReverse


#' L-R interaction score spatial display
#'
#' Generate a plot with scores at the spatial coordinates of the corresponding
#' sample locations. Not limited to BulkSignalR gene signature scores.
#'
#' @param v  A named vector containing the scores, names must be
#' the IDs of each location.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations as well as the unique IDs of spatial locations.
#' In case \code{ref.plot} is set to TRUE,
#' a label column is required additionally.
#' @param inter.name  Interaction name to display as plot title.
#' @param rev.y  A Boolean indicating whether low y coordinates should be
#' at the top of the plot.
#' @param ref.plot  A Boolean indicating whether a reference map of the tissue
#' with area labels should be plot aside.
#' @param ref.plot.only  A Boolean indicating that only the reference plot
#' should be output.
#' @param image.raster  Raster object image to plot raw tissue image as ref.
#' @param x.col  Column name in \code{areas} containing x coordinates.
#' @param y.col  Column name in \code{areas} containing y coordinates.
#' @param label.col  Column name in \code{areas} containing area labels.
#' @param idSpatial.col  Column name in \code{areas} containing the unique
#' IDs of spatial locations.
#' @param cut.p  Proportion of top and bottom values for thresholding.
#' @param low.color  Color for low score values.
#' @param mid.color  Color for score = 0.
#' @param high.color  Color for high score values.
#' @param title.fs Title font size.
#' @param legend.fs Legend items font size.
#' @param axis.fs Axis ticks font size.
#' @param label.fs Legend titles and axis names font size.
#' @param dot.size Dot size.
#' @return A single (scores) or side-by-side (reference tissue & scores) plot.
#' @export
#' @examples
#' print('spatialPlot')
#' if(FALSE){
#' img <- png::readPNG(path.to.file)
#'
#' spatialPlot(v,
#'     areas,
#'     inter.name="L->R",
#'     image.raster=img)
#' }
#' @import ggplot2
#' @import grid
#' @importFrom gridExtra grid.arrange
#' @importFrom scales rescale
spatialPlot <- function(v, areas, inter.name, rev.y=TRUE, ref.plot=FALSE,
                        ref.plot.only=FALSE, image.raster=NULL,
                        x.col="array_col", y.col="array_row",
                        label.col="label", idSpatial.col="idSpatial",
                        cut.p=0.01, low.color="royalblue3",
                        mid.color="white", high.color="orange",
                        title.fs=12, legend.fs=10, axis.fs=10,
                        label.fs=12, dot.size=0.5){

   x <- y <- label <- score <- NULL

  if (cut.p<0 || cut.p>0.1)
    stop("cut.p must lie in [0;0.1]")

  if (!all(c(x.col, y.col, label.col, idSpatial.col) %in% names(areas)))
    stop(paste0("One of x.col, y.col, idSpatial.col, or label.col is not in ",
                "names(areas)"))

  # put the scores in the right order for display at (x,y) coordinates
  v <- as.numeric(v[areas[[idSpatial.col]]])
  w <- .cutExtremeValues(v, cut.p)
  if (ref.plot || ref.plot.only)
    # include the labels
    tissue <- data.frame(x=areas[[x.col]], y=areas[[y.col]],
                         label=factor(areas[[label.col]]),
                         score=w)
  else
    # no label needed
    tissue <- data.frame(x=areas[[x.col]], y=areas[[y.col]],
                         score=w)

  if (rev.y){
    # revert y-axis
    ymin <- min(tissue$y)
    ymax <- max(tissue$y)
    tissue$y <- ymax-tissue$y+ymin
  }

  # reference tissue areas are required on the side
  if (ref.plot || ref.plot.only){
      if (is.null(image.raster))
        ref <- ggplot2::ggplot(data=tissue, ggplot2::aes(x=x, y=y)) +
        ggplot2::ggtitle("Reference tissue") +
        ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size) +
        ggplot2::theme_set(ggplot2::theme_bw(base_size = 10)) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=title.fs))
      
      else {ref <- grid::rasterGrob(image.raster)}
     
  }

  if (ref.plot.only){
    # do not plot L-R scores
    lr <- ref
    if (!is.null(image.raster))
        lr <- image.raster
    ref.plot <- FALSE
  }
  else{
    # plot L-R scores
    if (min(w) >= 0 || max(w) <= 0)
      # single gradient
      lr <- ggplot2::ggplot(data=tissue, ggplot2::aes(x=x, y=y)) +
        ggplot2::ggtitle(inter.name) +
        ggplot2::geom_point(ggplot2::aes(color=score), size=dot.size) +
        ggplot2::scale_color_gradient(low=low.color, high=high.color) +
        ggplot2::theme_set(ggplot2::theme_bw(base_size=10)) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=title.fs))
    else
      # dual gradient
      lr <- ggplot2::ggplot(data=tissue, ggplot2::aes(x=x, y=y)) +
        ggplot2::ggtitle(inter.name) +
        ggplot2::geom_point(ggplot2::aes(color=score), size=dot.size) +
        ggplot2::scale_color_gradientn(colors=c(low.color, mid.color, high.color),
values=scales::rescale(c(min(w)-1e-6, 0, max(w)+1e-6), c(0, 1))) +
        ggplot2::theme_set(ggplot2::theme_bw(base_size=10)) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=title.fs))
  }

  if (ref.plot)
    gridExtra::grid.arrange(ref, lr, ncol=2)
  else {
    if(!is.null(image.raster))
        plot(lr)
    if (is.null(image.raster))
        lr
    }

} # spatialPlot


#' Generate L-R interaction score spatial plots in a folder
#'
#' Generate a series of individual spatial score plots in a folder.
#' Not limited to BulkSignalR gene signature scores.
#'
#' @param scores  A matrix of scores, one L-R interaction per row and
#' spatial locations in the columns. This matrix is typically obtained
#' from BulkSignalR functions \code{scoreLRGeneSignatures} or \code{scScoring}.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations as well as the unique IDs of spatial locations.
#' In case \code{ref.plot} is set to TRUE,
#' a label column is required additionally.
#' @param plot.folder  The folder name in which the plot files will be written.
#' @param width  The width of each individual plot.
#' @param height  The height of each individual plot.
#' @param pointsize  PDF font point size.
#' @param rev.y  A Boolean indicating whether low y coordinates should be
#' at the top of the plot.
#' @param ref.plot  A Boolean indicating whether a reference map of the tissue
#' with area labels should be plot aside.
#' @param x.col  Column name in \code{areas} containing x coordinates.
#' @param y.col  Column name in \code{areas} containing y coordinates.
#' @param label.col  Column name in \code{areas} containing area labels.
#' @param idSpatial.col  Column name in \code{areas} containing the unique
#' IDs of spatial locations.
#' @param cut.p  Proportion of top and bottom values for thresholding.
#' @param low.color  Color for low score values.
#' @param mid.color  Color for score = 0.
#' @param high.color  Color for high score values.
#' @param title.fs Title font size.
#' @param legend.fs Legend items font size.
#' @param axis.fs Axis ticks font size.
#' @param label.fs Legend titles and axis names font size.
#' @param dot.size Dot size.
#' @return A set of PDF files are created in the provided folder.
#' @export
#' @examples
#' print('generateSpatialPlots')
#' if(FALSE){
#'
#' generateSpatialPlots(scores, areas, plot.folder)
#'
#' }
generateSpatialPlots <- function(scores, areas, plot.folder, width=5, height=3,
                                 pointsize=8, rev.y=TRUE, ref.plot=TRUE,
                                 x.col="array_col", y.col="array_row",
                                 label.col="label", idSpatial.col="idSpatial",
                                 cut.p=0.01, low.color="royalblue3",
                                 mid.color="white", high.color="orange",
                                 title.fs=12, legend.fs=10, axis.fs=10,
                                 label.fs=12, dot.size=0.5){

  for (i in seq_len(nrow(scores))){
    inter <- gsub("\\}","",gsub("\\{","",rownames(scores)[i]))
    fn <- gsub(" +/ +","-",inter,perl=TRUE)

    grDevices::pdf(paste0(plot.folder, "/interaction-plot-", fn), width=width,
        height=height, useDingbats=FALSE, pointsize=pointsize)
    spatialPlot(scores[i, areas[[idSpatial.col]]], areas, inter,
                rev.y=rev.y, ref.plot=ref.plot,
                x.col=x.col, y.col=y.col,
                label.col=label.col, idSpatial.col=idSpatial.col,
                cut.p=cut.p, low.color=low.color, mid.color=mid.color,
                high.color=high.color, title.fs=title.fs,
                legend.fs=legend.fs, axis.fs=axis.fs, label.fs=label.fs,
                dot.size=dot.size)
    grDevices::dev.off()
  }

} # generateSpatialPlots


#' Generate a visual index of spatial score distributions
#'
#' Generate an index made of series of small individual spatial score plots
#' in a PDF. Not limited to BulkSignalR gene signature scores.
#'
#' @param scores  A matrix of scores, one L-R interaction per row and
#' spatial locations in the columns. This matrix is typically obtained
#' from BulkSignalR functions \code{scoreLRGeneSignatures} or \code{scScoring}.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations, the unique IDs of spatial locations, and
#' a tissue label column.
#' @param out.file File name for the output PDF.
#' @param dot.size Dot size.
#' @param base.h  Width of each plot.
#' @param base.v  Height of each plot.
#' @param ratio the vertical/horizontal ratio.
#' @return A PDF file is created that contains the index.
#' @export
#' @examples
#' print('spatialIndexPlot')
#' if(FALSE){
#'
#' spatialIndexPlot(scores, areas, out.file)
#'
#' }
#' @importFrom gridExtra grid.arrange
spatialIndexPlot <- function(scores, areas, out.file,
                             dot.size=0.25, ratio=1.25,
                             base.v=2.5, base.h=3){

  # one reference plot at the beginning
  plots <- list(spatialPlot(scores[1,], areas, "",
                            ref.plot.only=TRUE, dot.size=dot.size))

  # actual plots
  for (i in seq_len(nrow(scores))){
    inter <- gsub("\\}", "", gsub("\\{", "", rownames(scores)[i]))
    plots <- c(plots, list(
      spatialPlot(scores[i,], areas, inter,
                  ref.plot=FALSE, dot.size=dot.size)
    ))
  }

  # file output at matrix dimension computation
  m <- length(plots)
  n <-  round(sqrt(m/ratio*base.v/base.h))
  l <- m %/% n
  if (l*n < m)
    l <- l+1
  grDevices::pdf(out.file, width=n*base.h, height=l*base.v, useDingbats=FALSE, pointsize=6)
  gridExtra::grid.arrange(grobs=plots, ncol=n, nrow=l)
  grDevices::dev.off()

} # spatialIndexPlot


#' Statistical association of scores with area labels
#'
#' Compute the statistical association of L-R interaction score spatial
#' distributions with tissue area labels.
#' Not limited to BulkSignalR gene signature scores.
#'
#' @param scores  A matrix of scores, one L-R interaction per row and
#' spatial locations in the columns. This matrix is typically obtained
#' from BulkSignalR functions \code{scoreLRGeneSignatures} or \code{scScoring}.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations, the unique IDs of spatial locations, and
#' a label column.
#' @param test  The chosen statistical test, parametric (normal) or
#' nonparametric (see also details below).
#' @param label.col  Column name in \code{areas} containing area labels.
#' @param idSpatial.col  Column name in \code{areas} containing the unique
#' IDs of spatial locations.
#' @param fdr.proc  Multiple hypothesis correction procedure, see
#' \code{multtest}.
#' @return A data.frame with the names of the interactions, the value of the
#' chosen statistics, and the corresponding Q-value. In case the
#' nonparametric Kruskal-Wallis test is chosen, additional columns are provided
#' testing each label for significantly larger scores (Kruskal-Wallis is global
#' and only says whether one or several labels show a bias). Individual
#' labels are tested with Wilcoxon and two columns are added *per* label,
#' one for the statistics and one for a Bonferroni-corrected P-value over
#' all the labels.
#'
#' In case the parametric model is chosen, a linear model followed by ANOVA
#' are used for the global association test. Individual labels are tested
#' with T-tests (Bonferroni-corrected).
#' @export
#' @examples
#' print('spatialAssociation')
#' if(FALSE){
#'
#' spatialAssociation(scores, areas)
#'
#' }
#' @import multtest
#' @importFrom foreach %do%
spatialAssociation <- function(scores, areas, test=c("Kruskal-Wallis","ANOVA"),
                               label.col="label", idSpatial.col="idSpatial",
                               fdr.proc=c("BH", "Bonferroni",
                               "Holm", "Hochberg", "SidakSS", "SidakSD", "BY",
                               "ABH", "TSBH")){

  test <- match.arg(test)
  fdr.proc <- match.arg(fdr.proc)
  if (!(label.col %in% names(areas)))
    stop("label.col is not in names(areas)")

  labels <- factor(areas[[label.col]])
  ul <- unique(labels)
  n <- length(ul)
  i <- lab <- NULL
  res <- foreach::foreach(i=seq_len(nrow(scores)), .combine=rbind) %do% {
    inter <- gsub("\\}","",gsub("\\{","",rownames(scores)[i]))
    v <- scores[i, areas[[idSpatial.col]]]
    if (test == "Kruskal-Wallis"){
      kw <- stats::kruskal.test(v,labels)

      # specific associations with each label
      pvals <- foreach::foreach(lab=ul, .combine=c) %do% {
        wt <- stats::wilcox.test(x=v[labels==lab], y=v[labels!=lab],
                          alternative="greater")
        list(min(c(wt$p.value*n, 1))) # Bonferroni correction
      }
      names(pvals) <- ul

      cbind(data.frame(interaction=inter, pval=kw$p.value, H=kw$statistic,
                       stringsAsFactors=FALSE),
            pvals)
    }
    else{
      # ANOVA
      df <- data.frame(score=v, label=labels)
      my.lm <- stats::lm(score ~ label, data=df)
      ano <- stats::anova(my.lm)

      # specific associations with each label
      pvals <- foreach::foreach(lab=ul, .combine=c) %do% {
        tt <- stats::t.test(x=v[labels==lab], y=v[labels!=lab],
                                 alternative="greater")
        list(min(c(tt$p.value*n, 1))) # Bonferroni correction
      }
      names(pvals) <- ul

      cbind(data.frame(interaction=inter, pval=ano$`Pr(>F)`[1],
                       F=ano$`F value`[1], stringsAsFactors=FALSE),
            pvals)
    }
  }

  # multiple hypothesis correction on the global association P-values
  rawp <- res$pval
  adj <- multtest::mt.rawp2adjp(rawp, fdr.proc)
  res$qval <- adj$adjp[order(adj$index), fdr.proc]

  rownames(res) <- res$interaction
 
  label.index.stop <- ncol(res)-1
  res[, c(1:3, ncol(res), 4:label.index.stop)] # put label columns at the end

} # spatialAssociation


#' Heatmap plot of association of scores with area labels
#'
#' Plot a heatmap featuring Q-values of statistical association between
#' L-R interaction score spatial distributions and tissue area labels.
#'
#' @param associations  A statistical association data.frame generated
#' by the function \code{spatialAssociation}.
#' @param qval.thres  The maximum Q-value to consider in the plot (a
#' L-R interaction must associate with one label at least with a Q-value
#' smaller or equal to this threshold).
#' @param  colors  A function returning a color for a given value such as
#' generated by \code{circlize::colorRamp2}.
#' @return Display a heatmap linking L-R interactions to labels.
#' @export
#' @examples
#' print('spatialAssociationPlot')
#' if(FALSE){
#'
#' spatialAssociationPlot(associations)
#'
#' }
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
spatialAssociationPlot <- function(associations, qval.thres=0.01, colors=NULL){

  # transform and filter data
  mat <- data.matrix(associations[, -(1:4)])
  mat[mat == 0] <- min(mat[mat > 0])
  mat <- -log10(mat)
  thres <- -log10(qval.thres)
  good <- apply(mat, 1, max) >= thres
  mat <- mat[good, ]

  if (is.null(colors))
    # create a color scale
    colscale <- circlize::colorRamp2(breaks=c(0, thres-1e-10,
                                    seq(thres, max(mat), length.out=10)),
                           colors=c("lightgray", "lightgray",
                                    grDevices::hcl.colors(10, "Viridis")))
  else
    colscale <- colors

  # plot heatmap
  di.lab <- stats::dist(t(mat))
  hc.lab <- stats::hclust(di.lab, method="ward.D")
  dend.lab <- stats::as.dendrogram(hc.lab)
  di.int <- stats::dist(mat)
  hc.int <- stats::hclust(di.int, method="ward.D")
  dend.int <- stats::as.dendrogram(hc.int)
  ComplexHeatmap::Heatmap(mat, col=colscale, cluster_rows=dend.int,
          cluster_columns=dend.lab, show_row_dend=FALSE,
          show_column_dend=FALSE)

} # spatialAssociationPlot


#' 2D-projection of spatial score distributions
#'
#' Use PCA or t-SNE to obtain a 2D-projection of a set of spatial scores.
#' This plot summarizes the diversity of patterns occuring in a spatial
#' dataset. Use the function \code{spatialIndexPlot} to create a large
#' visual index of many spatial distributions.
#' Not limited to BulkSignalR gene signature scores.
#'
#' @param scores  A matrix of scores, one L-R interaction per row and
#' spatial locations in the columns. This matrix is typically obtained
#' from BulkSignalR functions \code{scoreLRGeneSignatures} or \code{scScoring}.
#' @param associations  A statistical association data.frame generated
#' by the function \code{spatialAssociation}.
#' @param proj  Projection method : 'PCA' or 'tSNE' are available arguements.
#' @param qval.thres  The maximum Q-value to consider in the plot (a
#' L-R interaction must associate with one label at least with a Q-value
#' smaller or equal to this threshold).
#' @param with.names  A Boolean indicating whether L-R names should be plotted.
#' @param text.fs Point label font size in case \code{with.names} is TRUE.
#' @param legend.fs Legend items font size.
#' @param axis.fs Axis ticks font size.
#' @param label.fs Legend titles and axis names font size.
#' @param dot.size Dot size.
#' @return Display a 2D-projection of the score spatial distributions.
#' @export
#' @examples
#' print('spatialDiversityPlot')
#' if(FALSE){
#'
#' spatialDiversityPlot(scores,associations)
#'
#' }
#' @importFrom foreach %do%
#' @importFrom Rtsne Rtsne
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
spatialDiversityPlot <- function(scores, associations, proj=c("PCA","tSNE"),
                                 qval.thres=0.01, with.names=FALSE,
                                 text.fs=2.5, legend.fs=10, axis.fs=10,
                                 label.fs=12, dot.size=1){

  i <- PC1 <- PC2 <- name <- label <- tSNE1 <- tSNE2 <- NULL

  proj <- match.arg(proj)

  # find the strongest association for each L-R interaction
  labels <- names(associations)[-(1:4)]
  cols <- stats::setNames(c(grDevices::rainbow(length(labels), s=0.5), "lightgray"),
                          c(labels, "non_signif"))

  best.label <- foreach::foreach(i=seq_len(nrow(associations)),
                                 .combine=c) %do% {
    qvals <- as.numeric(associations[i, -(1:4)])
    if (min(qvals) > qval.thres)
      "non_signif"
    else
      labels[which.min(qvals)]
  }

  # the plot itself
  if (proj == "PCA"){
    pca <- stats::prcomp(scores, scale.=TRUE)
    dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], label=best.label, name=rownames(scores))
    if (with.names)
      p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=PC1, y=PC2, label=name)) +
      ggrepel::geom_text_repel(size=text.fs,max.overlaps = Inf) +
      ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size)
    else
      p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=PC1, y=PC2)) +
      ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size)

    p <- p + ggplot2::theme_set(ggplot2::theme_bw(base_size = 10)) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs))
    p
  }
  else{
    tsne <- Rtsne::Rtsne(scores, perplexity=10)
    dat <- data.frame(tSNE1=tsne$Y[,1], tSNE2=tsne$Y[,2], label=best.label, name=rownames(scores))
    if (with.names)
      p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=tSNE1, y=tSNE2, label=name)) +
      ggrepel::geom_text_repel(size=text.fs) +
      ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size)
    else
      p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=tSNE1, y=tSNE2)) +
      ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size)

    p <- p + ggplot2::theme_set(ggplot2::theme_bw(base_size = 10)) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs))
    p
  }

} # spatialDiversityPlot