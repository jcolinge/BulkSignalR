#' Smooth spatial expression data
#'
#' @param bsrdm  A BSRDataModel object containing the expression data to smooth.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations.
#' @param nnn    Number of nearest-neighbor locations to use for smoothing
#' each location. In case \code{radius} is set, then it is the maximum number
#' of nearest neighbors within the radius.
#' @param radius  A maximal distance to include neighbors in the smoothing.
#' @param weight.ratio  The weight given to the central location.
#' @param x.col  Column name in \code{areas} containing x coordinates.
#' @param y.col  Column name in \code{areas} containing y coordinates.
#' 
#' @details The expression data contained in a BSRDataModel object are
#' smoothed using a weighted average of nearby locations.
#' 
#' Two strategies are available to identify the neighbors. It is possible to
#' simply set the number of nearest-neighbors (parameter \code{nnn}). An
#' alternative consists in providing a distance radius (\code{radius}) along
#' with a a maximum number of nearest-neighbors within the radius
#' (\code{nnn.radius}). To properly define the radius, the user must know the
#' location coordinates. The strategy with the radius enables having corner
#' locations with two neighbors only and border locations with three
#' neighbors only, whereas to simply set a maximum of four neighbors for
#' instance would retrieve the four closest neighbors in every case.
#' 
#' For each location, its nearest-neighbors are found and a weighted average
#' computed with \code{weight.ratio} given to the central location itself
#' and a total weight of 1-\code{weight.ratio} shared within the neighbors
#' based on the inverse of their distances. In case \code{radius} is set,
#' some locations may have less than \code{nnn} neighbors (see above).
#' At such locations, the weight given to
#' the central location is augmented according to
#' 1-(1-\code{weight.ratio})*(number of neighbors)/\code{nnn}.
#' 
#' @export
#'
#' @return A BSRDataModel object containing the smoothed ncounts.
#'   
#' @importFrom RANN nn2
#'
#' @examples
#' print('smoothSpatialCounts')
#' if (FALSE){
#'
#' sm.bsrdm <- smoothSpatialCounts(bsrdm, areas, radius=1.2, nnn.radius=4)
#'
#' }
smoothSpatialCounts <- function(bsrdm, areas, nnn=4,
                                radius=NULL, weight.ratio=0.5,
                                x.col="array_col", y.col="array_row"){
  
  if (!is.null(param(bsrdm)$spatial.smooth) && (param(bsrdm)$spatial.smooth))
    warning("Spatial smoothing or ligand max has been already applied to these data")
  if (!all(c(x.col, y.col) %in% names(areas)))
    stop("One of x.col or y.col is not in names(areas)")
  if ((weight.ratio >= 1) || (weight.ratio <= 0))
    stop("weight.ratio must lie in ]0 ; 1[")
  
  # gets neighbor spots and their distances
  spots <- areas[, c(x.col, y.col)]
  if (is.null(radius))
    # RANN standard search
    neighb <- RANN::nn2(spots, k=nnn+1)
  else
    # RANN radius search
    neighb <- RANN::nn2(spots, k=nnn+1, radius=radius, searchtype="radius")
  
  # determine weights
  c <- neighb$nn.idx
  c[, 1] <- 0
  c[c > 0] <- 1
  ratios <- 1 - (1-weight.ratio)*rowSums(c)/nnn
  d <- 1 / (neighb$nn.dists * c)[, -1]
  d[is.infinite(d)] <- 0
  tot <- rowSums(d) / (1 - ratios)
  weights <- sweep(d, 1, tot, "/")
  weights[is.nan(weights)] <- 0
  weights <- cbind(ratios, weights)
  
  # smooth the transcriptomes
  orig <- ncounts(bsrdm)
  smooth <- orig
  for (i in seq_len(ncol(orig))){
    used <- neighb$nn.idx[i,] > 0
    if (sum(used) > 1)
      smooth[, i] <- orig[, neighb$nn.idx[i, used]] %*% weights[i, used]
    else
      smooth[, i] <- orig[, neighb$nn.idx[i, used]] * weights[i, used]
  }

  ncounts(bsrdm) <- smooth
  bsrdm@param$spatial.smooth <- TRUE
  bsrdm
  
} # smoothSpatialCounts


#' Get maximal ligand expression at nearby locations
#'
#' @param bsrdm  A BSRDataModel object containing the expression data to smooth.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations.
#' @param nnn    Number of nearest-neighbor locations to use for smoothing
#' each location. In case \code{radius} is set, then it is the maximum number
#' of nearest neighbors within the radius.
#' @param radius  A maximal distance to include neighbors in the smoothing.
#' @param x.col  Column name in \code{areas} containing x coordinates.
#' @param y.col  Column name in \code{areas} containing y coordinates.
#' 
#' @details Ligand expression data contained in a BSRDataModel object are
#' modified to consider the possibility that the ligand of a L-R interaction
#' might be expressed at nearby locations. This is achieved replacing each
#' ligand expression by its maximum over the central location and its
#' neighbors. Since ligands and receptors are never used as gene targets
#' in computing the receptor downstream signal correlations, this
#' substitution is compatible with our statistical model. Moreover,
#' the reciprocal configuration where the ligand is expressed at the
#' central location and hits a receptors at a neighbor location is
#' covered when the same ligand maximization scheme is applied to
#' the neighbor. L-R localization and gene signature scoring is defined
#' by the location at which the receptor is expressed after applying
#' this function.
#' 
#' Two strategies are available to identify the neighbors. It is possible to
#' simply set the number of nearest-neighbors (parameter \code{nnn}). An
#' alternative consists in providing a distance radius (\code{radius}) along
#' with a a maximum number of nearest-neighbors within the radius
#' (\code{nnn.radius}). To properly define the radius, the user must know the
#' location coordinates. The strategy with the radius enables having corner
#' locations with two neighbors only and border locations with three
#' neighbors only, whereas to simply set a maximum of four neighbors for
#' instance would retrieve the four closest neighbors in every case.
#' 
#' @export
#'
#' @return A BSRDataModel object containing the maximized ligand expressions.
#'   
#' @importFrom RANN nn2
#'
#' @examples
#' print('maxLigandSpatialCounts')
#' if (FALSE){
#'
#' max.bsrdm <- maxLigandSpatialCounts(bsrdm, areas, radius=1.2, nnn.radius=4)
#'
#' }
maxLigandSpatialCounts <- function(bsrdm, areas, nnn=4, radius=NULL,
                                   x.col="array_col", y.col="array_row"){
  
  if (!is.null(param(bsrdm)$spatial.smooth) && (param(bsrdm)$spatial.smooth))
    warning("Spatial smoothing or ligand max has been already applied to these data")
  if (!all(c(x.col, y.col) %in% names(areas)))
    stop("One of x.col or y.col is not in names(areas)")
  
  # gets neighbor spots and their distances
  spots <- areas[, c(x.col, y.col)]
  if (is.null(radius))
    # RANN standard search
    neighb <- RANN::nn2(spots, k=nnn+1)
  else
    # RANN radius search
    neighb <- RANN::nn2(spots, k=nnn+1, radius=radius, searchtype="radius")
  
  # maximize ligands in the transcriptomes
  orig <- ncounts(bsrdm)
  maxi <- orig
  ligands <- intersect(LRdb$ligand, rownames(orig))
  for (i in seq_len(ncol(orig))){
    used <- neighb$nn.idx[i,] > 0
    if (sum(used) > 1)
      maxi[ligands, i] <- apply(orig[ligands, neighb$nn.idx[i, used]], 1, max)
    else
      maxi[ligands, i] <- orig[ligands, i]
  }
  
  bsrdm@ncounts <- maxi
  bsrdm@param$spatial.smooth <- TRUE
  bsrdm
  
} # maxLigandSpatialCounts


#' Rasterize raw image from filepath 
#'
#' @param path.to.file Path to an image file.
#'
#' @return A raster image object.
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


#' (Counter-)clockwise 90 degree rotation
#'
#' @param rasterImage Raster image.
#' @param degrees Rotation angle in degrees. 
#'
#' @return A rotated raster image.
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
  
} # coordsFlip


#' Flip over horizontal/vertical axis
#'
#' @param rasterImage Raster image.
#' @param axis Axis along which the image is flipped.
#'
#' @return A flipped raster image
#'   
#' @export
#' @examples
#' print('scalesReverse')
#' if(FALSE){
#' img <- png::readPNG(path.to.file)
#' my.image.as.raster <- grDevices::as.raster(path.to.file)  
#' my.image.flipped <- scalesReverse(my.image.as.raster, axis="x")
#' 
#'  }
scalesReverse <- function(rasterImage, axis = "x") {
  
    rasterImage <- if (axis == "x") {
        apply(rasterImage, 2, rev)
    } else {
        t(apply(rasterImage, 1, rev))
    }
    grDevices::as.raster(rasterImage)
  
} # scalesReverse


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
#' @param image.raster  Raster object image to plot raw tissue image as
#' reference.
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
#' @param legend.dot.factor A factor applied to obtain the legend dot size.
#' @param ref.colors A vector of colors to bypass those automatically chosen
#' by ggplot2 for the tissue areas in the reference plot.
#' @details A single (scores) or side-by-side (reference tissue & scores) plot
#' is generated.
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
                        label.fs=12, dot.size=0.5, legend.dot.factor=10,
                        ref.colors=NULL){
  
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
    if (is.null(image.raster)){
      ref <- ggplot2::ggplot(data=tissue, ggplot2::aes(x=x, y=y)) +
        ggplot2::ggtitle("Reference tissue") +
        ggplot2::geom_point(ggplot2::aes(color=label), size=dot.size) +
        ggplot2::theme_set(ggplot2::theme_bw(base_size = 10)) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=axis.fs)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=label.fs)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legend.fs)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=title.fs)) +
        ggplot2::guides(colour=ggplot2::guide_legend(
                        override.aes=list(size=dot.size*legend.dot.factor)))
        if (!is.null(ref.colors)){
          if (length(unique(tissue$label)) != length(ref.colors))
            stop("The number of reference colors do not match the number of labels")
          ref <- ref + ggplot2::scale_color_manual(values=ref.colors)
        }
      }
    else
      ref <- grid::rasterGrob(image.raster)
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
#' @param image.raster  Raster object image to plot raw tissue image as
#' reference.
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
#' @param ref.colors A vector of colors to bypass those automatically chosen
#' by ggplot2 for the tissue areas in the reference plot.
#' @details A set of PDF files are created in the provided folder.
#' @export
#' @examples
#' print('generateSpatialPlots')
#' if(FALSE){
#'
#' generateSpatialPlots(scores, areas, plot.folder)
#'
#' }
generateSpatialPlots <- function(scores, areas, plot.folder, width=5, height=3,
                                 pointsize=8, rev.y=TRUE,ref.plot=TRUE,image.raster=NULL,
                                 x.col="array_col", y.col="array_row",
                                 label.col="label", idSpatial.col="idSpatial",
                                 cut.p=0.01, low.color="royalblue3",
                                 mid.color="white", high.color="orange",
                                 title.fs=12, legend.fs=10, axis.fs=10,
                                 label.fs=12, dot.size=0.5,ref.colors=NULL){
  
  for (i in seq_len(nrow(scores))){
    inter <- gsub("\\}","",gsub("\\{","",rownames(scores)[i]))
    fn <- gsub(" +/ +","-",inter,perl=TRUE)
    
    grDevices::pdf(paste0(plot.folder, "/interaction-plot-", fn), width=width,
                   height=height, useDingbats=FALSE, pointsize=pointsize)
    figure <- spatialPlot(scores[i, areas[[idSpatial.col]]], areas, inter,
                rev.y=rev.y, ref.plot=ref.plot, 
                image.raster=image.raster,
                x.col=x.col, y.col=y.col,
                label.col=label.col, idSpatial.col=idSpatial.col,
                cut.p=cut.p, low.color=low.color, mid.color=mid.color,
                high.color=high.color, title.fs=title.fs,
                legend.fs=legend.fs, axis.fs=axis.fs, label.fs=label.fs,
                dot.size=dot.size,ref.colors=ref.colors)
    if(ref.plot)
      figure
    else 
      print(figure)

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
#' @param image.raster  Raster object image to plot raw tissue image as
#' reference.
#' @param ref.plot  A Boolean indicating whether a reference map of the tissue
#' with area labels should be plot first.
#' @param cut.p  Proportion of top and bottom values for thresholding.
#' @param low.color  Color for low score values.
#' @param mid.color  Color for score = 0.
#' @param high.color  Color for high score values.
#' @param title.fs Title font size.
#' @param legend.fs Legend items font size.
#' @param axis.fs Axis ticks font size.
#' @param label.fs Legend titles and axis names font size.
#' @param dot.size Dot size.
#' @param ref.colors A vector of colors to bypass those automatically chosen
#' by ggplot2 for the tissue areas in the reference plot.
#' @param x.col  Column name in \code{areas} containing x coordinates.
#' @param y.col  Column name in \code{areas} containing y coordinates.
#' @param label.col  Column name in \code{areas} containing area labels.
#' @param idSpatial.col  Column name in \code{areas} containing the unique
#' IDs of spatial locations.
#' @param base.h  Width of each plot.
#' @param base.v  Height of each plot.
#' @param ratio the vertical/horizontal ratio.
#' @details A PDF file is created that contains the index.
#' @export
#' @examples
#' print('spatialIndexPlot')
#' if(FALSE){
#'
#' spatialIndexPlot(scores, areas, out.file)
#'
#' }
#' @import grid
#' @importFrom gridExtra grid.arrange
spatialIndexPlot <- function(scores, areas, out.file, ref.plot=TRUE,
                             image.raster = NULL,
                             x.col="array_col", y.col="array_row",
                             label.col="label", idSpatial.col="idSpatial",
                             cut.p=0.01,low.color="royalblue3",
                             mid.color="white", high.color="orange",
                             title.fs=12, legend.fs=10, axis.fs=10,
                             label.fs=12, dot.size=0.25,ratio=1.25,
                             base.v=2.5, base.h=3, ref.colors=NULL){
  
  # one reference plot at the beginning
  if (ref.plot)
    if(is.null(image.raster))
      plots <- list(spatialPlot(scores[1,], areas, "",
                    x.col=x.col, y.col=y.col,
                    label.col=label.col, idSpatial.col=idSpatial.col,
                    ref.plot.only=TRUE,cut.p=cut.p,
                    low.color=low.color,
                    mid.color=mid.color, high.color=high.color,
                    title.fs=title.fs, legend.fs=legend.fs, 
                    axis.fs=axis.fs, label.fs=label.fs,
                    dot.size=dot.size, ref.colors=ref.colors))
    else 
      plots <- list(grid::rasterGrob(image.raster))
  else
    plots <- list()

  # actual plots
  for (i in seq_len(nrow(scores))){
    inter <- gsub("\\}", "", gsub("\\{", "", rownames(scores)[i]))
    plots <- c(plots, list(
      spatialPlot(scores[i,], areas, inter,
                  x.col=x.col, y.col=y.col,
                  label.col=label.col, idSpatial.col=idSpatial.col,
                  ref.plot=FALSE,cut.p=cut.p,
                  low.color=low.color,
                  mid.color=mid.color, high.color=high.color,
                  title.fs=title.fs, legend.fs=legend.fs, 
                  axis.fs=axis.fs,
                  label.fs=label.fs,
                  dot.size=dot.size)
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


#' Generate separated plots for a L-R interaction
#'
#' Generate a detailed view related to a chosen interaction made of series of
#' small individual spatial plots: tissue organization (optional), gene
#' signature score, ligand and receptor expression.
#'
#' @param v  A named vector containing the gene signature scores for the
#' L-R interaction including the contribution of the pathway, names must be
#' the IDs of each location. Alternatively, v can be a gene signature score
#' matrix such as those returned by \code{scoreLRGeneSignatures} and the
#' row named "{\code{L}} / {\code{R}}" will be used.
#' @param L  The name of the ligand.
#' @param R  The name of the receptor.
#' @param ncounts  The (normalized) expression matrix with column names equal
#' to the IDs of each location.
#' @param areas  A data.frame containing at least the x and y
#' coordinates of the locations as well as the unique IDs of spatial locations.
#' In case \code{ref.plot} is set to TRUE,
#' a label column is required additionally.
#' @param inter.name  Interaction name to display as plot title,
#' equal to "\code{L} / \code{R}" unless specified.
#' @param rev.y  A Boolean indicating whether low y coordinates should be
#' at the top of the plot.
#' @param ref.plot  A Boolean indicating whether a reference map of the tissue
#' with area labels should be plot aside.
#' @param image.raster  Raster object image to plot raw tissue image as
#' reference.
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
#' @param ref.colors A vector of colors to bypass those automatically chosen
#' by ggplot2 for the tissue areas in the reference plot.
#' @param legend.dot.factor A factor applied to obtain the legend dot size.
#' @details A set of spatial plots are generated including an optional
#' reference tissue plot (image or areas represented), the gene signature
#' scores, the ligand expression values, and the receptor expression values.
#' @export
#' @examples
#' print('separatedLRPlot')
#' if(FALSE){
#' img <- png::readPNG(path.to.file)
#'
#' separatedLRPlot(v,
#'     L,
#'     R,
#'     ncounts,
#'     areas,
#'     image.raster=img)
#' }
#' @import grid
#' @importFrom gridExtra grid.arrange
separatedLRPlot <- function(v, L, R, ncounts, areas, inter.name=NULL, rev.y=TRUE,
                            ref.plot=TRUE, image.raster=NULL,
                            x.col="array_col", y.col="array_row",
                            label.col="label", idSpatial.col="idSpatial",
                            cut.p=0.01, low.color="royalblue3",
                            mid.color="white", high.color="orange",
                            title.fs=12, legend.fs=10, axis.fs=10,
                            label.fs=12, dot.size=0.5, legend.dot.factor=10,
                            ref.colors=NULL){

  if (is.matrix(v))
    v <- v[paste0("{",L,"} / {",R,"}"),]
  if (is.null(inter.name))
    inter.name <- paste(L,"/",R)
  
  # one reference plot at the beginning
  if (ref.plot)
    if(is.null(image.raster))
      plots <- list(spatialPlot(v, areas, "", ref.plot.only=TRUE,
                                x.col=x.col, y.col=y.col,
                                label.col=label.col, idSpatial.col=idSpatial.col,
                                dot.size=dot.size,
                                legend.dot.factor=legend.dot.factor,
                                ref.colors=ref.colors))
    else 
      plots <- list(grid::rasterGrob(image.raster))
  else
    plots <- list()
  
  # gene signature scores
  plots <- c(plots, list(
    spatialPlot(v, areas, inter.name, rev.y=rev.y, ref.plot=FALSE,
                x.col=x.col, y.col=y.col,
                label.col=label.col, idSpatial.col=idSpatial.col,
                cut.p=cut.p, low.color=low.color,
                mid.color=mid.color, high.color=high.color,
                title.fs=title.fs, legend.fs=legend.fs, axis.fs=axis.fs,
                label.fs=label.fs, dot.size=dot.size)
  ))
  
  # ligand and receptor plots
  plots <- c(plots, list(
    spatialPlot(ncounts[L,], areas, inter.name=L, rev.y=rev.y, ref.plot=FALSE,
                x.col=x.col, y.col=y.col,
                label.col=label.col, idSpatial.col=idSpatial.col,
                cut.p=cut.p, low.color=low.color,
                mid.color=mid.color, high.color=high.color,
                title.fs=title.fs, legend.fs=legend.fs, axis.fs=axis.fs,
                label.fs=label.fs, dot.size=dot.size),
    spatialPlot(ncounts[R,], areas, inter.name=R, rev.y=rev.y, ref.plot=FALSE,
                x.col=x.col, y.col=y.col,
                label.col=label.col, idSpatial.col=idSpatial.col,
                cut.p=cut.p, low.color=low.color,
                mid.color=mid.color, high.color=high.color,
                title.fs=title.fs, legend.fs=legend.fs, axis.fs=axis.fs,
                label.fs=label.fs, dot.size=dot.size)
  ))

  # display  
  gridExtra::grid.arrange(grobs=plots, ncol=2, nrow=2)

} # separatedLRPlot


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
#' @param test  The chosen statistical test or statistics
#' (see details below).
#' @param label.col  Column name in \code{areas} containing area labels.
#' @param idSpatial.col  Column name in \code{areas} containing the unique
#' IDs of spatial locations.
#' @param fdr.proc  Multiple hypothesis correction procedure, see
#' \code{multtest}.
#' @return A data.frame with the names of the interactions, the value of the
#' chosen statistics, and the corresponding Q-value.
#' @details In case the
#' nonparametric Kruskal-Wallis test is chosen, additional columns are provided
#' testing each label for significantly larger scores (Kruskal-Wallis is global
#' and only says whether one or several labels show a bias). Individual
#' labels are tested with Wilcoxon and two columns are added *per* label,
#' one for the statistics and one for a Bonferroni-corrected P-value over
#' all the labels.
#'
#' In case an actual statistical test is chosen, a parametric test (ANOVA) and
#' a non-parametric test (Kruskal-Wallis) are available for the global analysis.
#' Individual labels are tested with T-tests or Wilcoxon (Bonferroni-corrected)
#' accordingly.
#'
#' In case a statistics is preferred, Spearman correlation or explained variance
#' (r2 or coefficient of determination, through linear models) are are available.
#' They mesure the relationship
#' between each individual area and \code{scores}. For the explained variance,
#' a global value (R2) is also computed from a multi-linear model (the same as
#' what is used for the ANOVA).
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
spatialAssociation <- function(scores, areas, test=c("Kruskal-Wallis","ANOVA","Spearman","r2"),
                               label.col="label", idSpatial.col="idSpatial",
                               fdr.proc=c("BH", "Bonferroni",
                                          "Holm", "Hochberg", "SidakSS", "SidakSD", "BY",
                                          "ABH", "TSBH")){
  
  test <- match.arg(test)
  fdr.proc <- match.arg(fdr.proc)
  if (!(label.col %in% names(areas)))
    stop("label.col is not in names(areas)")
  if (!(idSpatial.col %in% names(areas)))
    stop("idSpatial.col is not in names(areas)")
  
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
    else if (test == "ANOVA"){
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
    else if (test=="Spearman") {
      # Spearman correlation
      
      # specific correlations with each label
      corrs <- foreach::foreach(lab=ul, .combine=c) %do% {
        local <- rep(0, length(labels))
        local[labels == lab] <- 1
        co <- stats::cor(v, local, method="spearman")
        list(co)
      }
      names(corrs) <- ul
      
      cbind(data.frame(interaction=inter, stringsAsFactors=FALSE), corrs)
    }
    else{
      # r2 from linear regressions
      
      df <- data.frame(score=v, label=labels)
      my.lm <- stats::lm(score ~ label, data=df)
      ano <- stats::anova(my.lm)
      R2 <- ano$`Sum Sq`[1]/sum(ano$`Sum Sq`)

      # specific r2 for each label
      r2s <- foreach::foreach(lab=ul, .combine=c) %do% {
        local <- rep(0, length(labels))
        local[labels == lab] <- 1
        dfl <- data.frame(score=v, local=local)
        my.llm <- stats::lm(score ~ local, data=dfl)
        lano <- stats::anova(my.llm)
        r2 <- lano$`Sum Sq`[1]/sum(lano$`Sum Sq`)
        list(r2)
      }
      names(r2s) <- ul
      
      cbind(data.frame(interaction=inter, global.R2=R2, stringsAsFactors=FALSE),
            r2s)
    }
  }
  
  # multiple hypothesis correction on the global association P-values
  if (test %in% c("Kruskal-Wallis","ANOVA")){
    rawp <- res$pval
    adj <- multtest::mt.rawp2adjp(rawp, fdr.proc)
    res$qval <- adj$adjp[order(adj$index), fdr.proc]
    label.index.stop <- ncol(res)-1
    res <- res[, c(1:3, ncol(res), 4:label.index.stop)] # put Q-values in column 4
  }

  rownames(res) <- res$interaction
  res  
  
} # spatialAssociation


#' Heatmap plot of association of scores with area labels
#'
#' Plot a heatmap featuring Q-values or values of statistical association between
#' L-R interaction score spatial distributions and tissue area labels.
#'
#' @param associations  A statistical association data.frame generated
#' by the function \code{spatialAssociation}.
#' @param qval.thres  The maximum Q-value to consider in the plot (a
#' L-R interaction must associate with one label at least with a Q-value
#' smaller or equal to this threshold).
#' @param absval.thres  The minimum value to consider in the plot (a
#' L-R interaction must associate with one label at least with an absolute
#' value larger or equal to this threshold).
#' @param  colors  A function returning a color for a given value such as
#' generated by \code{circlize::colorRamp2}.
#' @details Display a heatmap linking L-R interactions to labels.
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
spatialAssociationPlot <- function(associations, qval.thres=0.01, absval.thres=0,
                             colors=NULL){
  
  # transform and filter data
  if (sum(c("pval","qval") %in% names(associations)) == 2){
    # log-scale on Q-values
    mat <- data.matrix(associations[, -(1:4)])
    mat[mat == 0] <- min(mat[mat > 0])
    mat <- -log10(mat)
    thres <- -log10(abs(qval.thres))
    good <- apply(mat, 1, max) >= thres
    mat <- mat[good, ]
  }
  else{
    # linear scale
    if ("global.R2" %in% names(associations))
      mat <- data.matrix(associations[, -(1:2)])
    else
      mat <- data.matrix(associations[, -1])
    thres <- abs(absval.thres)
    good <- apply(abs(mat), 1, max) >= thres
    mat <- mat[good, ]
  }
  
  if (is.null(colors))
    # create a color scale
    if (min(mat) >= 0)
      colscale <- circlize::colorRamp2(breaks=c(0, thres-1e-10,
                                              seq(thres, max(mat), length.out=10)),
                                     colors=c("lightgray", "lightgray",
                                              grDevices::hcl.colors(10, "Viridis")))
    else
      if (thres == 0)
        colscale <- circlize::colorRamp2(breaks=c(min(mat), 0, max(mat)),
                                         colors=c("royalblue", "white", "orangered"))
      else
        colscale <- circlize::colorRamp2(breaks=c(min(mat), -thres, thres, max(mat)),
                                         colors=c("royalblue", "white", "white", "orangered"))
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
#' Use PCA or t-SNE to obtain a 2D-projection of a set of spatial scores
#' or associations.
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
#' @param score.based  A logical indicating whether the plot should be
#' based on scores or the associations directly.
#' @param qval.thres  The maximum Q-value to consider in the plot (a
#' L-R interaction must associate with one label at least with a Q-value
#' smaller or equal to this threshold). Relevant for Kruskal-Wallis and
#' ANOVA tests in \code{spatialAssociation}.
#' @param val.thres  The minimum value to consider in the plot (a
#' L-R interaction must associate with one label at least with a value
#' larger or equal to this threshold). Relevant for Spearman and r2
#' associations in \code{spatialAssociation}.
#' @param with.names  A logical indicating whether L-R names should be plotted.
#' @param text.fs Point label font size in case \code{with.names} is TRUE.
#' @param legend.fs Legend items font size.
#' @param axis.fs Axis ticks font size.
#' @param label.fs Legend titles and axis names font size.
#' @param dot.size Dot size.
#' @param perplexity  Perplexity parameter for t-SNE.
#' @details Display a 2D-projection of the score spatial distributions.
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
                           score.based=FALSE,
                           qval.thres=0.01, val.thres=0, with.names=FALSE,
                           text.fs=2.5, legend.fs=10, axis.fs=10,
                           label.fs=12, dot.size=1,
                           perplexity=10){
  
  i <- PC1 <- PC2 <- name <- label <- tSNE1 <- tSNE2 <- log.scale <- remove.col <- NULL
  
  proj <- match.arg(proj)
  
  # adapt to the type of associations
  if (sum(c("pval","qval") %in% names(associations)) == 2){
    remove.col <- 1:4
    log.scale <- TRUE
  }
  else{
    # linear scale
    if ("global.R2" %in% names(associations))
      remove.col <- 1:2
    else
      remove.col <- 1
    log.scale <- FALSE
  }
  labels <- names(associations)[-remove.col]
  
  cols <- stats::setNames(c(grDevices::rainbow(length(labels), s=0.5), "lightgray"),
                          c(labels, "non_signif"))
  
  # find the strongest association for each L-R interaction
  best.label <- foreach::foreach(i=seq_len(nrow(associations)),
                                 .combine=c) %do% {
                                   if (log.scale){
                                     qvals <- as.numeric(associations[i, -remove.col])
                                     if (min(qvals) > qval.thres)
                                       "non_signif"
                                     else
                                       labels[which.min(qvals)]
                                   }
                                   else{
                                     values <- as.numeric(associations[i, -remove.col])
                                     if (max(values) < val.thres)
                                       "non_signif"
                                     else
                                       labels[which.max(values)]
                                   }
                                 }
  
  # the plot itself
  if (proj == "PCA"){
    if (score.based)
      pca <- stats::prcomp(scores, scale.=TRUE)
    else
      if (log.scale)
        pca <- stats::prcomp(-log10(data.matrix(associations[, -remove.col])),
                             scale.=TRUE)
      else
        pca <- stats::prcomp(data.matrix(associations[, -remove.col]), scale.=TRUE)
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
    if (score.based)
      tsne <- Rtsne::Rtsne(scores, perplexity=perplexity)
    else
      if (log.scale)
        tsne <- Rtsne::Rtsne(-log10(data.matrix(associations[, -remove.col])),
                             perplexity=perplexity, check_duplicates=FALSE)
      else
        tsne <- Rtsne::Rtsne(data.matrix(associations[, -remove.col]),
                             perplexity=perplexity, check_duplicates=FALSE)
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
