# Pure display functions (for convenience) =====================================

#' Internal function to cut extreme values from a matrix
#'
#' @param m         A matrix.
#' @param p          Proportion of top and bottom values for thresholding.
#' @return A matrix with values beyond top and bottom thresholds repaced by
#' the latter thresholds.
#' @keywords internal
.cutExtremeValues <- function(m, p) {
    thres.lo <- stats::quantile(m, prob = p)
    thres.hi <- stats::quantile(m, prob = 1 - p)
    m[m > thres.hi] <- thres.hi
    m[m < thres.lo] <- thres.lo
    m
}

#' Bubble Plot to explore LR & Pathways
#'
#' Quick check to observe LR - Pathways association
#' with their respective correlation
#' and Q-values.
#'
#' @param bsrinf     BulkSignalR inference object.
#' @param pathways  Vector of pathway names to keep.
#' @param qval.thres     Maximum Q-value.
#' @param filter.L     Vector of ligands to keep.
#' @param filter.R     Vector of receptors to keep.
#' @param color     Main color used for the gradient.
#' @param pointsize Global point size.
#'
#' @return  A bubble plot displayed in the current viewport.
#'
#' This is a convenience function to propose a simple way
#' of representing LR - Pathways association
#' with their respective correlation
#' and Q-values.
#'
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' pathways <- LRinter(bsrinf)[1,c("pw.name")]
#' bubblePlotPathwaysLR(bsrinf,
#' pathways = pathways,
#' qval.thres = 0.1,
#' color = "red",
#' pointsize = 8
#' )
#' @import ggplot2
bubblePlotPathwaysLR <- function(
    bsrinf,
    pathways,
    qval.thres = 1,
    filter.L = NULL,
    filter.R = NULL,
    color = "#16a647",
    pointsize = 6) {
    filtered.brinf <- LRinter(bsrinf)
    filtered.brinf <- filtered.brinf[filtered.brinf$qval < qval.thres, ]

    if (!is.null(filter.R) | !is.null(filter.L)) {
        filtered.brinf <- filtered.brinf[filtered.brinf$L %in% filter.L | 
        filtered.brinf$R %in% filter.R, ]
    }

    filtered.brinf$LR <- paste(filtered.brinf$L, filtered.brinf$R, sep = " / ")

    filtered.brinf <- filtered.brinf[, c("LR", "pw.name", "LR.corr", "qval")]
    filtered.brinf$log10.qval <- -log10(filtered.brinf$qval)
    filtered.brinf$log10.LR.corr <- -log10(abs(filtered.brinf$LR.corr)) 

    filtered.brinf <- filtered.brinf[filtered.brinf$pw.name %in% pathways, ]

    if (dim(filtered.brinf)[1] == 0) {
        stop("The pathways you have selected do no exist.")
    }

    limit.P <- 8
    if (length(pathways) >= limit.P) {
        message("We recommend less than ", limit.P, " pathways.")
        message("Too many pathways were given in input.")
    }
    limit.LR <- 50
    if (length(unique(filtered.brinf$LR)) > limit.LR) {
        message("Too many LR interactions detected (",
            length(unique(filtered.brinf$LR)), ").")
        message("We recommend less than ", limit.LR,
            " LR interactions to visualize.")
        message("Try to reduce (Qval-Threshold, number of pathways...).")
    }

    message(length(unique(filtered.brinf$LR)), " LR interactions detected.")

    plot(ggplot2::ggplot(
        filtered.brinf,
        ggplot2::aes_string(x = "LR", y = "pw.name")
    ) +
        ggplot2::geom_point(ggplot2::aes_string(
            size = "log10.LR.corr",
            fill = "log10.qval"
        ), alpha = 0.75, shape = 21) +
        labs(
            x = "", y = "", size = "-log10 (LR.corr)",
            fill = "-log10 (Qval)"
        ) +
        ggplot2::theme(
            legend.key = ggplot2::element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.position = "right", legend.box = "horizontal",
            axis.text.x = ggplot2::element_text(
                colour = "black",
                size = pointsize, angle = 90, vjust = 0.3, hjust = 1
            ),
            axis.text.y = ggplot2::element_text(
                colour = "black",
                face = "bold", size = pointsize
            ),
            legend.text = ggplot2::element_text(
                size = pointsize - 1,
                face = "bold", colour = "black"
            ),
            legend.title = ggplot2::element_text(
                size = pointsize,
                face = "bold"
            ),
            panel.background = ggplot2::element_blank(),
            # panel.grid.minor = ggplot2::element_line(size = 0.25,
            # linetype = 'solid',colour = "grey"),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA, size = 1.2
            ),
            panel.grid.major = ggplot2::element_line(colour = "grey95")
        ) +
        ggplot2::scale_fill_gradient(
            low = "white",
            high = color, space = "Lab", na.value = "grey50",
            guide = "colourbar", aesthetics = "fill"
        ) +
        ggplot2::scale_y_discrete(limits = rev(levels(filtered.brinf$pw.name))))

} # bubblePlotPathwaysLR


#' Heatmap function for gene signature expression
#'
#' Generate one heatmap used by
#' \code{signatureHeatmaps}.
#'
#' @param counts  Matrix of counts.
#' @param name    Name of the heatmap to generate.
#' @param h.height  Heatmap  height in cm.
#' @param fontsize  Font size for row (gene) names.
#' @param col.fontsize  Font size for column (sample) names.
#' @param legend.fontsize  Font size for the legends.
#' @param annot.fontsize   Font size for column annotation names.
#' @param scoring Vector of sample scores for a
#' chosen pathway. If NULL, then no column annotation is produced.
#' @param cols.scoring   Fixed colorRamp2 object.
#' @param hcl.palette   Palette from
#' HCL colormaps supported by ComplexHeatmap.
#'
#' @return A ComplexHeatmap object.
#'
#' This is a convenience function that relies
#' on the \code{ComplexHeatmap}
#' package.
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @keywords internal
.customheatmap <- function(
    counts,
    name,
    fontsize = 6,
    col.fontsize = 6,
    legend.fontsize = 8,
    annot.fontsize = 8,
    h.height = 3,
    scoring = NULL,
    cols.scoring = NULL,
    hcl.palette = "Blues 3",
    show_column_names = FALSE) {
  
  counts.scaled <- t(scale(t(counts)))
  
  cols <- circlize::colorRamp2(
    breaks = c(-1, 0, 3),
    hcl_palette = hcl.palette, reverse = TRUE
  )
  
  if (!is.null(scoring)){
    top.annotation <- ComplexHeatmap::HeatmapAnnotation(
      LR.scoring = as.vector(scoring),
      border = c(LR.scoring = TRUE),
      show_legend = TRUE,
      simple_anno_size = grid::unit(2.5, "mm"),
      show_annotation_name = TRUE,
      col = list(LR.scoring = cols.scoring),
      annotation_name_gp= gpar(fontsize = annot.fontsize)
    )
  }
  
  di.gene <- stats::dist(counts.scaled)
  hc.gene <- stats::hclust(di.gene, method = "ward.D")
  dend.row <- stats::as.dendrogram(hc.gene)
  
  di.spl <- stats::dist(t(counts.scaled))
  hc.spl <- stats::hclust(di.spl, method = "ward.D")
  dend.spl <- stats::as.dendrogram(hc.spl)
  
  ComplexHeatmap::ht_opt(
    heatmap_border = TRUE,
    annotation_border = FALSE
  )
  
  if (!is.null(scoring)){
    ComplexHeatmap::Heatmap(counts.scaled,
                            name = name,
                            cluster_rows = dend.row, cluster_columns = dend.spl,
                            show_row_dend = FALSE, show_column_dend = TRUE,
                            col = cols, show_row_names = TRUE,
                            show_column_names = show_column_names,
                            use_raster = TRUE, raster_device = "png",
                            raster_quality = 8, raster_by_magick = FALSE,
                            rect_gp = grid::gpar(col = "white"),
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_names_gp = grid::gpar(fontsize = col.fontsize),
                            top_annotation = top.annotation,
                            height = grid::unit(h.height, "cm"),
                            column_gap = grid::unit(0.5, "mm"),
                            heatmap_legend_param = list(
                              labels_gp = gpar(fontsize = legend.fontsize)
                            )
    )
  }
  else{
    ComplexHeatmap::Heatmap(counts.scaled,
                            name = name,
                            cluster_rows = dend.row, cluster_columns = dend.spl,
                            show_row_dend = FALSE, show_column_dend = TRUE,
                            col = cols, show_row_names = TRUE,
                            show_column_names = show_column_names,
                            use_raster = TRUE, raster_device = "png",
                            raster_quality = 8, raster_by_magick = FALSE,
                            rect_gp = grid::gpar(col = "white"),
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_names_gp = grid::gpar(fontsize = col.fontsize),
                            height = grid::unit(h.height, "cm"),
                            column_gap = grid::unit(0.5, "mm"),
                            heatmap_legend_param = list(
                              labels_gp = gpar(fontsize = legend.fontsize)
                            )
    )
  }
  
} # .customheatmap


#' Heatmap function to dissect one pathway signature
#'
#' Plots a stack of three heatmaps to assess the expression of
#' the target genes or proteins in a chosen pathway, the receptor expressions,
#' and the ligand expressions.
#'
#' @param pathway        The chosen pathway name.
#' @param bsrdm     BulkSignalR data model object.
#' @param bsrsig     BulkSignalR signature object.
#' @param heights    A vector of 3 heights (in cm) for the 3 heatmaps.
#' @param fontsize    Font size for the gene names.
#' @param legend.fontsize  Font size for the legends.
#' @param title.fontsize   Font size for the pathway name as plot title.
#' @param col.fontsize  Font size for column (sample) names.
#' @param annot.fontsize   Font size for column annotation names.
#' @param ht.gap   Space between heatmaps (in mm).
#' @param show_column_names   Add column names in the heatmaps.

#' @return  A plot is created.
#'
#' @export
#' @examples
#' data(bsrdm, package = "BulkSignalR")
#' data(bsrinf, package = "BulkSignalR")
#' if(FALSE){
#' bsrinf.redP <- reduceToPathway(bsrinf)
#' bsrinf.redPBP <- reduceToBestPathway(bsrinf)
#' bsrsig.redPBP <- BSRSignature(bsrinf, qval.thres = 1)
#' pathway1 <- pathways(bsrsig.redPBP)[1]
#' signatureHeatmaps(
#' pathway = pathway1,
#' bsrdm = bsrdm,
#' bsrsig = bsrsig.redPBP,
#' show_column_names = TRUE)
#' }
#' @import ComplexHeatmap
#' @importFrom ComplexHeatmap %v%
#' @importFrom circlize colorRamp2
#' @import grid
signatureHeatmaps <- function(pathway,
                              bsrdm,
                              bsrsig,
                              heights = c(4,2,4),
                              fontsize = 6,
                              legend.fontsize = 8,
                              title.fontsize = 8,
                              col.fontsize = 6,
                              annot.fontsize = 8,
                              ht.gap = 3,
                              show_column_names = TRUE) {
  
  idx.path.sig <- which(pathways(bsrsig) == pathway)
  
  if (rlang::is_empty(idx.path.sig)) {
    stop("Pathway is not defined in signature.")
  }
  if (!is.numeric(heights)){
    stop("heights must be provided as a numerical vector.")
  }
  if (length(heights) != 3){
    stop("Three heights must be provided exactly.")
  }
  
  scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig)
  
  counts <- ncounts(bsrdm)
  
  filter.L <- unique(unlist(ligands(bsrsig)[idx.path.sig]))
  
  counts.L <- counts[filter.L, ]
  
  palette.L <- "RdPu"
  cols.L <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                 hcl_palette = palette.L, reverse = TRUE)
  
  filter.R <- unique(unlist(receptors(bsrsig)[idx.path.sig]))
  filter.T <- unique(unlist(tgGenes(bsrsig)[idx.path.sig]))
  
  # Remove in receptors, genes that are potential targets.
  filter.R <- filter.R[!filter.R %in% filter.T]
  
  counts.R <- counts[filter.R, ]
  palette.R <- "YlGn"
  cols.R <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                 hcl_palette = palette.R, reverse = TRUE)
  
  counts.T <- counts[filter.T, ]
  palette.T <- "Blues 3"
  cols.T <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                 hcl_palette = palette.T, reverse = TRUE)
  
  cols.scoring <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                       colors = c("blue", "white", "red"))
  
  abundance.samples <- dim(counts.T)[2]
  abundance.genes <- c(dim(counts.T)[1], dim(counts.R)[1], dim(counts.L)[1])
  
  if (length(idx.path.sig) == 1){
    scoring <- as.vector(scoresPathway[idx.path.sig[1],])
  }
  else{
    scoring <- colMeans(scoresPathway[idx.path.sig,])
  }
  
  p.T <- .customheatmap(
    counts = counts.T,
    name = "Targets",
    h.height = heights[1],
    scoring = scoring,
    hcl.palette = palette.T,
    cols.scoring = cols.scoring,
    show_column_names = FALSE,
    fontsize = fontsize,
    legend.fontsize = legend.fontsize,
    annot.fontsize = annot.fontsize
  )
  
  p.R <- .customheatmap(
    counts = counts.R,
    name = "Receptors",
    h.height = heights[2],
    hcl.palette = palette.R,
    show_column_names = FALSE,
    fontsize = fontsize,
    legend.fontsize = legend.fontsize,
  )
  
  p.L <- .customheatmap(
    counts = counts.L,
    name = "Ligands",
    h.height = heights[3],
    hcl.palette = palette.L,
    show_column_names = show_column_names,
    fontsize = fontsize,
    legend.fontsize = legend.fontsize,
    col.fontsize = col.fontsize,
  )
  
  ComplexHeatmap::draw(p.T %v% p.R %v% p.L,
                       ht_gap = grid::unit(ht.gap, "mm"),
                       column_title = pathway,
                       column_title_gp = gpar(fontsize=title.fontsize))
  
} # signatureHeatmaps


#' Heatmap function for LR scores
#'
#' Generate a heatmap representing ligand-receptor gene
#' signature scores.
#'
#' @param mat.c         A matrix with the signature scores such as output by
#' \code{scoreLRGeneSignatures()}.
#' @param dend.row       A precomputed row dendrogram.
#' @param dend.spl       A precompute sample (column) dendrogram.
#' @param cols           A vector of colors to use for the heatmap.
#' @param pointsize      Heatmap font point size
#' @param bottom.annotation  \code{ComplexHeatmap} package bottom annotations.
#' @param n.col.clust    Number of column clusters.
#' @param n.row.clust    Number of row clusters.
#' @param gap.size       Gap size between clusters.
#' @param cut.p          Proportion of top and bottom values for thresholding.
#' @param row.names      A logical to turn on/off the display of row names.
#' @param column.names      A logical to turn on/off the display
#' of column (sample) names.
#' @param hcl.palette    support for HCL colormaps in ComplexHeatmap
#' using color mapping function with circlize::colorRamp2().
#' palettes are listed in grDevides::hcl.pals().
#' of row (gene) names.
#' @param reverse    A logical to reverse or not colors in hcl.palette.
#'
#' @return A heatmap. Since heatmap plotting tend to be slow on the screen,
#' it is advisable to plot in a file instead.
#'
#' If hcl.palette is set, the colors parameter won't be used.
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
#' data(bsrdm, package = "BulkSignalR")
#' data(bsrinf, package = "BulkSignalR")
#' 
#' bsrinf.redBP <- reduceToBestPathway(bsrinf)
#' bsrsig.redBP <- BSRSignature(bsrinf,
#'     qval.thres = 0.001
#' )
#'
#' scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
#'     name.by.pathway = FALSE
#' )
#' simpleHeatmap(scoresLR[1:3, ],
#'     column.names = TRUE,
#'     hcl.palette = "Cividis")
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
simpleHeatmap <- function(mat.c,
    dend.row = NULL,
    dend.spl = NULL, cols = NULL, pointsize = 4,
    bottom.annotation = NULL, n.col.clust = 0,
    n.row.clust = 0, gap.size = 0.5, cut.p = 0.01, 
    row.names = TRUE,
    column.names = TRUE, hcl.palette = NULL,
    reverse = FALSE) {

    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop(
            "Package \"ComplexHeatmap\" needed for this function ",
            "to work. Please install it."
        )
    }
    if (cut.p < 0 || cut.p > 0.1) {
        stop("cut.p must lie in [0;0.1]")
    }

    if (cut.p != 0) {
        mat.c.cut <- .cutExtremeValues(mat.c, cut.p)
    }

    if (is.null(cols)) {
        if (!requireNamespace("circlize", quietly = TRUE)) {
            stop(
                "Package \"circlize\" needed for this function to ",
                "work (generation of color scale). Please install it."
            )
        }

        cols <- circlize::colorRamp2(
            breaks = c(min(mat.c.cut), 0, max(mat.c.cut)),
            colors = c("royalblue3", "white", "orange")
        )

        if (!is.null(hcl.palette)) {
            cols <- circlize::colorRamp2(
                breaks = c(min(mat.c.cut), 0, max(mat.c.cut)), ,
                hcl_palette = hcl.palette, reverse = reverse
            )
        }
    }

    if (is.null(dend.spl)) {
        di.spl <- stats::dist(t(mat.c))
        hc.spl <- stats::hclust(di.spl, method = "ward.D")
        dend.spl <- stats::as.dendrogram(hc.spl)
    }
    if (is.null(dend.row)) {
        di.gene <- stats::dist(mat.c)
        hc.gene <- stats::hclust(di.gene, method = "ward.D")
        dend.row <- stats::as.dendrogram(hc.gene)
    }

    
    if (n.row.clust > 0) {
        if (n.col.clust > 0) {
            plot(ComplexHeatmap::Heatmap(mat.c,
                cluster_rows = dend.row,
                cluster_columns = dend.spl, col = cols,
                show_row_names = row.names,
                show_column_names = column.names, 
                use_raster = TRUE, raster_device = "png",
                raster_quality = 8, raster_by_magick = FALSE,
                row_names_gp = grid::gpar(fontsize = pointsize),
                show_row_dend = TRUE, bottom_annotation = bottom.annotation,
                split = n.row.clust, gap = grid::unit(gap.size, "mm"),
                column_split = n.col.clust,
                column_gap = grid::unit(gap.size, "mm")
            ))
        } else {
            plot(ComplexHeatmap::Heatmap(mat.c,
                cluster_rows = dend.row,
                cluster_columns = dend.spl, col = cols,
                show_row_names = row.names,
                show_column_names = column.names, 
                use_raster = TRUE, raster_device = "png",
                raster_quality = 8, raster_by_magick = FALSE,
                row_names_gp = grid::gpar(fontsize = pointsize),
                show_row_dend = TRUE, bottom_annotation = bottom.annotation,
                split = n.row.clust, gap = grid::unit(gap.size, "mm")
            ))
        }
    } else if (n.col.clust) {
        plot(ComplexHeatmap::Heatmap(mat.c,
            cluster_rows = dend.row,
            cluster_columns = dend.spl, col = cols, show_row_names = row.names,
            show_column_names = column.names, 
            use_raster = TRUE, raster_device = "png",
            raster_quality = 8, raster_by_magick = FALSE,
            row_names_gp = grid::gpar(fontsize = pointsize),
            column_names_gp = grid::gpar(fontsize = pointsize),
            show_row_dend = TRUE, bottom_annotation = bottom.annotation,
            column_split = n.col.clust, column_gap = grid::unit(gap.size, "mm")
        ))
    } else {
        plot(ComplexHeatmap::Heatmap(mat.c,
            cluster_rows = dend.row,
            cluster_columns = dend.spl, col = cols, 
            show_row_names = row.names,
            show_column_names = column.names,
            row_names_gp = grid::gpar(fontsize = pointsize),
            column_names_gp = grid::gpar(fontsize = pointsize),
            use_raster = TRUE, raster_device = "png",
            raster_quality = 8, raster_by_magick = FALSE,
            show_row_dend = TRUE, bottom_annotation = bottom.annotation
        ))
    }

    
} # simpleHeatmap


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
#' such as CIBERSORT or BisqueRNA.
#' It is provided for convenience.
#' @return A matrix containing the scores of
#' each gene signature in each sample.
#' Note that ligand-receptor gene
#' signature scores should be computed with
#' \code{"\link[=BSRDataModel-class]{scoreLRGeneSignatures}"}
#' instead.
#' @export
#' @examples
#' data(sdc, package = "BulkSignalR")
#' data(bsrdm, package = "BulkSignalR")
#' 
#' data(immune.signatures, package = "BulkSignalR")
#' imm.scores <- scoreSignatures(bsrdm, immune.signatures)
#' @importFrom methods is
#' @importFrom matrixStats rowMeans2
scoreSignatures <- function(ds, ref.signatures, robust = FALSE) {
    if (!is(ds, "BSRDataModel")) {
        stop("ds must be BSRDataModel object")
    }

    ref.signatures <- ref.signatures[ref.signatures$gene %in%
        rownames(ncounts(ds)), ]
    pool <- unique(ref.signatures$gene)

    # convert normalized counts into z-scores
    if (logTransformed(ds)) {
        ncounts <- 2**ncounts(ds)[pool, ]
    } else {
        ncounts <- ncounts(ds)
    }
    if (robust) {
        z <- (ncounts - apply(ncounts, 1, stats::median)) 
        z <- z / apply(ncounts, 1, stats::mad)
    } else {
        z <- (ncounts - matrixStats::rowMeans2(ncounts)) / 
        apply(ncounts, 1, stats::sd)
    }

    # compute the gene signature scores
    pop <- unique(ref.signatures$signature)
    sig <- matrix(0,
        nrow = length(pop), ncol = ncol(ncounts),
        dimnames = list(pop, colnames(ncounts))
    )
    for (p in pop) {
        n <- sum(ref.signatures$signature == p)
        if (n > 1) {
            sig[p, ] <- as.numeric(colSums(
                z[ref.signatures$gene[ref.signatures$signature == p], ]
            ) / n)
        } else {
            sig[p, ] <- as.numeric(
                z[ref.signatures$gene[ref.signatures$signature == p], ]
            )
        }
    }

    sig
    
} # scoreSignatures


#' Alluvial plot
#'
#' @description Representation of the links
#' between ligands, receptors, and pathways.
#'
#' @param bsrinf A BSRInference object.
#' @param keywords vector of keywoprds to filter pathways.
#' @param type filter on Ligand, Receptor or pathway id.
#' @param qval.thres threshold over Q-value.
#' @return NULL
#'
#' This is a convenience function that relies on the \code{ggalluvial}
#' package to propose a simple way
#' of representing ligands, receptors,
#  and downstream pathway associations.
#' @import ggplot2
#' @import ggalluvial
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' alluvialPlot(bsrinf,
#'     keywords = c("LAMC1"),
#'     type = "L",
#'     qval.thres = 0.01)
alluvialPlot <- function(bsrinf, keywords, type = c("L", "R", "pw.id"),
                        qval.thres = 0.01) {
    interactions <- data.frame(
        L = unlist(ligands(bsrinf)),
        R = unlist(receptors(bsrinf)),
        pw.name = LRinter(bsrinf)$pw.name,
        pw.id = LRinter(bsrinf)$pw.id,
        qval = LRinter(bsrinf)$qval,
        targets = vapply(tgGenes(bsrinf),
            FUN = function(x) paste(x, collapse = "  "),
            character(1)
        )
    )

    type <- match.arg(type)

    if (type == "L") {
        subset.interactions <- interactions[interactions$L %in% keywords, ]
    }
    if (type == "R") {
        subset.interactions <- interactions[interactions$R %in% keywords, ]
    }
    if (type == "pw.id") {
        subset.interactions <- interactions[interactions$pw.id %in% keywords, ]
    }

    if (dim(subset.interactions)[1] == 0) {
        message(paste(keywords, collapse = " "),
            " for ", type, " not found.")
        stop("Try another value for filtering.")
    }
    subset.interactions <- subset.interactions[
    subset.interactions$qval <= qval.thres, ]
    subset.interactions$count <- 1
    subset.interactions <- subset.interactions[, 
    c("L", "R", "pw.name", "count")]

    stratum <- ggalluvial::StatStratum

    pl <- ggplot2::ggplot(
        subset.interactions,
        ggplot2::aes_string(y = "count", axis1 = "L", 
            axis2 = "R", axis3 = "pw.name")
    ) +
        ggalluvial::geom_alluvium(ggplot2::aes_string(fill = "R"),
            width = 1 / 12) +
        ggalluvial::geom_stratum(width = 1 / 12, 
            fill = "black", color = "grey") +
        ggplot2::geom_label(stat = stratum, 
            aes(label = ggplot2::after_stat(stratum))) +
        ggplot2::scale_x_discrete(limits = c("L", "R", "pw.name"),
            expand = c(0.5, 0.5)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::ggtitle("Ligand-Receptor Interactions & Underlying Pathways")

    pl <- pl + ggplot2::theme_bw()
    pl <- pl + ggplot2::theme(legend.position = "none")
    pl <- pl + ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = .5),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
    )

    plot(pl)

} # alluvialPlot


#' Chord Diagram of LR interactions with correlations
#'
#' @description Chord diagram.
#'
#' @param bsrinf A BSRInference object
#' @param pw.id.filter One Pathway ID accepted only to
#  retrieve the respective LR interactions.
#' @param qval.thres Threshold over Q-values.
#' @param ligand Ligand
#' of the LR pair that you want to
#' highlight in the chord diagram.
#' @param receptor Receptor
#' of the LR pair that you want to highlight
#' in the chord diagram.
#' @param limit Number of interactions you can visualize.
#  Maximum set to 30.
#' @return Circos Plot on the screen or a file
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2 circos.par chordDiagramFromDataFrame
#' @importFrom circlize circos.trackPlotRegion circos.text circos.axis
#' @importFrom circlize get.cell.meta.data circos.clear
#'
#' @export
#' @examples
#' data(bsrinf, package = "BulkSignalR")
#' chordDiagramLR(bsrinf,
#' pw.id.filter = "R-HSA-3000178",
#' limit = 20,
#' ligand="ADAM15", 
#' receptor="ITGAV"
#' )
chordDiagramLR <- function(
    bsrinf,
    pw.id.filter = NULL, qval.thres = 1,
    ligand = NULL, receptor = NULL,
    limit = 20) {

    if (limit > 40) {
        message("Number of selected interactions is too large", limit, ".")
        message("Number of visualised interactions sould be less than 40.")
    }

    if (!is.null(ligand) && !is.null(receptor)) {
        pair.to.highlight <- paste(ligand, receptor, sep = "-")
    } else {
        pair.to.highlight <- NULL
    }

    dataframe.bsrinf <- data.frame(
        ligands = unlist(ligands(bsrinf)),
        receptors = unlist(receptors(bsrinf)),
        corr = LRinter(bsrinf)$LR.corr,
        pw.id = LRinter(bsrinf)$pw.id,
        pathways = LRinter(bsrinf)$pw.name,
        pval = LRinter(bsrinf)$pval,
        qval = LRinter(bsrinf)$qval
    )

    dataframe.bsrinf$pair <- paste(dataframe.bsrinf$ligands,
        dataframe.bsrinf$receptors,
        sep = "-"
    )

    if (!is.null(pair.to.highlight) && 
        !(pair.to.highlight %in% dataframe.bsrinf$pair)) {
        stop(
            "Highlighted LR pair ", pair.to.highlight, " was not found for ",
            pw.id.filter, ".\n"
        )
    }

    # Filters
    if (!is.null(pw.id.filter)) {
        dataframe.bsrinf <- dataframe.bsrinf[
        dataframe.bsrinf$pw.id %in% pw.id.filter, ]
    }

    dataframe.bsrinf <- dataframe.bsrinf[
    dataframe.bsrinf$qval < qval.thres, ]


    if (dim(dataframe.bsrinf)[1] == 0) {
        stop("Pathway ID was not found.\n")
    }

    if (dim(dataframe.bsrinf)[1] < limit) {
        limit <- dim(dataframe.bsrinf)[1]
        message("Only ", limit, " interactions were found.\n")
    }

    dataframe.bsrinf <- dataframe.bsrinf[order(dataframe.bsrinf$qval), ]
    dataframe.bsrinf <- dataframe.bsrinf[seq_len(limit), ]
    dataframe.bsrinf <- unique(dataframe.bsrinf[,
        c("ligands", "receptors", "corr", "pair")])

    if (length(unique(dataframe.bsrinf$corr)) == 1){
        cr <- "#febd17" # for BSRInferenceComp class with all corr = 1
    }
    else{
        cr <- circlize::colorRamp2(c(
            min(dataframe.bsrinf$corr),
            max(dataframe.bsrinf$corr)
        ), c("white", "#febd17"))
    }

    myList.ligands <- rep("gray25", 
        times = length(dataframe.bsrinf$ligands))
    names(myList.ligands) <- as.list(dataframe.bsrinf$ligands)

    myList.receptors <- rep("#7fbb00", 
        times = length(dataframe.bsrinf$receptors))
    names(myList.receptors) <- as.list(dataframe.bsrinf$receptors)

    myList <- c(myList.receptors, myList.ligands)

    link.col <- rep("dodgerblue3", nrow(dataframe.bsrinf))

    link.lwd <- rep(1, nrow(dataframe.bsrinf))

    link.width <- rep(0.12, nrow(dataframe.bsrinf))

    if (!is.null(pair.to.highlight)) {
        index.filter <- which(pair.to.highlight == dataframe.bsrinf$pair)
        link.col[index.filter] <- "#d40000"
        link.lwd[index.filter] <- 3
        link.width[index.filter] <- 0.15
    }


    interactions <- data.frame(
        from = dataframe.bsrinf$ligands,
        to = dataframe.bsrinf$receptors,
        value = dataframe.bsrinf$corr
    )

    circlize::circos.par(points.overflow.warning = FALSE)

    circlize::chordDiagramFromDataFrame(interactions,
        # grid.col = grid.col,
        col = cr,
        annotationTrack = "grid",
        grid.col = myList,
        transparency = 0.7,
        preAllocateTracks = 1,
        directional = 1,
        direction.type = "arrows",
        link.arr.length = link.width,
        link.arr.width = link.width,
        link.arr.type = "triangle",
        link.arr.lty = "solid",
        link.arr.lwd = link.lwd,
        link.arr.col = link.col,
        big.gap = 2,
        small.gap = 1
    )

    circlize::circos.trackPlotRegion(
        track.index = 2,
        panel.fun = function(x, y) {
            xlim <- circlize::get.cell.meta.data("xlim")
            ylim <- circlize::get.cell.meta.data("ylim")
            sector.name <- circlize::get.cell.meta.data("sector.index")

            circlize::circos.text(mean(xlim), ylim[1] + 1.9, sector.name,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.7
            )

            circlize::circos.axis(
                h = "top", labels = FALSE, minor.ticks = FALSE,
                major.tick.length = 1,
                major.at = c(xlim),
                sector.index = sector.name,
                track.index = 2
            )
        }
    )

    # LEGEND

    lgd_points <- Legend(
        labels = c("Ligands", "Receptors"),
        type = "points", pch = 16,
        legend_gp = grid::gpar(col = c("gray25", "#7fbb00")),
        title_position = "topleft",
        labels_gp = grid::gpar(font = 6),
        title = "LR"
    )

    lgd_links <- Legend(
        at = c(round(min(interactions$value),
            digits = 2
        ), round(max(interactions$value),
            digits = 2
        )), col_fun = cr,
        title = "Correlation", direction = "horizontal",
        grid_width = unit(0.9, "mm"),
        grid_height = unit(1.3, "mm"),
        labels_gp = grid::gpar(font = 6),
        title_position = "topcenter",
    )

    lgd_list_vertical <- packLegend(lgd_points)

    draw(lgd_list_vertical,
        x = unit(2, "mm"),
        y = unit(2, "mm"),
        just = c("left", "bottom")
    )
    draw(lgd_links,
        x = unit(2.7, "inch"),
        y = unit(2, "mm"),
        just = c("left", "bottom")
    )

    circlize::circos.clear()


} # chordDiagramLR
