#' Get correlated ligand-receptor pairs.
#'
#' Internal function to compute the Spearman correlations
#' of all the ligand-receptor
#' pairs in LRdb and return those above a minimum value.
#'
#' @param ds              A BSRDataModel object.
#' @param min.cor         The minimum correlation required.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#'
#' @return A data frame containing putative ligand-receptor pairs along
#'   with their correlations above \code{min.cor}. This table is the first step
#'   of a ligand-receptor analysis.
#'
#'
#' @details The \code{restrict.genes} parameter is used for special cases where
#'   LRdb must be further restricted to a subset.
#'   The putative ligand-receptor pairs has 3 columns : R, L and corr.
#'
#' @importFrom methods is
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @keywords internal
.getCorrelatedLR <- function(ds, min.cor = 0.25, restrict.genes = NULL) {

    # local binding
    i <- NULL

    if ((min.cor < -1) || (min.cor > 1))
        stop("min.cor must lie in [-1;+1]")
    if (!is(ds, "BSRDataModel"))
        stop("ds must be an object of class BSRDataModel")

    lrgenes <- intersect(c(LRdb$ligand,
                           LRdb$receptor), rownames(ncounts(ds)))
    if (!is.null(restrict.genes))
        lrgenes <- intersect(lrgenes, restrict.genes)

    # compute all the correlations at once
    corlr <- stats::cor(t(ncounts(ds)[lrgenes, ]), method = "spearman")

    # get the pairs
    pairs <- foreach::foreach(i = seq_len(nrow(LRdb)),
                              .combine = rbind) %do% {
              if (LRdb$ligand[i] %in% rownames(corlr) &&
                  LRdb$receptor[i] %in% rownames(corlr))
                  data.frame(L = LRdb$ligand[i],
                             R = LRdb$receptor[i],
                             corr = corlr[LRdb$ligand[i],
                                          LRdb$receptor[i]],
                             stringsAsFactors = FALSE)
              else
                  NULL
    }

    good <- pairs$corr >= min.cor
    pairs[good,]

}  # .getCorrelatedLR


#' Internal function to check receptor signaling downstream
#'
#' @param lr              A data frame as returned by
#'   \code{.getCorrelatedLR()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway
#'   IDs).
#' @param rncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param id.col          Column index or name in \code{pw} for the pathway IDs.
#' @param gene.col        Column index or name in \code{pw} for the gene
#'   symbols.
#' @param pw.col          Column index or name in \code{pw} for the pathway
#'   names.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#'
#' @return A table reporting all the ligand-receptor pairs provided in \code{lr}
#'   along with the pathways found and data about target gene correlations with
#'   the receptor.
#'
#' @importFrom foreach %do% %dopar%
#' @keywords internal
.downstreamSignaling <- function(lr, pw, pw.size, rncounts, id.col, gene.col,
                                 pw.col, min.positive, with.complex = TRUE) {
    if (!is.matrix(rncounts))
        stop("rncounts must be a matrix")

    # local binding
    r <- p <- pl <- id <- NULL

    # define interaction types
    control.int <- "controls-expression-of"
    incomplex.int <- c("in-complex-with","interacts-with")
    directed.int <- c("controls-state-change-of", "catalysis-precedes",
                      "controls-expression-of", "controls-transport-of",
                      "controls-phosphorylation-of")
    if (with.complex)
        correlated.int <- union(control.int, incomplex.int)
    else
        correlated.int <- control.int

    # compute downstream correlations
    corrg <- stats::cor(t(rncounts), method = "spearman")

    # the global computation above is faster than restricted to the receptors
    corrg <- corrg[unique(lr$R), ]

    # check each putative LR pair, loop over the receptors
    reg.proc <- foreach::foreach(r=unique(lr$R),.combine=rbind) %do% {
        # reg.proc <- NULL
        # for (r in unique(lr$putative.pairs$R)){

        # loop over the pathways containing the receptor r
        pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
        if (length(pa)>0){
            receptor.ligands <- unique(lr$L[lr$R==r])
            best.2nd <- foreach::foreach(p=pa,.combine=rbind) %do% {
                # best.2nd <- NULL
                # for (p in pa){
                int <- SingleCellSignalR::PwC_ReactomeKEGG[
                    SingleCellSignalR::PwC_ReactomeKEGG$a.gn %in% pw[pw[[id.col]]==p,gene.col] &
                    SingleCellSignalR::PwC_ReactomeKEGG$b.gn %in% pw[pw[[id.col]]==p,gene.col],
                ]
                directed <- int$type %in% directed.int

                # double the undirected interactions and generate a directed graph
                ret <- int[!directed,c("a.gn", "b.gn")]
                from <- ret$a.gn
                ret$a.gn <- ret$b.gn
                ret$b.gn <- from
                d.int <- unique(rbind(int[,c("a.gn", "b.gn")],ret))
                g <- igraph::graph_from_data_frame(d.int, directed=TRUE)

                # extract the target genes of receptor r
                if (r %in% d.int$a.gn || r %in% d.int$b.gn){
                    # putative targets in the pathway
                    target.genes <- setdiff(c(
                        int[int$type %in% correlated.int & int$a.gn==r, "b.gn"],
                        int[int$type %in% correlated.int & int$b.gn==r, "a.gn"],
                        int[int$type %in% directed.int, "b.gn"]),
                        r
                    )

                    # reduce putative to reachable from the receptor
                    sp <- igraph::shortest.paths(g, r, target.genes)
                    target.genes <- colnames(sp)[!is.infinite(sp[r,])]

                    # eliminate ligands of the receptor if present
                    target.genes <- setdiff(target.genes, receptor.ligands)

                    if (length(target.genes) >= min.positive){
                        # if all conditions are met, list all target genes with
                        # their correlations to the receptor in a data frame
                        # row. Target genes are sorted wrt correlations.
                        c <- corrg[r, target.genes]
                        o <- order(c)
                        c <- c[o]
                        target.genes <- target.genes[o]
                        data.frame(pathway=p, target.corr=paste(c,collapse=";"),
                                   target.genes=paste(target.genes,collapse=";"),
                                   len=length(c),
                                   stringsAsFactors=FALSE)
                    }
                    else
                        NULL
                }
                else
                    NULL
            }
            if (!is.null(best.2nd))
                # one or several pathways containing the receptor r were found,
                # combine them in |-separated strings
                data.frame(R=r, pathways=paste(best.2nd$pathway, collapse="|"),
                           target.corr=paste(best.2nd$target.corr, collapse='|'),
                           target.genes=paste(best.2nd$target.genes, collapse='|'),
                           len=paste(best.2nd$len, collapse='|'),
                           stringsAsFactors=FALSE)
            else
                NULL
        }
        else
            NULL
    }

    # combine LR pair correlations with R-target gene correlations
    rownames(reg.proc) <- reg.proc$R
    conf.pairs <- lr[lr$R %in% reg.proc$R,]
    conf.pairs$pwid <- reg.proc[conf.pairs$R, "pathways"]
    conf.pairs$target.corr <- reg.proc[conf.pairs$R, "target.corr"]
    conf.pairs$len <- reg.proc[conf.pairs$R, "len"]
    conf.pairs$corr.pw <- reg.proc[conf.pairs$R, "corr"]
    conf.pairs$target.genes <- reg.proc[conf.pairs$R, "target.genes"]
    pw.name <- unique(pw[,c(id.col, pw.col)])
    pw2name <- stats::setNames(pw.name[[2]], pw.name[[1]])
    conf.pairs$pwname <- foreach::foreach(pl=conf.pairs$pwid,.combine=c) %do% {
        paste(foreach::foreach(id=unlist(strsplit(pl,"\\|")),
                               .combine=c) %do% {
                                   pw2name[id]
                               },
              collapse="|"
        )
    }

    conf.pairs[,c("L", "R", "corr", "pwid", "pwname", "len", "target.genes",
                  "target.corr")]

}  # .downstreamSignaling


#' Internal function to check receptor signaling downstream
#'
#' Assess the existence of correlations between a receptor,
#' part of a ligand-receptor pair, and
#' genes coding for proteins forming a complex with the receptor or genes
#' regulated by the receptor downstream signaling.
#'
#' @param ds              A BSRDataModel object.
#' @param lr              A table as returned by \code{.getCorrelatedLR()}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @return A data frame extending \code{lr} content with the pathways found to
#' contain the receptors and data about target gene correlations with those
#' receptors. Strings in semi-colon-separated format are used to report
#' target genes and their Spearman correlations with the receptor in the
#' data frame. The target genes are sorted according to the correlation
#' coefficient.
#'
#' In a pathway of the reference, i.e., a Reactome pathway or the genes of a
#' GOBP term, the target genes are the
#' genes coding for proteins forming a complex with the receptor and the
#' genes in the pathway downstream the receptor,
#' which are given as regulated by the pathway. If \code{with.complex} is
#' set to \code{FALSE}, then only the
#' regulated genes are considered. Participation to a complex and being
#' regulated as well as the pathway directed topologies
#' are defined by Reactome and KEGG pathways as provided by PathwayCommons.
#'
#' The maximum pathway size is used to limit the redundancy inherent to GOBP
#' and Reactome. The minimum pathway size is
#' used to avoid overspecific, noninformative results.
#'
#' @importFrom methods is
#' @keywords internal
.checkReceptorSignaling <- function(ds, lr, reference=c("REACTOME-GOBP",
                                                     "REACTOME","GOBP"),
                                    max.pw.size=200, min.pw.size=5,
                                    min.positive=4, restrict.pw=NULL,
                                    with.complex=TRUE){

    if (!is(ds, "BSRDataModel"))
        stop("ds must be a BSRDataModel object")

    reference <- match.arg(reference)
    results <- list()

    # Reactome pathways
    if (reference %in% c("REACTOME-GOBP","REACTOME")){
        react <- reactome[reactome$`Gene name` %in% rownames(ncounts(ds)),]
        if (!is.null(restrict.pw))
            react <- react[react$`Reactome ID` %in% restrict.pw,]
        pw.size <- table(react$`Reactome ID`)
        pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
        contains.receptors <- react[react$`Gene name` %in% lr$R, "Reactome ID"]
        pw.size <- pw.size[names(pw.size) %in% contains.receptors]
        corgenes <- unique(c(lr$R,
                react[react$`Reactome ID` %in% names(pw.size), "Gene name"])
        )
        results$reactome.pairs <- .downstreamSignaling(lr, react, pw.size,
                ncounts(ds)[corgenes,], "Reactome ID", "Gene name",
                "Reactome name", min.positive, with.complex=with.complex)
    }

    # GOBP
    if (reference %in% c("REACTOME-GOBP","GOBP")){
        go <- gobp[gobp$`Gene name` %in% rownames(ncounts(ds)),]
        if (!is.null(restrict.pw))
            go <- go[go$`GO ID` %in% restrict.pw,]
        pw.size <- table(go$`GO ID`)
        pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
        contains.receptors <- go[go$`Gene name` %in% lr$R, "GO ID"]
        pw.size <- pw.size[names(pw.size) %in% contains.receptors]
        corgenes <- unique(c(lr$R,
                        go[go$`GO ID` %in% names(pw.size), "Gene name"])
        )
        results$gobp.pairs <- .downstreamSignaling(lr, go, pw.size,
                ncounts(ds)[corgenes,], "GO ID", "Gene name", "GO name",
                min.positive, with.complex=with.complex)
    }

    # merge
    if (reference == "REACTOME-GOBP"){
        pairs <- unique(rbind(results$reactome.pairs[,1:2],
                              results$gobp.pairs[,1:2])
        )
        react.keys <- paste(results$reactome.pairs[[1]],
                            results$reactome.pairs[[2]], sep="|")
        gobp.keys <- paste(results$gobp.pairs[[1]],
                           results$gobp.pairs[[2]], sep="|")
        results$merged.pairs <- rbind(results$reactome.pairs,
                    results$gobp.pairs[!(gobp.keys %in% react.keys),]
        )
    }
    else if (reference == "REACTOME")
        results$merged.pairs <- results$reactome.pairs
    else
        results$merged.pairs <- results$gobp.pairs

    results$merged.pairs

} # .checkReceptorSignaling


#' Internal function to assign P-values to LR interactions
#'
#' Estimate the P-value of each ligand-receptor pair based
#' on the data frame output by \code{\link{.checkReceptorSignaling}}.
#'
#' @param pairs         A data frame output by \code{checkReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @return A BSRInference object.
#' @keywords internal
.pValuesLR <- function(pairs, param, rank.p = 0.75,
                      fdr.proc = c("BH", "Bonferroni", "Holm", "Hochberg",
                                   "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {

    if (rank.p < 0 || rank.p > 1)
        stop("rank.p must lie in [0;1]")
    fdr.proc <- match.arg(fdr.proc)

    # prepare the chosen model CDF
    LR.par <- param$LR.0$model
    RT.par <- param$RT.0$model
    if (LR.par$distrib != RT.par$distrib)
        stop("Distinct statistical models for LR and RT nulls are not allowed")
    if (LR.par$distrib == 'censored_normal')
        cdf <- .cdfGaussian
    else if (LR.par$distrib == 'censored_mixed_normal')
        cdf <- .cdfMixedGaussian
    else if (LR.par$distrib == 'censored_stable')
        cdf <- .cdfAlphaStable
    else if (LR.par$distrib == 'empirical')
        cdf <- .cdfEmpirical
    else if (LR.par$distrib == 'kernel_empirical')
        cdf <- .cdfKernelEmpirical
    else
        stop(paste0("Unknown statistical model: ", LR.par$LR.0$model$distrib))

    # estimate P-values
    res <- NULL
    for (i in 1:nrow(pairs)){
        # all the data related to each pathway containing a given
        # receptor were collapsed separated by |
        # we need to split those pathways
        pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
        pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
        tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
        spear <- unlist(strsplit(pairs$target.corr[i],split="\\|"))
        len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))

        # estimate the LR correlation P-value
        if (pairs$corr[i] >= 0)
            # normal case
            p.lr <- 1 - cdf(pairs$corr[i], LR.par)
        else
            # to enable searching for inhibitory L-R interactions
            p.lr <- cdf(pairs$corr[i], LR.par)

        # estimate the target gene correlation P-value based on rank statistics
        # for the individual correlation Gaussian model
        for (k in 1:length(len)){
            spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
            r <- min(max(1,trunc(rank.p*len[k])),len[k])
            rank.corr <- spears[r]
            p.rt <- stats::pbinom(r-1, len[k], cdf(rank.corr, RT.par))
            res <- rbind(res,data.frame(pairs[i,c("L","R")],
                        LR.corr=pairs[i,"corr"], pw.id=pwid[k],
                        pw.name=pwname[k], rank=r, len=len[k],
                        rank.corr=rank.corr, target.genes=tg[k],
                        target.corr=spear[k], pval=p.lr*p.rt,
                        stringsAsFactors=FALSE))
        }
    }
    
    # avoid the impossible
    key <- paste(res$L, res$R, res$pw.id, sep="||")
    bad <- duplicated(key)
    res <- res[!bad,]

    # multiple hypothesis correction
    rawp <- res$pval
    adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
    res$qval <- adj$adjp[order(adj$index),fdr.proc]

    res

}  # .pValuesLR
