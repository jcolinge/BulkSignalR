#' Assign cell types to L-R interactions
#'
#' Generate a data.frame linking interactions to cell types.
#'
#' @param bsrdm   A BSRDataModel object.
#' @param bsrinf  A BSRInference object.
#' @param ct.scores   A matrix of cell type signature scores.
#' @param normalize.scores  A logical indicating whether scores should be
#' normalized before assigning cell types.
#' @param min.weight  Minimum weight to keep in the linear model (cell types
#' with lower weights will be discarded) if \code{lasso==TRUE}. Otherwise,
#' minimum correlation coefficient of each individual cell type.
#' @param min.r2  Minimum r2 between a candidate cell type and a L-R gene
#' signature score.
#' @param min.r2.after  Minimum r2 between the proposed linear model and
#' a L-R gene signature score to retain the model.
#' @param lasso  Logical indicating that the LASSO (or linear regression if
#' only one cell type satisfies the \code{min.r2} criterion) should be
#' used. Otherwise, Spearman linear correlation is used.
#' @param qval.thres  Maximum Q-value of the L-R pairs to be considered.
#' 
#' @return A data.frame containing the cell type assignments for each
#' L-R interaction. Unique interactions are considered only (thanks to
#' \code{"\link[=BSRInference-class]{reduceToBestPathway}"} that is applied
#' internally). An interaction can be associated with several cell types
#' or none. In case it is associated with a single cell type, it is labelled
#' autocrine (indicative only).
#' 
#' Cell type signature scores must be provided. They can be computed
#' with BulkSignalR utility function \code{"\link{scoreSignatures}"}, but
#' also any other external tool such as CIBERSORT or BisqueRNA. In case
#' such a tool would score cell types in a nonlinear fashion, we
#' recommend to transform the score matrix to restore a linear relationship
#' cell type abundance/score. By default, cell type (and L-R gene
#' signature) scores are normalized between 0 and 1 to make the weights of each
#' cel type in the linear models as comparable as possible.
#' @export
#' @examples
#' print('assignCellTypesToInteractions')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' 
#' # microenvironment cell populations
#' data(immune.signatures, package="BulkSignalR")
#' immune.signatures <- immune.signatures[immune.signatures$signature %in%
#'                                c("B cells","Dentritic cells","Macrophages",
#'                                  "NK cells","T cells","T regulatory cells"),]
#' data("tme.signatures", package="BulkSignalR")
#' signatures <- rbind(immune.signatures,tme.signatures[
#'           tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
#' tme.scores <- scoreSignatures(bsrdm, signatures)
#' 
#' # assignment
#' lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
#' @import glmnet
#' @importFrom foreach %do%
assignCellTypesToInteractions <- function(bsrdm, bsrinf, ct.scores,
                                          normalize.scores=TRUE, min.weight=0.1,
                                          min.r2=0.25, min.r2.after=0.35,
                                          lasso=TRUE, qval.thres=1e-3){
  
  # local binding
  i <- NULL

  # L-R interaction table and scores
  bsrinf.red <- reduceToBestPathway(bsrinf)
  inter <- LRinter(bsrinf.red)
  inter <- inter[inter$qval <= qval.thres,]
  bsrsig <- getLRGeneSignatures(bsrinf.red, qval.thres=qval.thres)
  lr.scores <- scoreLRGeneSignatures(bsrdm, bsrsig)
  
  # check order
  inter.keys <- paste(inter[["L"]], inter[["R"]], sep="-")
  score.keys <- gsub("\\{", "", rownames(lr.scores))
  score.keys <- gsub("\\}", "", score.keys)
  score.keys <- gsub(" / ", "-", score.keys)
  if (!all(inter.keys == score.keys))
    stop("LT table and scores not in the same order")
  
  # CT and LR score matrices preparation
  x <- t(ct.scores)
  if (normalize.scores){
    mins <- apply(x, 2, min)
    x <- sweep(x, 2, mins, "-")
    maxs <- apply(x, 2, max)
    x <- sweep(x, 2, maxs, "/")
    
    mins <- apply(lr.scores, 1, min)
    lr.scores <- sweep(lr.scores, 1, mins, "-")
    maxs <- apply(lr.scores, 1, max)
    lr.scores <- sweep(lr.scores, 1, maxs, "/")
  }
  
  # assignment of CTs
  a <- foreach::foreach(i=seq_len(nrow(inter)), .combine=rbind) %do% {
    y <- lr.scores[i,]
    c <- suppressWarnings(stats::cor(y, x, method="spearman"))
    c2 <- c**2
    good <- c2 >= min.r2 & c > 0
    if (sum(good) > 0)
      if (lasso){
        # actual LASSO since >1 CTs
        if (sum(good) > 1){
          xL <- x[,good]
          lb <- rep(0, length(colnames(xL)))
          ub <- rep(Inf, length(colnames(xL)))
          cvfit <- glmnet::cv.glmnet(xL, y, grouped=FALSE,
                             lower.limits=lb, upper.limits=ub)
          lasso.fit.l <- glmnet::glmnet(x=xL, y=y, lambda=cvfit$lambda.min,
                                lower.limits=lb, upper.limits=ub)
          MAE <- as.vector(glmnet::assess.glmnet(lasso.fit.l, newx=xL, newy=y)$mae)
          p <- stats::predict(lasso.fit.l, newx=xL)
          gc2 <- suppressWarnings(stats::cor(y,p))**2
          lcoef <- stats::coef(lasso.fit.l)[-1,1]
          if (gc2 > min.r2.after && sum(lcoef>0) > 0)
            data.frame(L=inter$L[i], R=inter$R[i], cell.type=names(lcoef)[lcoef>0],
                       weight=lcoef[lcoef>0], MAE=MAE, r2=gc2[1,1],
                       alg="LASSO", stringsAsFactors=FALSE)
          else
            NULL
        }
        else{
          # linear model since 1 CT
          dat <- data.frame(y=y, x=x[,good])
          lmfit <- stats::lm(y~x, data=dat)
          p <- stats::predict(lmfit)
          gc2 <- suppressWarnings(stats::cor(y,p))**2
          lcoef <- stats::coef(lmfit)[-1]
          if (gc2 > min.r2.after && lcoef > 0)
            data.frame(L=inter$L[i], R=inter$R[i], cell.type=colnames(x)[good],
                       weight=lcoef, MAE=mean(abs(p-y)), r2=gc2, alg="lm",
                       stringsAsFactors=FALSE)
          else
            NULL
        }
      }
    else
      # Spearman correlation only
      data.frame(L=inter$L[i], R=inter$R[i], cell.type=colnames(good)[good], r2=c2[good],
                 alg="Spearman.cor", stringsAsFactors=FALSE)
    else
      NULL
  }
  
  a[a$weight >= min.weight, ]
  rownames(a) <- NULL
  a
  
} # assignCellTypesToInteractions


#' Build a table describing a cellular network
#'
#' Generate a data.frame including all the links between cell types
#' mediated by L-R interactions with their respective weights.
#'
#' @param lr  The data.frame output by \code{"\link{assignCellTypesToInteractions}"}.
#' @param autocrine  A logical indicating whether autocrine interactions should
#' be included.
#' 
#' @return A data.frame containing all the links in the cellular network. A
#' link is created between two cell types as soon as there was a L-R interaction
#' that was associated with both cell types. The link is given a score
#' equal to the geometric mean of each cell type assignment r2.
#' @export
#' @examples
#' print('assignCellTypesToInteractions')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' 
#' # microenvironment cell populations
#' data(immune.signatures, package="BulkSignalR")
#' immune.signatures <- immune.signatures[immune.signatures$signature %in%
#'                                c("B cells","Dentritic cells","Macrophages",
#'                                  "NK cells","T cells","T regulatory cells"),]
#' data("tme.signatures", package="BulkSignalR")
#' signatures <- rbind(immune.signatures,tme.signatures[
#'           tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
#' tme.scores <- scoreSignatures(bsrdm, signatures)
#' 
#' # assignment
#' lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
#' 
#' # cellular network
#' g.table <- cellularNetworkTable(lr2ct)
#' @importFrom foreach %do%
cellularNetworkTable <- function(lr, autocrine=FALSE){
  
  # global binding
  key <- j <- NULL

  # paracrine
  keys <- paste(lr$L, lr$R, sep="-")
  t <- table(keys)
  if (sum(t > 1) > 1)
    nt <- foreach::foreach (key=names(t)[t>1], .combine=rbind) %do% {
      rows <- which(keys == key)
      L <- lr$L[rows[1]]
      R <- lr$R[rows[1]]
      r2 <- lr$r2[rows[1]]
      foreach::foreach (i=seq_len(length(rows)-1), .combine=rbind) %do% {
        row.1 <- rows[i]
        CT1 <- lr$cell.type[row.1]
        w.1 <- lr$weight[row.1]
        foreach::foreach (j=i+seq_len(length(rows)-i), .combine=rbind) %do% {
          row.2 <- rows[j]
          CT2 <- lr$cell.type[row.2]
          w.2 <- lr$weight[row.2]
          data.frame(L=L, R=R, inter=key, CT1=CT1, CT2=CT2, r2=r2, weight.1=w.1,
                     weight.2=w.2, score=sqrt(w.1*w.2), type="paracrine",
                     stringsAsFactors=FALSE)
        }
      }
    }
  else
    nt <- NULL
  
  # paracrine
  if (autocrine){
    if (sum(t == 1) > 1)
      nt <- rbind(nt, foreach::foreach (key=names(t)[t==1], .combine=rbind) %do% {
        i <- which(keys == key)
        L <- lr$L[i]
        R <- lr$R[i]
        r2 <- lr$r2[i]
        CT1 <- lr$cell.type[i]
        w.1 <- lr$weight[i]
        data.frame(L=L, R=R, inter=key, CT1=CT1, CT2=CT1, r2=r2, weight.1=w.1,
                   weight.2=w.1, score=w.1, type="autocrine",
                   stringsAsFactors=FALSE)
      })
  }
  
  nt
  
} # cellularNetworkTable


#' Build a cellular network
#'
#' Generate a igraph object including all the links between cell types.
#'
#' @param tab  The data.frame output by \code{"\link{cellularNetworkTable}"}.
#' 
#' @return A igraph object containing all the links in the cellular network.
#' @export
#' @examples
#' print('assignCellTypesToInteractions')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' 
#' # microenvironment cell populations
#' data(immune.signatures, package="BulkSignalR")
#' immune.signatures <- immune.signatures[immune.signatures$signature %in%
#'                                c("B cells","Dentritic cells","Macrophages",
#'                                  "NK cells","T cells","T regulatory cells"),]
#' data("tme.signatures", package="BulkSignalR")
#' signatures <- rbind(immune.signatures,tme.signatures[
#'           tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
#' tme.scores <- scoreSignatures(bsrdm, signatures)
#' 
#' # assignment
#' lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
#' 
#' # cellular network
#' g.table <- cellularNetworkTable(lr2ct)
#' gCN <- cellularNetwork(g.table)
#' #plot(gCN, edge.width=5*E(gCN)$score)
#' @import igraph
cellularNetwork <- function(tab){
  
  df <- tab[,c("CT1","CT2","inter","L","R","r2","score","type")]
  igraph::graph_from_data_frame(df, directed=FALSE)
  
} # cellularNetwork


#' Build a summary cellular network
#'
#' Generate a igraph object with one link between each cell type.
#'
#' @param tab  The data.frame output by \code{"\link{cellularNetworkTable}"}.
#' 
#' @return A igraph object containing a summary cellular network with
#' edge weights proportional to the sum of individual link scores. Edge
#' weight are normalized to a total of one.
#' @export
#' @examples
#' print('assignCellTypesToInteractions')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' 
#' # microenvironment cell populations
#' data(immune.signatures, package="BulkSignalR")
#' immune.signatures <- immune.signatures[immune.signatures$signature %in%
#'                                c("B cells","Dentritic cells","Macrophages",
#'                                  "NK cells","T cells","T regulatory cells"),]
#' data("tme.signatures", package="BulkSignalR")
#' signatures <- rbind(immune.signatures,tme.signatures[
#'           tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
#' tme.scores <- scoreSignatures(bsrdm, signatures)
#' 
#' # assignment
#' lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
#' 
#' # cellular network
#' g.table <- cellularNetworkTable(lr2ct)
#' gSummary <- summarizedCellularNetwork(g.table)
#' #plot(gSummary, edge.width=1+30*E(gSummary)$score)
#' @import igraph
#' @importFrom foreach %do%
summarizedCellularNetwork <- function(tab){
  i <- NULL
  sscore <- sum(tab$score)
  ct <- unique(tab[,c("CT1","CT2")])
  sum.tab <- foreach::foreach (i=seq_len(nrow(ct)),.combine=rbind) %do% {
    parsum <- sum(tab[tab$CT1==ct[i,"CT1"] & tab$CT2==ct[i,"CT2"], "score"])
    data.frame(CT1=ct[i,"CT1"], CT2=ct[i,"CT2"], score=parsum/sscore)
  }
  igraph::graph_from_data_frame(sum.tab, directed=FALSE)
  
} # summarizedCellularNetwork


# Relate to a gene set ========================================================


#' Relate ligands to a gene set
#'
#' Finds ligands related to a gene set by following receptor, and
#' receptor downstream pathway targets.
#'
#' @param bsrinf  BSRInference object.
#' @param gs   The gene set.
#' @param min.cor  Minimum Spearman correlation between the receptor of a triple
#' (L,R,pw) and a gene of the gene set.
#' @param qval.thres  Maximum Q-value imposed to the (L,R,pw) triples to be
#' considered.
#' @return A data.frame listing all the (L,R,pathway) triples that lead
#' to at least one gene in the gene set. The number of genes found
#' by each triple is indicated in the column \code{n.genes}.
#' 
#' @export
#' @examples
#' print('relateToGeneSet')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' data(p.EMT, package='BulkSignalR')
#' p.EMT <- p.EMT$gene
#' triggers <- relateToGeneSet(bsrinf, p.EMT) # Here in SDC instead of HNSCC!
#' @importFrom foreach %do%
relateToGeneSet <- function(bsrinf, gs, min.cor=0.25, qval.thres=0.001){
  
  i <- NULL
  # get target gens and L-R interactions
  inter <- LRinter(bsrinf)
  tg <- tGenes(bsrinf)
  tcor <- tgCorr(bsrinf)
  
  # impose threshold on FDR
  good <- inter$qval <= qval.thres
  inter <- inter[good,]
  tg <- tg[good]
  tcor <- tcor[good]
  
  # select L-R interactions that relate to the gene set
  foreach::foreach(i=seq_len(nrow(inter)), .combine=rbind) %do% {
    k <- (tg[[i]] %in% gs) & (tcor[[i]] >= min.cor)
    if (sum(k) > 0)
      data.frame(L=inter$L[i], R=inter$R[i], pw.name=inter$pw.name[i],
                 qval=inter$qval[i], n.genes=sum(k),
                 genes=paste(tg[[i]][k], collapse=";"),
                 cor=paste(tcor[[i]][k], collapse=";"), stringsAsFactors=F)
    else
      NULL
  }
  
} # relateToGeneSet 



#' Cell type frequencies in relations to gene sets
#'
#' Count how many times and with which weights cell types were involved
#' in the (L,R,pathway) triples that targeted genes in a gene set.
#'
#' @param rel  The data.frame output by
#' \code{"\link[=BSRInference-class]{relateToGeneSet}"}.
#' @param lr   The data.frame output by
#' \code{"\link{assignCellTypesToInteractions}"}.
#' @param min.n.genes  Minimum number of genes in the gene set for
#' one (L,R,pathway) triple.
#' 
#' @return A list of two slots: t for counting how many times each cell type
#' is involved; s for summing the weights of each involved cell type.
#' @export
#' @examples
#' print('assignCellTypesToInteractions')
#' data(sdc,package='BulkSignalR')
#' bsrdm <- prepareDataset(counts = sdc)
#' bsrdm <- learnParameters(bsrdm, 
#'          null.model = "normal",
#'          quick = FALSE, 
#'          plot.folder = "./",
#'          filename = "sdc",
#'          verbose = TRUE)
#' bsrinf <- initialInference(bsrdm)
#' 
#' # microenvironment cell populations
#' data(immune.signatures, package="BulkSignalR")
#' immune.signatures <- immune.signatures[immune.signatures$signature %in%
#'                                c("B cells","Dentritic cells","Macrophages",
#'                                  "NK cells","T cells","T regulatory cells"),]
#' data("tme.signatures", package="BulkSignalR")
#' signatures <- rbind(immune.signatures,tme.signatures[
#'           tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
#' tme.scores <- scoreSignatures(bsrdm, signatures)
#' 
#' # assignment
#' lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
#' 
#' # relate to p-EMT (should be done in HNSCC normally, not in SDC)
#' data(p.EMT, package='BulkSignalR')
#' p.EMT <- p.EMT$gene
#' triggers <- relateToGeneSet(bsrinf, p.EMT)
#' 
#' # counts
#' cf <- cellTypeFrequency(triggers, lr2ct)
cellTypeFrequency <- function(rel, lr, min.n.genes=1){
  
  rel <- rel[rel$n.genes >= min.n.genes,]
  rel.keys <- paste(rel$L, rel$R, sep="-")
  lr.keys <- paste(lr$L, lr$R, sep="-")
  already <- NULL
  CTs <- NULL
  weights <- NULL
  for (key in unique(rel.keys)){
    rows <- which(lr.keys == key)
    if (length(rows) > 0)
      for (j in rows){
        CTs <- c(CTs, lr$cell.type[j])
        weights <- c(weights, lr$weight[j])
      }
  }
  
  t <- table(CTs)
  s <- NULL
  for (ct in names(t))
    s <- c(s, sum(weights[CTs==ct]))
  names(s) <- names(t)
  list(t=t, s=s)
  
} # cellTypeFrequency
