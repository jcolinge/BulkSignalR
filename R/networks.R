#' Generate a ligand-receptor network
#'
#' Generate a ligand-receptor network from a ligand-receptor table.
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @param LLR.thres     Log-likelihood threshold.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param node.size     Default node size in the network.
#' @return An \code{igraph} object featuring the ligand-receptor network. Default colors and node sizes are assigned,
#' which can be changed afterwards if necessary.
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
#' gLR <- getLRNetwork(pp,qval.thres=0.01)
#' plot(gLR)
#' write.graph(gLR,file="SDC-LR-network.graphml",format="graphml")
#' }
getLRNetwork <- function(pairs,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL,node.size=5){
    
  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")
  red.mode.R <- (sum(regexpr("^{",pairs$R,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$R,perl=TRUE)==nchar(pairs$R))==nrow(pairs))
  red.mode.L <- (sum(regexpr("^{",pairs$L,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$L,perl=TRUE)==nchar(pairs$L))==nrow(pairs))
  if (red.mode.L || red.mode.R)
    stop("Does not work with ligand-receptor pairs already reduced to the ligands or the receptors")
  
  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>=LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<=pval.thres,]
    else
      pairs <- pairs[pairs$qval<=qval.thres,]
    
  pairs <- reduceToBestPathway(pairs)
  ref <- paste(pairs$L,pairs$R,sep="||")
  dup <- duplicated(ref)
  pairs <- pairs[!dup,] # different pathways for the same LR pair can have the same P-value
  g <- igraph::graph_from_data_frame(pairs[,c("L","R")],directed=TRUE) # different pathways for the same LR pair can have the same P-value

  g.names <- igraph::vertex_attr(g,"name")
  g <- igraph::set_vertex_attr(g,name="size",value=node.size)
  g <- igraph::set_vertex_attr(g,name="label",value=g.names)
  g.types <- stats::setNames(rep("ligand",length(g.names)),g.names)
  g.types[g.names%in%pairs$R] <- "receptor"
  g <- igraph::set_vertex_attr(g,name="node.type",value=g.types[g.names])
  g.colors <- rep("green",length(g.names))
  g.colors[g.names%in%pairs$R] <- "red"
  g <- igraph::set_vertex_attr(g,name="color",value=g.colors)
  g.shapes <- rep("circle",length(g.names))
  g.shapes[g.names%in%pairs$R] <- "square"
  g <- igraph::set_vertex_attr(g,name="shape",value=g.shapes)
  g <- igraph::set_edge_attr(g,name="edge.type",value="LR")
  el <- as_edgelist(g)
  pval <- NULL
  qval <- NULL
  corr <- NULL
  for (i in 1:nrow(el)){
    j <- which(pairs$L==el[i,1] & pairs$R==el[i,2])
    pval <- c(pval,pairs$pval[j])
    qval <- c(qval,pairs$qval[j])
    corr <- c(corr,pairs$LR.corr[j])
  }
  g <- igraph::set_edge_attr(g,name="corr",value=corr)
  g <- igraph::set_edge_attr(g,name="pval",value=pval)
  g <- igraph::set_edge_attr(g,name="log10.pval",value=-log10(pval))
  g <- igraph::set_edge_attr(g,name="qval",value=qval)
  g <- igraph::set_edge_attr(g,name="log10.qval",value=-log10(qval))
  
} # getLRNetwork


#' Generate multiple ligand-receptor networks
#'
#' Generate a ligand-receptor network for each cluster of samples based on ligand-receptor gene signatures.
#'
#' @param ds         A BulkSignalR data set.
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @param n.clusters    The number of clusters.
#' @param min.score     The minimum required average z-score of a gene signature in a cluster.
#' @param LLR.thres     Log-likelihood threshold in \code{pairs}.
#' @param pval.thres    P-value threshold in \code{pairs}.
#' @param qval.thres    Q-value threshold in \code{pairs}.
#' @param cut.p          Proportion of top and bottom values for thresholding.
#' @return A list containing the main elements of the clustering analysis as well as an inner list
#' containing the \code{igraph} objects created for each cluster.
#' 
#' The elements of clustering are the retained global ligand-receptor table, the correspondinf gene signatures and their
#' scores across \code{ds} samples, the resulting hierarchical cultering and sample membership. 
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
#' gLR <- getLRNetwork(pp,qval.thres=0.01)
#' plot(gLR)
#' write.graph(gLR,file="SDC-LR-network.graphml",format="graphml")
#' }
getMultipleLRNetworks <- function(ds,pairs,n.clusters,min.score=0,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL,cut.p=0.01){
  
  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")
  red.mode.R <- (sum(regexpr("^{",pairs$R,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$R,perl=TRUE)==nchar(pairs$R))==nrow(pairs))
  red.mode.L <- (sum(regexpr("^{",pairs$L,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$L,perl=TRUE)==nchar(pairs$L))==nrow(pairs))
  if (red.mode.L || red.mode.R)
    stop("Does not work with ligand-receptor pais already reduced to the ligands or the receptors")
  if (n.clusters<1 || n.clusters>ncol(ds$ncounts))
    stop("n.clusters must be > 1 and <= ncol(ds$ncounts)")

  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>=LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<=pval.thres,]
  else
    pairs <- pairs[pairs$qval<=qval.thres,]
  
  pairs <- reduceToBestPathway(pairs)
  ref <- paste(pairs$L,pairs$R,sep="||")
  dup <- duplicated(ref)
  pairs <- pairs[!dup,] # different pathways for the same LR pair can have the same P-value
  signatures <- getLRGeneSignatures(pairs,pval.thres=1) # already filtered
  scores <- scoreLRGeneSignatures(ds,signatures)
  mat.c <- .cutExtremeValues(scores,cut.p)
  di.spl <- stats::dist(t(mat.c))
  hc.spl <- stats::hclust(di.spl,method="ward.D")
  cluster <- cutree(hc.spl,n.clusters)
  networks <- foreach::foreach(i=1:n.clusters,.combine=c) %do% {
    LRscores <- rowMeans(scores[,cluster==i])
    sel.pairs <- pairs[LRscores>=min.score,]
    list(getLRNetwork(sel.pairs,pval.thres=1))
  }
  
  list(pairs=pairs,hclust.spl=hc.spl,signatures=signatures,scores=scores,n.clusters=n.clusters,clusters=cluster,networks=networks)

} # getMultipleLRNetworks


.edgesLRIntracell <- function(pairs,pw,id.col,gene.col,min.cor=0.3,signed=FALSE){
  
  directed.int <- c("controls-state-change-of","catalysis-precedes","controls-expression-of","controls-transport-of","controls-phosphorylation-of")

  arcs <- NULL
  for (i in 1:nrow(pairs)){
    r <- pairs$R[i]
    p <- pairs$pw.id[i]
    tg <- strsplit(pairs$target.genes[i],split=";")[[1]]
    corr <- as.numeric(strsplit(pairs$all.corr[i],split=";")[[1]])
    if (signed)
      targets <- tg[corr>=min.cor]
    else
      targets <- tg[abs(corr)>=min.cor]
    
    # build a hybrid directed/undirected graph
    a.iter <- data.frame(from=pairs$L[i],to=r,edge.type="LR",stringsAsFactors=FALSE)
    genes.in.pw <- pw[pw[[id.col]]==p,gene.col]
    int <- SingleCellSignalR::PwC_ReactomeKEGG[SingleCellSignalR::PwC_ReactomeKEGG$a.gn%in%genes.in.pw & SingleCellSignalR::PwC_ReactomeKEGG$b.gn%in%genes.in.pw,]
    directed <- int$type%in%directed.int
    d.int <- unique(rbind(int[directed,c("a.gn","b.gn")],int[!directed,c("a.gn","b.gn")],int[!directed,c("b.gn","a.gn")]))
    g <- igraph::graph_from_data_frame(d.int,directed=TRUE)
    targets <- intersect(targets,c(d.int$a.gn,d.int$b.gn))

    # keep shortest paths from the receptor to the targets only
    if ((r%in%d.int$a.gn || r%in%d.int$b.gn) && length(targets)>0){
      paths <- suppressWarnings(igraph::shortest_paths(g,from=V(g)[V(g)$name==r],to=V(g)[V(g)$name%in%targets],output="vpath"))
      if (length(paths$vpath)>0)
        for (j in 1:length(paths$vpath)){
          vertices <- paths$vpath[[j]]
          if (length(vertices)>0)
            for (k in 2:length(vertices)){
              from <- V(g)$name[vertices[k-1]]
              to <- V(g)$name[vertices[k]]
              edge.type <- SingleCellSignalR::PwC_ReactomeKEGG[SingleCellSignalR::PwC_ReactomeKEGG$a.gn==from & SingleCellSignalR::PwC_ReactomeKEGG$b.gn==to,"type"][1]
              a.iter <- rbind(a.iter,data.frame(from=from,to=to,edge.type=edge.type,stringsAsFactors=FALSE))
            }
        }
    }
    arcs <- rbind(arcs,unique(a.iter))
  }
  
  unique(arcs)

} # .edgesLRIntracell


getLRIntracellNetwork <- function(pairs,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL,min.cor=0.3,signed=FALSE,restrict.pw=NULL,node.size=5){
  
  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")
  red.mode.R <- (sum(regexpr("^{",pairs$R,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$R,perl=TRUE)==nchar(pairs$R))==nrow(pairs))
  red.mode.L <- (sum(regexpr("^{",pairs$L,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$L,perl=TRUE)==nchar(pairs$L))==nrow(pairs))
  if (red.mode.L || red.mode.R)
    stop("Does not work with ligand-receptor pairs already reduced to the ligands or the receptors")
  
  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>=LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<=pval.thres,]
  else
    pairs <- pairs[pairs$qval<=qval.thres,]

  pool <- unique(c(pairs$L,pairs$R))
  all.edges <- NULL

  # reactome pathways
  pairs.react <- pairs[grep("^R-",pairs$pw.id),]
  if (nrow(pairs.react)>0){
    ids <- unique(reactome[reactome$`Gene name`%in%pool,"Reactome ID"])
    react <- reactome[reactome$`Reactome ID`%in%ids,]
    if (!is.null(restrict.pw))
      react <- react[react$`Reactome ID`%in%restrict.pw,]
    all.edges <- .edgesLRIntracell(pairs.react,react,"Reactome ID","Gene name",min.cor)
  }
  
  # GOBP
  pairs.go <- pairs[grep("^GO:",pairs$pw.id),]
  if (nrow(pairs.go)>0){
    ids <- unique(gobp[gobp$`Gene name`%in%pool,"GO ID"])
    go <- gobp[gobp$`GO ID`%in%ids,]
    if (!is.null(restrict.pw))
      go <- go[go$`GO ID`%in%restrict.pw,]
    all.edges <- unique(rbind(all.edges,.edgesLRIntracell(pairs.go,go,"GO ID","Gene name",min.cor)))
  }
  
  # generate igraph object -------------------
  
  # some LR interactions are also given as in-complex-with interactions, get rid of the latter
  ref <- paste(all.edges[[1]],all.edges[[2]],sep="||")
  dup <- ref[duplicated(ref)]
  all.edges <- all.edges[!(ref%in%dup) | (ref%in%dup & all.edges$edge.type=="LR"),] 
  ref <- paste(all.edges[[1]],all.edges[[2]],sep="||")
  duppl <- duplicated(ref)
  if (sum(duppl)>0)
    all.edges <- all.edges[!duppl,] # some duplicates remain, eliminate

  # build graph
  directed.int <- c("controls-state-change-of","catalysis-precedes","controls-expression-of","controls-transport-of","controls-phosphorylation-of","LR")
  directed <- all.edges$edge.type%in%directed.int
  ret <- all.edges[!directed,]
  from <- ret$from
  ret$from <- ret$to
  ret$to <- from
  d.int <- unique(rbind(all.edges,ret))
  g <- igraph::graph_from_data_frame(d.int,directed=TRUE)

  g.names <- igraph::vertex_attr(g,"name")
  g <- igraph::set_vertex_attr(g,name="size",value=node.size)
  g <- igraph::set_vertex_attr(g,name="label",value=g.names)
  g.types <- stats::setNames(rep("downstream.gene",length(g.names)),g.names)
  g.types[g.names%in%pairs$R] <- "receptor"
  g.types[g.names%in%pairs$L] <- "ligand"
  g <- igraph::set_vertex_attr(g,name="node.type",value=g.types[g.names])
  g.colors <- rep("gray80",length(g.names))
  g.colors[g.names%in%pairs$R] <- "red"
  g.colors[g.names%in%pairs$L] <- "green"
  g <- igraph::set_vertex_attr(g,name="color",value=g.colors)
  g.shapes <- rep("rectangle",length(g.names))
  g.shapes[g.names%in%pairs$R] <- "square"
  g.shapes[g.names%in%pairs$L] <- "circle"
  g <- igraph::set_vertex_attr(g,name="shape",value=g.shapes)

} # getLRIntracellNetwork


getMultipleLRIntracellNetworks <- function(){
  
  
} # getMultipleLRIntracellNetworks