getLRNetwork <- function(pairs,LLR.thres=NULL,pval.thres=NULL,qval.thres=NULL){
    
  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P/Q-values and LLR present in LR table")
  if (t == 0)
    stop("P/Q-values or LLR must be available in LR table")
  if (is.null(LLR.thres)+is.null(pval.thres)+is.null(qval.thres) != 2)
    stop("One selection criterion out of P-value, Q-value, or LLR only")
  
  if (!is.null(LLR.thres))
    pairs <- pairs[pairs$LLR>LLR.thres,]
  else
    if (!is.null(pval.thres))
      pairs <- pairs[pairs$pval<pval.thres,]
    else
      pairs <- pairs[pairs$qval<qval.thres,]
    
  pairs <- unique(pairs[,c("L","R")])
  g <- igraph::graph_from_data_frame(pairs,directed=TRUE)
  
  g.names <- vertex_attr(gLR,"name")
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

} # getLRNetwork
