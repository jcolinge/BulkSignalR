% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.R
\name{getLRNetwork}
\alias{getLRNetwork}
\title{Generate a ligand-receptor network}
\usage{
getLRNetwork(
  bsrinf,
  pval.thres = NULL,
  qval.thres = NULL,
  node.size = 5,
  red.pairs = NULL
)
}
\arguments{
\item{bsrinf}{A BSRInference object.}

\item{pval.thres}{P-value threshold.}

\item{qval.thres}{Q-value threshold.}

\item{node.size}{Default node size in the network.}

\item{red.pairs}{A data frame with columns L (ligands) and R
(receptors) that restrict LR pairs to those listed.}
}
\value{
An \code{igraph} object featuring the ligand-receptor network.
Default colors and node sizes are assigned,
which can be changed afterwards if necessary.
}
\description{
Generate a ligand-receptor network from a ligand-receptor table.
}
\examples{
print('getLRNetwork')
data(sdc,package='BulkSignalR')
bsrdm <- prepareDataset(counts = sdc)
bsrdm <- learnParameters(bsrdm, 
         null.model = "normal",
         quick = FALSE, 
         plot.folder = "./",
         filename = "sdc",
         verbose = TRUE)
bsrinf <- initialInference(bsrdm)
gLR <- getLRNetwork(bsrinf, qval.thres=1e-4)
# plot(gLR)
# write.graph(gLR, file="SDC-LR-network.graphml", format="graphml")

}
