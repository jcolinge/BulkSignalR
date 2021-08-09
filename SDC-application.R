library(data.table)
library(BulkSignalR)
library(igraph)

# activate parallel computing for faster model training [optional]
library(doParallel)
n.proc <- 2
cl <- makeCluster(n.proc)
registerDoParallel(cl)

# prepare data and train the statistical model
data(sdc,package="BulkSignalR")
sample.types <- rep("tumor",ncol(sdc)-2)
ds <- prepareDataset(sdc[,-grep("^N",names(sdc))],sample.types)
ds <- learnParameters(ds,verbose=TRUE,quick=TRUE)
#save(ds,file="ds.rda",compress="bzip2")
#load("ds.rda")

# score ligand-receptor interactions
ds.LR <- getCorrelatedLR(ds,min.cor=0.25)
ds.LR <- checkReceptorSignaling(ds,ds.LR)
#save(ds.LR,file="ds-LR.rda",compress="bzip2")
#load("ds-LR.rda")
pp <- pValuesLR(ds.LR,ds$param)
pp <- pp[order(pp$qval),]
pp[1:20,c("L","R","pval","qval","pw.id","pw.name")] # top 20 interaction with pathway redundancy

# extract gene signatures to report ligand-receptor-downstream pathway scores
pp.red <- reduceToBestPathway(pp)
pp.red[1:20,c("L","R","pval","qval","pw.id","pw.name")] # top 20 interactions without pathway redundancy
sum(pp.red$qval<0.01) # number of significant interactions
p.red.P <- reduceToPathway(pp.red)
signatures <- getLRGeneSignatures(p.red.P,qval.thres=1e-4)
scores <- scoreLRGeneSignatures(ds,signatures,rename.by.pathway=T)
simpleHeatmap(scores,pointsize=8)
simpleHeatmap(scores,"SDC-LR-heatmap.pdf",width=6,height=4,pointsize=4)

# correlate with the microenvironment
data(tme.signatures,package="BulkSignalR")
tme.scores <- scoreSignatures(ds,tme.signatures)
dualHeatmap(scores,tme.scores,pointsize=8,vert.p=0.82)

# generate a ligand-receptor network and export it in .graphML for Cytoscape or similar tools
gLR <- getLRNetwork(pp,qval.thres=0.01)
write.graph(gLR,file="SDC-LR-network.graphml",format="graphml")

# play around with igraph functions as an alternative to Cytoscape
plot(gLR)
plot(gLR,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
# community detection
u.gLR <- as.undirected(gLR) # most community finding algorithms work for undirected graphs only
comm <- cluster_edge_betweenness(u.gLR)
plot(comm,u.gLR,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
table(comm$membership)
V(u.gLR)$label[comm$membership==28]
# cohesive blocks
cb <- cohesive_blocks(u.gLR)
plot(cb,u.gLR,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")
h <- hierarchy(cb)
h <- set_vertex_attr(h,name="size",value=7)
plot(h)
bl <- blocks(cb)
bl[sapply(bl,function(x) length(x)<20 && any(V(u.gLR)[x]$label%in%c("CTLA-4","JAG1","CD274")))]

# generate different ligand-receptor networks depending on the gene signature clusters
mult.net <- getMultipleLRNetworks(ds,pp,n.clusters=4,qval.thres=0.01)
plot(mult.net$hclust.spl)
table(mult.net$clusters)
simpleHeatmap(mult.net$scores,file.name="SDC-clusters-LR-heatmap.pdf",dend.spl=as.dendrogram(mult.net$hclust.spl),n.col.clust=4,row.names=FALSE) # cluster numbers are different in the ComplexHeatmap output
lay.1 <- layout_with_kk(mult.net$networks[[1]])
plot(mult.net$networks[[1]],
     layout=lay.1,
     edge.arrow.size=0.3,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
ug.1 <- as.undirected(mult.net$networks[[1]])
cb.1 <- cohesive_blocks(ug.1)
plot(cb.1,ug.1,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")
lay.2 <- layout_with_kk(mult.net$networks[[2]])
plot(mult.net$networks[[2]],
     layout=lay.2,
     edge.arrow.size=0.3,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
ug.2 <- as.undirected(mult.net$networks[[2]])
cb.2 <- cohesive_blocks(ug.2)
plot(cb.2,ug.2,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")
lay.3 <- layout_with_kk(mult.net$networks[[3]])
plot(mult.net$networks[[3]],
     layout=lay.3,
     edge.arrow.size=0.3,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
ug.3 <- as.undirected(mult.net$networks[[3]])
cb.3 <- cohesive_blocks(ug.3)
plot(cb.3,ug.3,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")
lay.4 <- layout_with_kk(mult.net$networks[[4]])
plot(mult.net$networks[[4]],
     layout=lay.4,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
ug.4 <- as.undirected(mult.net$networks[[4]])
cb.4 <- cohesive_blocks(ug.4)
plot(cb.4,ug.4,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")

# generate a ligand-receptor network complemented with intracellular, receptor downstream pathways [computations are a bit longer here]
gLRintra <- getLRIntracellNetwork(pp,qval.thres=0.01)
write.graph(gLRintra,file="SDC-LR-intracellular-network.graphml",format="graphml")
lay <- layout_with_kk(gLRintra)
plot(gLRintra,
     layout=lay,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# reduce complexity by focusing on strongly targeted pathways
top <- unique(pp[pp$pval<1e-10,c("pw.id","pw.name")])
top
gLRintra.res <- getLRIntracellNetwork(pp,qval.thres=0.01,restrict.pw=top$pw.id)
lay <- layout_with_fr(gLRintra.res)
plot(gLRintra.res,
     layout=lay,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# generate different ligand-receptor-intracellular downstream pathway networks depending on the gene signature clusters [lengthy computation]
mult.net.intra <- getMultipleLRIntracellNetworks(ds,pp,n.clusters=4,qval.thres=0.01)
plot(mult.net$hclust.spl)
table(mult.net$clusters)
simpleHeatmap(mult.net$scores,dend.spl=as.dendrogram(mult.net$hclust.spl),n.col.clust=4,row.names=FALSE) # cluster numbers are different in the ComplexHeatmap output
lay.4 <- layout_with_fr(mult.net.intra$networks[[4]])
plot(mult.net.intra$networks[[4]],
     layout=lay.4,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# stop cluster if parallel computation was used [optional]
stopCluster(cl)

