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
ds <- learnParameters(ds,normal="normal",induced="tumor",verbose=TRUE)
save(ds,file="ds.rda",compress="bzip2")
#load("ds.rda")

# stop cluster, no parallel computations beyond this point [optional]
stopCluster(cl)

# score ligand-receptor interactions
ds.LR <- getCorrelatedLR(ds,min.cor=0.25)
ds.LR <- checkReceptorSignaling(ds,ds.LR)
save(ds.LR,file="ds-LR.rda",compress="bzip2")
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
simpleHeatmap(scores,"SDC-LR-heatmap.pdf",width=6,height=4,pointsize=4)

# correlate with the microenvironment
data(tme.signatures,package="BulkSignalR")
tme.scores <- scoreSignatures(ds,tme.signatures)
dualHeatmap(scores,tme.scores,"SDC-LR-TME-heatmap.pdf",width=6,height=4.5,pointsize=4,vert.p=0.8)

# generate ligand-receptor networks and export in .graphML for Cytoscape or similar tools
gLR <- getLRNetwork(pp,qval.thres=0.01)
write.graph(gLR,file="SDC-LR-network.graphml",format="graphml")

# play around with igraph functions as an alternative to Cytoscape
plot(gLR)
gLR <- set_vertex_attr(gLR,name="size",value=5)
lay <- layout_with_kk(gLR)
plot(gLR,layout=lay)
# community detection
u.gLR <- as.undirected(gLR) # most community finding algorithms work for undirected graphs only
comm <- cluster_edge_betweenness(u.gLR)
plot(comm,u.gLR)
table(comm$membership)
V(u.gLR)$label[comm$membership==28]
# cohesive blocks
cb <- cohesive_blocks(u.gLR)
plot(cb,u.gLR)
h <- hierarchy(cb)
h <- set_vertex_attr(h,name="size",value=7)
plot(h)
bl <- blocks(cb)
bl[sapply(bl,function(x) length(x)<20 && any(V(u.gLR)[x]$label%in%c("CTLA-4","TGFB1","JAG1","CD274")))]


