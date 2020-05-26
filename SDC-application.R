library(data.table)
library(BulkSignalR)

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

# stop cluster, no paralle computations beyond this point [optional]
stopCluster(cl)

# score ligand-receptor interactions
ds.LR <- getCorrelatedLR(ds,min.cor=0.25)
ds.LR <- checkReceptorSignaling(ds,ds.LR)
pp <- pValuesLR(ds.LR,ds$param)

# extract gene signatures to report ligand-receptor-downstream pathway scores
p.red.P <- reduceToPathway(pp)
signatures <- getLRGeneSignatures(p.red.P,qval.thres=1e-5)
scores <- scoreLRGeneSignatures(ds,signatures,rename.by.pathway=T)
simpleHeatmap(scores,"SDC-LR-heatmap.pdf",width=6,height=4.5,pointsize=4)

# correlate with the microenvironment
data(tme.signatures,package="BulkSignalR")
tme.scores <- scoreSignatures(ds,tme.signatures)
dualHeatmap(scores,tme.scores,"SDC-LR-TME-heatmap.pdf",width=6,height=5,pointsize=4,vert.p=0.85)
