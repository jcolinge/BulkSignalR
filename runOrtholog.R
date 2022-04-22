library(glue)
library(data.table)
library(devtools)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(paxtoolsr)
library(doParallel)
library(stringr)
library(viridis)
library(ggridges)
library(ggrepel)
library(biomaRt)
library(rjson)
library(stringi)

switch <- '/data/villemin/'
#devtools::load_all(glue("{switch}/code/refactoringPackagesR/JC/BulkSignalR/"),TRUE)
#remove.packages(pkgs="BulkSignalR",lib="/data/USERS/villemin/anaconda3/envs/r4.1.3/lib/R/library")
devtools::install(glue("{switch}code/refactoringPackagesR/JC/BulkSignalR"))
library(BulkSignalR)

set.seed(123)

# activate parallel computing for faster model training [optional]
n.proc <- 2
cl     <- makeCluster(n.proc)
registerDoParallel(cl)
#registerDoMC(cores=n.proc)
#registerDoSEQ()

##############################################################################################
############################    Load Config Parameters                   #####################
##############################################################################################

file <- c("GSE133079.Mouse_Bulk_norm_counts.txt.annoted.csv") # Systematic Identification of Cell-Cell Communication Networks in the Developing Brain. iScience 2019 Nov 22; (Mouse Brain)
#file <- c("DnReioTPM.csv") # A high-resolution mRNA expression time course of embryonic development in zebrafish.eLife 2017

config        <- fromJSON(file=glue("{switch}code/refactoringPackagesR/JC/BulkSignalR/data/config.json"))
config.string <- stri_join_list(config[file], sep = "_", collapse = NULL)

threshold  <- as.numeric(config[[file]]$threshold)
min.cor    <- as.numeric( config[[file]]$min.cor)
min.count  <- as.numeric( config[[file]]$min.count)
normalize  <- config[[file]]$normalize
ortholog  <- config[[file]]$ortholog

##############################################################################################
############################    Functions                   ##################################
##############################################################################################

findOrthoGenesCustom <- function(from_organism ="mmusculus",from_values=c("TP53"),BIOMART_CACHE="") {
          print("findOrthoGenesCustom")

          #Remove error due fo file system => Error: “No Lock available”
          if (BIOMART_CACHE!= "")
               Sys.setenv(BIOMART_CACHE=BIOMART_CACHE) 
          #Remove recurrent error => Error: “SSL certificate problem”
          httr::set_config(httr::config(ssl_verifypeer = FALSE))

          from_organism_gene_ensembl <- paste0(from_organism, c("_gene_ensembl"))

          to_organism <- "hsapiens"
          to_organism_features <- paste0(to_organism, c("_homolog_ensembl_gene","_homolog_associated_gene_name"))
            
          ensembl <- useEnsembl(biomart = "ensembl",verbose=TRUE, mirror = "useast")

          # list the available datasets in this Mart
          ensembl_datasets <- listDatasets(mart = ensembl)

          # Check organism definition is ok
          if (from_organism_gene_ensembl %in% ensembl_datasets$dataset ) {
          
               human <- useEnsembl("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
               organism_from <- useEnsembl("ENSEMBL_MART_ENSEMBL", dataset = from_organism_gene_ensembl)

               orthologs_dictionnary <- getLDS(
                      attributes = c("external_gene_name"),
                      filters = "external_gene_name", values = from_values,
                      mart = organism_from,
                      attributesL = c("hgnc_symbol"), 
                      martL = human)

               # It's not possible to keep the same human symbol several times.
               # Then you are not able to make a dictionnary (key->value) and convert back to initial symbol in another species.
               # That's the pitfall of gene symbols. 
               orthologs_dictionnary <- orthologs_dictionnary[!duplicated(orthologs_dictionnary[, c("HGNC.symbol")]),]
               # Also you can  have several human genes corresponding to one ortholog in another species.
               # In order to avoid spurious result in converting back from the hsapiens organism ,
               # we keep also only one unique gene for this other considered species.
               # You will miss stuffs but results will be cleaner and easier to interpret.
               orthologs_dictionnary <- orthologs_dictionnary[!duplicated(orthologs_dictionnary[, c("Gene.name")]),] #Lefty2

               rownames(orthologs_dictionnary) <- orthologs_dictionnary$HGNC.symbol
          } 

          else { stop("Organism not found !") }
    
    print(head(orthologs_dictionnary,10))

    orthologs_dictionnary


} #findOrthoGenesCustom 


##############################################################################################
############################                       ###########################################
##############################################################################################
     bench.dir         <- glue("{switch}/data/public/LR/bench/")
     bsrdm.file        <- glue("{bench.dir}/{file}_{config.string}_bsrdm.rds")
     bsrinf.file       <- glue("{bench.dir}/{file}_{config.string}_bsrinf.rds")

     print(glue("{bench.dir}/{file}"))

     matrix.expression.not.human <- fread(glue("{bench.dir}/{file}"),data.table=FALSE)
     matrix.expression.not.human <- distinct(matrix.expression.not.human, V1, .keep_all = TRUE)

     # mandatory
     rownames(matrix.expression.not.human) <- matrix.expression.not.human[,1]
     matrix.expression.not.human           <- matrix.expression.not.human[,-1]
     matrix.expression.not.human <- matrix.expression.not.human[rowSums(matrix.expression.not.human)>0, ]     # optional

     print ("==> Find Orthologous Genes & Convert to Human: ") 

     # Use internal package function to recover 
     # a dictionnary of associated genes between human and a distinct organism

     ortholog.dict    <- findOrthoGenes (from_organism = ortholog,from_values = rownames(matrix.expression.not.human))
     

     # As an alternative you can submit your own dictionnary of orthologous
     # ortholog.dict    <- findOrthoGenesCustom(
     #    from_organism = ortholog,
     #     from_values = rownames(matrix.expression.not.human),
     #     BIOMART_CACHE = "/data2/USERS/shared/cache")
     #)
     write.table(ortholog.dict, glue("{bench.dir}/{file}_{config.string}_dict.orthogene.tsv"), col.names=NA , sep="\t", quote=FALSE)

     # Use internal package function to convert any distinct organism from hsapiens to hsapiens
     # based on dictionnary of associated genes between human and a distinct organism

     matrix.expression.human <- convertToHuman(counts = matrix.expression.not.human,dictionnary = ortholog.dict)
 
     print ("==> PrepareDataset & LearnParameters: ") 

     bsrdm <- prepareDataset(
          counts = matrix.expression.human,
          min.count = min.count,
          normalize = normalize,
          species = ortholog, 
          conversion.dict = ortholog.dict)

     str(bsrdm)
     show(bsrdm)
     bsrdm <- learnParameters(bsrdm, quick = TRUE, verbose = TRUE)

     saveRDS(bsrdm, bsrdm.file)
     bsrdm <- readRDS(bsrdm.file)
     str(bsrdm)

     print ("==> InitialInference : Scoring of ligand-receptor interactions : ") # folder is created in current dir

     bsrinf <- initialInference(bsrdm)
     saveRDS(bsrinf, bsrinf.file)
     bsrinf <- readRDS(bsrinf.file)

     str(bsrinf)
     # Convert back bsrinf to original species
     # Will then be used with normal workflow.
     bsrinf <- resetToInitialOrganism(bsrinf,conversion.dict=ortholog.dict)
     #str(bsrinf)
     # Corrected : Necessary if you don't want a bug in scoreLRGeneSignatures
     #bsrdm@ncounts<- as.matrix(matrix.expression.not.human) 
     str(bsrinf)

     # Reducing the LR analysis --------------------------     
     write.table(as.data.frame(LRinter(bsrinf)), glue("{bench.dir}/{file}_{config.string}_LR.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

     print ("==> ReduceToLigand : ") 
     bsrinf.redL <- reduceToLigand(bsrinf)
     write.table(as.data.frame(LRinter(bsrinf.redL)), glue("{bench.dir}/{file}_{config.string}_LR.reduceToLigand.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

     print ("==> ReduceToReceptor : ") 
     bsrinf.redR <- reduceToReceptor(bsrinf)
     write.table(as.data.frame(LRinter(bsrinf.redR)), glue("{bench.dir}/{file}_{config.string}_LR.reduceToReceptor.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
     
     print ("==> ReduceToPathway : ") 
     bsrinf.redP <- reduceToPathway(bsrinf) # Reduced to only report one row per pathway.
     # For a given pathway, the reported P-values and target genes are those of the best ligand-receptor pair that was in this pathway.
     # All the ligands and all the receptors forming pairs related to a certain pathway, are combined.
     write.table(as.data.frame(LRinter(bsrinf.redP)), glue("{bench.dir}/{file}_{config.string}_LR.reduceToPathway.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

     print ("==> ReduceToBestPathway : ") # Keep  best ligand-receptor pair.
     bsrinf.red    <- reduceToBestPathway(bsrinf) #  # Patwhay can be on several rows.
     # Only one Ligand & Receptor per line.
     write.table(as.data.frame(LRinter(bsrinf.red)), glue("{bench.dir}/{file}_{config.string}_LR.reduceToBestPathway.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

     print ("==> ReduceToBestPathway after reduceToPathway : ")  # To avoid redundancy.
     bsrinf.redPP <- reduceToBestPathway(bsrinf.redP) # Keep one pathway per ligand-receptor pair.
     # All the ligands and all the receptors forming pairs related to a certain pathway, are combined.
     write.table(as.data.frame(LRinter(bsrinf.redPP)), glue("{bench.dir}/{file}_{config.string}_LR.reduceToPathway.reduceToBestPathway.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
     
     # Gene signatures --------------------------------------------

     # extract gene signatures to report combined ligand-receptor and
     # receptor downstream pathway scores
     sum(LRinter(bsrinf.redPP)$qval < threshold) # number of significant interactions
     sum(LRinter(bsrinf.redP)$qval < threshold)

     print ("==> getLRGeneSignatures : ") 
     bsrsig.red <- getLRGeneSignatures(bsrinf.red, qval.thres=threshold)

     print ("==> scoreLRGeneSignatures : scoresLR ") 
     scoresLR <- scoreLRGeneSignatures(bsrdm,bsrsig.red,name.by.pathway=FALSE)
     write.table(scoresLR, glue("{bench.dir}/{file}_{config.string}_scoreLR.tsv"),  col.names=NA,sep="\t",   quote=FALSE)

     print ("==> getLRGeneSignatures : ") 
     bsrsig.redPP <- getLRGeneSignatures(bsrinf.redPP,qval.thres=threshold)

     print ("==> scoreLRGeneSignatures : scorepathways  ") 
     scoresPathway <- scoreLRGeneSignatures(bsrdm,bsrsig.redPP,name.by.pathway=TRUE)
     write.table(scoresPathway, glue("{bench.dir}/{file}_{config.string}_scorePathways.tsv"),  col.names=NA,sep="\t",   quote=FALSE)
   
stopCluster(cl)

print("Done")

stop()
