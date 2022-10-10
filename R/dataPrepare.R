#' Modify PwC_ReactomeKEGG database
#'
#' User can define a dataframe with 3 columns named :
#' respectively a.gn, b.gn and  type
#' This can be used to replace the existing
#' PwC_ReactomeKEGG database.
#'
#' @param db     A dataframe with 3 column names :
#' a.gn, b.gn and type.
#'
#' @return NULL
#'
#' @export
#' @examples
#' print('resetInteractions')
#' 
resetInteractions <- function(db=data.frame(a.gn="A2M",b.gn="LRP1",type="controls-expression-of")) {

    if (all(colnames(db)==  c('a.gn','b.gn','type'))){
            # envir = package.BulkSignalR
            assign("PwC_ReactomeKEGG", db, envir = .GlobalEnv)
    } else {
      stop(paste0("db should be a dataframe with ",
            "3 columns named : 'a.gn' and 'b.gn'",
            " and 'type'."))
    }
} # resetInteractions

#' Modify LRdb database
#'
#' User can define a dataframe with 2 columns named
#' respectively ligand and receptor.
#' This can be used to extand or replace the existing
#' LRdb.
#'
#' @param db     A dataframe with 2 column names :
#' ligand and receptor.
#' @param switch  By default set to FALSE, it extends the
#' existing LRdb database. If TRUE, LRdb is replaced by user 
# provided database.
#'
#' @return NULL
#'
#' @export
#' @examples
#' print('resetLRdb')
#' data(sdc,package='BulkSignalR')
#'resetLRdb(db=data.frame(ligand="A2M",receptor="LRP1"),switch=FALSE)
resetLRdb <- function(db=data.frame(ligand="A2M",receptor="LRP1"),switch=FALSE ) {

    if(colnames(db)[1]=='ligand' &  colnames(db)[2]=='receptor'){
      
        if(switch){
            assign("LRdb", unique(db[,c('ligand','receptor')]), envir = .GlobalEnv)
        }
        else {  
            updated.LRdb <- rbind(LRdb[,c('ligand','receptor')],db[,c('ligand','receptor')])
            assign("LRdb", unique(updated.LRdb), envir = .GlobalEnv)
        }
    } else {
      stop(paste0("db should be a dataframe with ",
            "2 columns named : 'ligand' and 'receptor'."))
    }

}

#' Transform gmt file to dataframe 
#' in order to update core databases
#' of BulkSignalR.
#'
#' Pathways are defined in Reactome and
#' GoBP databases.
#' This can be replaced by more recent versions of
#' these databases using gmt files.
#'
#' \code{resetDownstreamPathways} is a function
#' we provide to user to refresh REACTOME 
#' and GO-BP content included in BulkSignalR.
#'
#' Then dataframe produced here can be loaded by
#' \code{resetDownstreamPathways}.
#'
#' We note discrepancy between format avaibable 
#' over internet. For example, GSEA-MSigDB.org
#' does not provide GO ID in its gmt format.
#'
#' Function is designed to work with these ressources
#' - For GO-BP. (Bader Lab provides updated versions
# of gmt files every month)
#' http://download.baderlab.org/EM_Genesets/
#' - For Reactome. (Directly from their website)
#' https://reactome.org/download/current/
#'
#' The code is inspired from read.gmt function
#' from gsa R package.
#'
#' @param filename    Path to GMT file
#' @param nameDB      Two options "GO-BP" or "REACTOME"
#'
#' @return gmt file as dataframe
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#'
#' @export
#' @examples
#' print('gmtToDataframe')
#' if(FALSE)
#'    gmtToDataframe(filename,"GO-BP")
#'
gmtToDataframe<- function(filename,nameDB=c("GO-BP","REACTOME")){

        if (! nameDB %in% c("GO-BP","REACTOME")){
           stop("GO-BP and REACTOME are the only keywords alllowed.")
        }

        if (! file.exists(filename)){
           stop("File doesn't exist.")
        }

        a<-scan(filename,what=list("",""),sep="\t", 
            quote=NULL, fill=T, flush=T,multi.line=F)

        if(nameDB=="GO-BP"){
            a <- setNames(a, c("ID","Description"))
        }

        if(nameDB=="REACTOME"){
            # Reactome is reversed compared to GO
            a <- setNames(a, c("Description","ID"))
            # Reorder
            a <- a[c("ID","Description")]
        } 
       
        geneset.ids<-a[1][[1]]
        geneset.descriptions<-a[2][[1]]
      
        dd<-scan(filename,what="",sep="\t", quote=NULL)
        dd<- dd[dd!=""]

        nn<-length(geneset.ids)
        n<-length(dd)
        ox<-rep(NA,nn)

        # Clean the two scans 
        # Only for GO where you need to split to get the ID
        if(nameDB=="GO-BP"){
            for(u in 1:nn){
                if (grepl("%", geneset.ids[u])){
                    geneset.ids[u]<-strsplit(geneset.ids[u] , "%")[[1]][3]
                }
            }    
            for(v in 1:n){
                    if( grepl("%", dd[v])){
                      dd[v]<-strsplit(dd[v] , "%")[[1]][3]
                    }
            } 
        }  
        # Compute indice of pathway
        ii<-1
        for(i in 1:nn){
            if(nameDB=="GO-BP"){
             while((dd[ii]!=geneset.ids[i]) | (dd[ii+1]!=geneset.descriptions[i]) ){
               ii=ii+1
             }
            }
            if(nameDB=="REACTOME"){
             while((dd[ii]!=geneset.descriptions[i]) | (dd[ii+1]!=geneset.ids[i]) ){
               ii=ii+1
             }
            }
         ox[i]=ii   
         ii=ii+1
        }
   
        genesets=vector("list",nn)

        for(i in 1:(nn-1)){
          i1<-ox[i]+2
          i2<-ox[i+1]-1
          geneset.descriptions[i]<-dd[ox[i]+1]
          geneset.ids[i]<-dd[ox[i]] 
          genesets[[i]]<-dd[i1:i2]
        }

        geneset.ids[nn]<-dd[ox[nn]] 
        geneset.descriptions[nn]=dd[ox[nn]+1]
        genesets[[nn]]=dd[(ox[nn]+2):n]

        data=list(geneset.ids=geneset.ids,
                  geneset.descriptions=geneset.descriptions,
                  genesets=genesets
                 )
          
        dataframeFromGmt <-foreach(i= 1:length(data$geneset.ids), 
                    .combine = rbind) %dopar% { 
   
           data.frame(a=rep(data$geneset.ids[[i]],length(data$genesets[[i]])),
                      b=data$genesets[[i]][1:length(data$genesets[[i]])],
                      c=rep(data$geneset.descriptions[[i]],length(data$genesets[[i]])))
        }

        # Due to the fact React and Go are organized differently
        if(nameDB=="REACTOME"){
            names(dataframeFromGmt)<-c('Reactome name','Gene name','Reactome ID')
            dataframeFromGmt <- dataframeFromGmt[c('Reactome ID','Gene name','Reactome name')]
        }

        if(nameDB=="GO-BP")
            names(dataframeFromGmt)<-c('GO ID','Gene name','GO name')

return(dataframeFromGmt)

} #gmtToDataframe

#' Reset database
#'
#' Pathways are defined in Reactome and
#' Gobp databases.
#' This can be replaced by more recent versions.
#' User can parse a gmt file and export
#' a dataframe with necessary information. 
#' See \code{gmtToDataframe}.
#' 
#' @param db     A dataframe with 3 columns names :
#' "Reactome ID" , "Gene name","Reactome name".
#' "GO ID", "Gene name",  "GO name"
#' @param nameDB    "GO-BP" or "REACTOME"
#' @return NULL
#'
#' @export
#' @examples
#' print('resetDownstreamPathways')
#' 
resetDownstreamPathways <- function(
    db=data.frame("Reactome ID"="R-HSA-109582" ,
     "Gene name"="IGKV2-28",
     "Reactome name"="Hemostasis"),
    nameDB=c("GO-BP","REACTOME")) {

    if (! nameDB %in% c("GO-BP","REACTOME")){
        stop("GO-BP and REACTOME are the only keywords alllowed.")
    }

    if(nameDB=="REACTOME"){
        if (all(colnames(db)==  c('Reactome ID','Gene name','Reactome name')))
             reactome <- db 
            
       else  {   
          stop(paste0("db should be a dataframe with "
              ,"3 columns named 'Reactome ID', ",
               ,"'Gene name' and 'Reactome name'."))
        }
    }

    if(nameDB=="GO-BP"){
        if (all(colnames(db)== c('GO ID','Gene name','GO name')))
            gobp <- db 
        else {
          stop(paste0("db should be a dataframe with "
                ,"3 columns named 'GO ID' "
                ,"'Gene name' and 'GO name'."))
        }
        
    }    
 

} # resetDownstreamPathways

#' Prepare a BSRDataModel object from expression data
#'
#' Take a matrix or data frame containing RNA sequencing,
#' microarray, or expression proteomics data and returns a BSRDataModel
#' object ready for subsequent training. Normally, BSRDataModel objects
#' are not instantiated directly, but through this function.
#'
#' @param counts     A table or matrix of read counts.
#' @param species    Data were obtained for this organism.
#' @param normalize  A logical indicating whether \code{counts} should be
#'   normalized according to \code{method} or if it was normalized beforehand.
#' @param symbol.col The index of the column containing the gene symbols in case
#'   those are not the row names of \code{counts} already.
#' @param min.count  The minimum read count of a gene to be considered expressed
#'   in a sample.
#' @param prop       The minimum proportion of samples where a gene must be
#'   expressed higher than \code{min.count} to keep that gene.
#' @param method     The normalization method ('UQ' for upper quartile or 'TC'
#'   for total count).
#' @param UQ.pc      Percentile for upper-quartile normalization, number
#' between 0 and 1 (in case the default 0.75 - hence the name - is not
#' appropriate).
#' @param log.transformed  A logical indicating whether expression data were
#'   already log2-transformed, e.g., some microarray data.
#' @param min.LR.found  The minimum number of ligands or receptors found in
#'   \code{count} row names after eliminating the rows containing too many
#'   zeros according to \code{min.count} and \code{prop}.
#' @param conversion.dict  Correspondence table of HUGO gene symbols
#' human/nonhuman. Not used unless the organism is not human.
#'
#' @return A BSRModelData object with empty model parameters.
#'
#' @details The \code{counts} matrix or table should be provided with expression
#'   levels of protein coding genes in each samples (column) and
#'   \code{rownames(counts)} set to HUGO offcial gene symbols. For commodity, it
#'   is also possible to provide \code{counts} with the
#'   gene symbols stored in one of its columns. This column must be specified
#'   with \code{symbol.col}. In such a case, \code{prepareDataset} will extract
#'   this column and use it to set the row names. Because row names must be
#'   unique, \code{prepareDataset} will eliminate rows with duplicated gene
#'   symbols by keeping the rows with maximum average expression. Gene symbol
#'   duplication may occur in protein coding genes after genome alignment
#'   due to errors in genome feature annotation files (GTF/GFF), where a handful
#'   of deprecated gene annotations might remain, or
#'   some genes are not given their fully specific symbols. If your read count
#'   extraction pipeline does not take care of this phenomenon, the maximum mean
#'   expression selection strategy implemented here should solve this difficulty
#'   for the sake of inferring ligand-receptor interactions.
#'
#'   If \code{normalize} is \code{TRUE} then normalization is performed
#'   according to \code{method}. If those two simple methods are not satisfying,
#'   then it is possible to provide a pre-normalized matrix setting
#'   \code{normalize} to \code{FALSE}.
#'
#'   In case proteomic or microarray data are provided, \code{min.count} must be
#'   understood as its equivalent with respect to those data.
#' @export
#' @examples
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
prepareDataset <- function(counts, normalize = TRUE, symbol.col = NULL, min.count = 10,
    prop = 0.1, method = c("UQ", "TC"), log.transformed = FALSE, min.LR.found = 80, 
    species = "hsapiens",conversion.dict = data.frame(Gene.name="A",row.names = "B"),
     UQ.pc = 0.75 ) {

 

    if (prop < 0 || prop > 1)
        stop("prop must lie in [0;1]")
    if (min.count < 0)
        stop("min.count must be positive")
    if (normalize)
        method <- match.arg(method)
    else
        if (nchar(method) == 0)
            stop(paste0("In case of user-normalized counts, the name of the ",
            "normalization must be documented through the parameter 'method'"))

    if (!is.null(symbol.col)) {
        if (!is.numeric(symbol.col))
            stop("symbol.col must be the index of the column containing the gene symbols")

        # simple but desperately slow counts <-
        # aggregate(.~symbol,data=counts,FUN=max)

        # home-made but fast
        symbols <- as.character(counts[, symbol.col])
        d <- symbols[duplicated(symbols)]
        bad <- NULL
        for (s in d) {
            i <- which(symbols == s)
            t <- rowSums(counts[i, -symbol.col])
            bad <- c(bad, i[-which.max(t)])
        }

        # remove duplicates and the gene symbol column
        if (!is.null(bad)){
            counts <- counts[-bad, -symbol.col]
            rownames(counts) <- symbols[-bad]
        }
        else{
            counts <- counts[,-symbol.col]
            rownames(counts) <- symbols
        }
    }

    if (is.null(rownames(counts)) || typeof(rownames(counts)) != "character")
        stop("The read count matrix must be provided with gene symbols as row names")

    # as of now we ensure that counts is a matrix
    if (!is.matrix(counts))
        counts <- data.matrix(counts)

    # avoid empty rows even if no normalization is performed here
    counts <- counts[rowSums(counts) > 0, ]

    if (normalize) {
        good.c <- rowSums(counts >= min.count) >= prop * ncol(counts)
        counts <- counts[good.c, ]
        if (method == "UQ"){
            tot <- apply(counts, 2, function(x) stats::quantile(x[x > 0],
                                                                prob=UQ.pc))
            if (sum(tot == 0) > 0)
                stop(paste0("Cannot perform UQ normalization (percentile=",
                            UQ.pc," ), not enough signal in sample(s) ",
                            paste(colnames(counts)[tot==0], collapse=", ")))
        }
        else
            tot <- colSums(counts)
        ncounts <- sweep(counts, 2, tot/stats::median(tot), "/")
    }
    else
        ncounts <- counts

    homolog.genes <- list()
    if (species != "hsapiens"){
          ncounts <- as.data.frame(ncounts) 
          ncounts$human.gene.name <- rownames(ncounts)
          conversion.dict$human.gene.name  <- rownames(conversion.dict) 
          ncounts$id <- 1:nrow(ncounts) 

          counts.transposed <- merge(ncounts, conversion.dict,
                                     by.x='human.gene.name',
                                     all=FALSE, sort=FALSE)
          counts.transposed <- counts.transposed[order(counts.transposed$id), ]

          homolog.genes <- list(counts.transposed$Gene.name)
          
          counts.transposed$id <- NULL
          ncounts$id <- NULL
          ncounts$human.gene.name <- NULL
          ncounts <- data.matrix(ncounts) 
          rm(counts.transposed)
    }
       
    nLR <- length(intersect(
        c(LRdb$ligand, LRdb$receptor),
        rownames(ncounts)))
    if (nLR < min.LR.found)
        stop(paste0("Not enough LR genes (",nLR," < ", min.LR.found,
                    " were found).\n"))

    new("BSRDataModel", ncounts=ncounts, log.transformed=log.transformed,
        normalization=toupper(method),initial.organism=species,
        initial.orthologs=homolog.genes)

}  # prepareDataset


#' @title Orthologs Gene Names 
#'
#' @description By default, BulkSignalR is designed to work with Homo sapiens.
#' In order to work with other organisms, gene names need to be first converted
#' to human following an orthology mapping process.
#' @param from_organism    An organism defined as in Ensembl: 
#' drerio, mmusculus, celegans, dmelanogaster, etc. This is the source organism 
#' from which you want to convert the gene names to Homo sapiens.
#' @param from_values   A vector of gene names from the current species studied.
#' @param method  Ortholog mapping method.
#' @return Return a data frame with 2 columns containing the gene names
#' for the two species.  
#' First column is the gene name from the source organism 
#' and the second column corresponds to the  homologous gene name
#' in  Homo sapiens.
#' @importFrom orthogene convert_orthologs
#'
#' @export
#' @examples
#' print('findOrthoGenes')
#' data(bodyMap.mouse)
#' ortholog.dict    <- findOrthoGenes (from_organism = "mmusculus", 
#'                                     from_values = rownames(bodyMap.mouse))
#'
findOrthoGenes<- function(from_organism ="mmusculus",
        from_values=c("TP53"),
        method = c("gprofiler","homologene","babelgene")) {

        method <- match.arg(method)
        if (!method  %in% c("gprofiler","homologene","babelgene"))
                stop("Method selected should be gprofiler,homologene or babelgene")

          orthologs_dictionary <- orthogene::convert_orthologs(gene_df = from_values,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = from_organism,
                                        output_species = "human",
                                        non121_strategy = "drop_both_species", # assure 1.1 
                                        method = method,
                                        verbose = FALSE) 
           
          orthologs_dictionary$index <- NULL  
          names(orthologs_dictionary)[1] <- paste("Gene.name")
    
    cat("Dictionary Size: ", 
        dim(orthologs_dictionary)[1],
         " genes \n", sep="") 

    nL <- length(intersect(
        LRdb$ligand,
        rownames(orthologs_dictionary)) )
    cat("-> ",nL, " : Ligands \n", sep="") 

    nR <- length(intersect(
        LRdb$receptor,
        rownames(orthologs_dictionary))) 
    cat("-> ", nR, " : Receptors \n", sep="") 
      
    orthologs_dictionary


} #findOrthoGenes 


#' @title Transpose to Human Gene Names
#'
#' @description By default, BulkSignalR is designed to work with Homo sapiens.
#' In order to work with other organisms, gene names need to be first converted
#' to human following an orthology mapping process.
#' @param counts     A table or matrix of read counts.
#' @param dictionary   A data frame where first column belong to 
#'  organism of study & row names are the human gene names.
#'
#' @return Return a counts matrix transposed for Human.
#'
#' @export
#' @examples
#' print('convertToHuman')
#' data(bodyMap.mouse)
#' 
#' ortholog.dict    <- findOrthoGenes (from_organism = "mmusculus", 
#'                                     from_values = rownames(bodyMap.mouse))
#' 
#' matrix.expression.human <- convertToHuman(counts = bodyMap.mouse,   
#' dictionary = ortholog.dict)
#'
convertToHuman <- function(counts,dictionary=data.frame(Gene.name="A",row.names = "B")) {

          # Should test counts have rownames.
          if(all(row.names(counts)==seq(1, nrow(counts))))
            stop("Rownames should be set as human gene names for counts.", call. = FALSE)
         if(all(row.names(dictionary)==seq(1, nrow(dictionary))))
            stop("Rownames should be set ashuman gene names dictionary.", call. = FALSE)
          if(dim(dictionary)[2]!=1)
            stop("Unique column must be set for dictionary.", call. = FALSE)
         if(! all(apply(counts, 2, function(x) is.numeric(x)))) 
            stop("Some variables are not defined as numerics.", call. = FALSE)

          # Transform Matrice using orthologs_dictionary
          counts$Gene.name  <- rownames(counts)
          dictionary$human.gene.name  <- rownames(dictionary) 
      
          counts$id <- 1:nrow(counts) 

          counts.transposed <- merge(counts,dictionary, by.x='Gene.name',all=FALSE,sort=FALSE)
          counts.transposed <- counts.transposed[order(counts.transposed$id), ]
          counts.transposed$id <- NULL

          # aesthetics only
          counts.transposed <-counts.transposed[,c(which(colnames(counts.transposed)=="human.gene.name"),which(colnames(counts.transposed)!="human.gene.name"))]
          #counts.transposed <-counts.transposed[c("human.gene.name", setdiff(names(counts.transposed), "human.gene.name"))]

          counts.transposed$Gene.name <- NULL

          rownames(counts.transposed) <- counts.transposed[,1]
          counts.transposed           <- counts.transposed[,-1]
        
          counts.transposed

 } # convertToHuman
