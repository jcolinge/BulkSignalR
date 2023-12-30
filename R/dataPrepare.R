#' Modify LRdb database
#'
#' User can provide a data frame with 2 columns named
#' ligand and receptor.
#' This can be used to extend or replace the existing
#' LRdb.
#'
#' @param db     A dataframe with 2 columns named
#' ligand and receptor.
#' @param switch  A logical indicating whether LRdb should be extended only
#' (FALSE, default) or completely replaced (TRUE).
#'
#' @return NULL
#'
#' @export
#' @examples
#' print('resetLRdb')
#' data(sdc,package='BulkSignalR')
#' resetLRdb(db=data.frame(ligand="A2M", receptor="LRP1"), switch=FALSE)
resetLRdb <- function(db, switch=FALSE) {

    if(colnames(db)[1]=='ligand' & colnames(db)[2]=='receptor'){
      
        if(switch){
            # envir = .GlobalEnv
            assign("LRdb", unique(db[, c('ligand', 'receptor')]),
                   envir=as.environment("LRdbEnv"))
        }
        else {  
            updated.LRdb <- rbind(LRdb[, c('ligand', 'receptor')],
                                  db[, c('ligand', 'receptor')])
            # envir = .GlobalEnv
            assign("LRdb", unique(updated.LRdb), envir=as.environment("LRdbEnv"))
        }
    } else {
      stop(paste0("db should be a data frame with ",
            "2 columns named 'ligand' and 'receptor'."))
    }
    
} # resetLRdb


#' Prepare a BSRDataModel object from expression data
#'
#' Take a matrix or data frame containing RNA sequencing,
#' microarray, or expression proteomics data and return a BSRDataModel
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
#'   for total count). If \code{normalize==FALSE}, then method must be
#'   used to document the name of the normalization method applied by the user.
#' @param UQ.pc      Percentile for upper-quartile normalization, number
#' between 0 and 1 (in case the default 0.75 - hence the name - is not
#' appropriate).
#' @param log.transformed  A logical indicating whether expression data were
#'   already log2-transformed, e.g., some microarray data.
#' @param min.LR.found  The minimum number of ligands or receptors found in
#'   \code{count} row names after eliminating the rows containing too many
#'   zeros according to \code{min.count} and \code{prop}.
#' @param conversion.dict  Correspondence table of HUGO gene symbols
#' human/nonhuman. Not used unless the organism is different from human.
#'
#' @return A BSRModelData object with empty model parameters.
#'
#' @details The \code{counts} matrix or table should be provided with expression
#'   levels of protein coding genes in each samples (column) and
#'   \code{rownames(counts)} set to HUGO official gene symbols. For commodity, it
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
#'   \code{normalize} to \code{FALSE}. In such a case, the parameter
#'   \code{method} must be used to document the name of the normalization
#'   algorithm used.
#'
#'   In case proteomic or microarray data are provided, \code{min.count} must be
#'   understood as its equivalent with respect to those data.
#' 
#' 
#' @export
#' @examples
#' print('prepareDataset')
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#'
prepareDataset <- function(counts, normalize = TRUE, symbol.col = NULL, min.count = 10,
    prop = 0.1, method = c("UQ", "TC"), log.transformed = FALSE, min.LR.found = 80, 
    species = "hsapiens", conversion.dict = NULL,
    UQ.pc = 0.75) {

    if ((species != "hsapiens") && is.null(conversion.dict))
        stop("Non-human species but no conversion.dict provided")
    
    if (normalize){
        if (prop < 0 || prop > 1)
            stop("prop must lie in [0;1]")
        if (UQ.pc <= 0 || UQ.pc > 1)
            stop("UQ.pc must lie in ]0;1]")
        if (min.count < 0)
            stop("min.count must be positive")
        method <- match.arg(method)
    }
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
    counts <- counts[rowSums(abs(counts)) > 0, ]

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
findOrthoGenes<- function(from_organism, from_values,
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
#' to human by orthology.
#' @param counts  A table or matrix of read counts.
#' @param dictionary   A data frame where the first column is made of gene symbols
#' for the actual organism and row names are the ortholog human gene symbols.
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
convertToHuman <- function(counts, dictionary) {

          # we need counts to be a data.frame
          if (is.matrix(counts)){
            was.matrix <- TRUE
            counts <- as.data.frame(counts)
          }
          else
            was.matrix <- FALSE
          # counts should have row names
          if(all(row.names(counts)==seq(1, nrow(counts))))
            stop("Rownames should be set as human gene names for counts.", call. = FALSE)
         if(all(row.names(dictionary)==seq(1, nrow(dictionary))))
            stop("Rownames should be set as human gene names dictionary.", call. = FALSE)
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

          if (was.matrix)
            # convert back to a matrix
            data.matrix(counts.transposed)
          else
            counts.transposed

 } # convertToHuman
