####################################################
###     Create / Get Ressources                  ###
####################################################

#' Create all resources
#'
#' Create a cache for all resources (pathways, lrdb & network)
#' downloaded from the web when the library is first loaded.
#' This functionality is handled with BiocFileCache.
#'
#' @param onRequest logical TRUE if you want to force
#' downloading again. This will overwrite the
#' pre-existing local database. Default is TRUE.
#' @param verbose Default is FALSE
#' @return Returns `NULL`, invisibly. 
#' @importFrom curl has_internet
#' @importFrom cli cli_alert_danger cli_alert
#' @export
#' @examples
#' createResources(onRequest=FALSE)
createResources <- function(onRequest = TRUE, verbose = FALSE) {
    cacheDir <- .SignalR$BulkSignalR_CACHEDIR
    resourcesCacheDir <- paste(cacheDir, "resources", sep = "/")

    hasInternet <- tryCatch(expr={curl::has_internet()}, 
        error = FALSE)
    
    if (!hasInternet & 
    !dir.exists(resourcesCacheDir)) {
        cli::cli_alert_danger("Your internet connection is off:")
        stop(
        "- Remote resources cannot be downloaded.\n"
        )   
    }

    if (!hasInternet & 
        onRequest) {
        cli::cli_alert_danger("Your internet connection is off:")
        stop(
        "- Remote resources cannot be downloaded.\n"
        )   
    }

    # Do it once, onLoad
    if (!dir.exists(resourcesCacheDir) | onRequest) {
        .cacheAdd(fpath = .SignalR$BulkSignalR_GO_URL,
            cacheDir = resourcesCacheDir,
            resourceName = "GO-BP", 
            verbose = verbose, download = TRUE)
        .cacheAdd(fpath = .SignalR$BulkSignalR_Reactome_URL,
            cacheDir = resourcesCacheDir, 
            resourceName = "Reactome",
            verbose = verbose, download = TRUE)
        .cacheAdd(fpath = .SignalR$BulkSignalR_Network_URL,
            cacheDir = resourcesCacheDir, 
            resourceName = "Network",
            verbose = verbose, download = TRUE)
        .cacheAdd(fpath = .SignalR$BulkSignalR_LRdb_URL,
            cacheDir = resourcesCacheDir, 
            resourceName = "LRdb",
            verbose = verbose, download = TRUE)

    }

    cacheVersion(dir="resources")
    
    return(invisible(NULL))
}


#' Get resources from the cache
#'
#' Get resources (pathways, lrdb or network
#' stored in the cache.
#'
#' @param resourceName   Resource name.
#' @param cache   Logical. Default is FALSE
#' If TRUE, you will use environment variables.
#' @return Returns a data frame containing the requested
#' resource.
#' @importFrom cli cli_alert_danger
#' @export
#' @examples
#' reactome <- getResource(resourceName = "Reactome",cache=TRUE)
getResource <- function(resourceName = NULL, cache = FALSE) {
    if (!resourceName %in% c("GO-BP", "Reactome", "Network","LRdb")) {
        cli::cli_alert_danger(
        ".val {GO-BP, Reactome, LRdb & Network} are the only keywords alllowed.\n"
        )
        stop()
    }

    if (cache == TRUE) {
        cacheDir <- .SignalR$BulkSignalR_CACHEDIR
        resourcesCacheDir <- paste(cacheDir, "resources", sep = "/")

        # safeguard
        if (!dir.exists(resourcesCacheDir)) {
            cli::cli_alert_danger("Resources repository does not exist.\n")
            stop()
        }

        bfc <- BiocFileCache::BiocFileCache(resourcesCacheDir, ask = FALSE)

        if(resourceName=="LRdb") {
            cacheHits <- BiocFileCache::bfcquery(bfc, query = resourceName,
            field = "rname")
            dataframe <- utils::read.csv(bfc[[cacheHits$rid]],
            stringsAsFactors = FALSE, sep = "\t", check.names = FALSE) 
        } else {
            dataframe <- .readRDSFromCache(bfc = bfc,
            resourceName = resourceName)
        }
        
        if (resourceName == "LRdb") {
            if (!all(
                c("ligand","receptor") %in% 
                colnames(dataframe))) {
                cli::cli_alert_danger("Colnames are not well defined.\n")
                stop()
            }
        }
        if (resourceName == "Reactome") {
            if (!all(
                c("Reactome ID", "Gene name", "Reactome name") %in% 
                colnames(dataframe))) {
                cli::cli_alert_danger("Colnames are not well defined.\n")
                stop()
            }
            dataframe <- dataframe[, 
            c("Reactome ID", "Gene name", "Reactome name")]
        }

        if (resourceName == "GO-BP") {
            if (!all(
                c("GO ID", "Gene name", "GO name") %in% 
                colnames(dataframe))) {
                cli::cli_alert_danger("Colnames are not well defined.\n")
                stop()
            }
            dataframe <- dataframe[, c("GO ID", "Gene name", "GO name")]
        }

        if (resourceName == "Network") {
            if (!all(
                c("a.gn", "type", "b.gn") %in% colnames(dataframe))) {
                cli::cli_alert_danger("Colnames are not well defined.\n")
                stop()
            }
            dataframe <- dataframe[, c("a.gn", "type", "b.gn")]
        }
    } else if (cache == FALSE) {
        if (resourceName == "Reactome") {
            dataframe <- .SignalR$BulkSignalR_Reactome
        }

        if (resourceName == "GO-BP") {
            dataframe <- .SignalR$BulkSignalR_Gobp
        }

        if (resourceName == "Network") {
            dataframe <- .SignalR$BulkSignalR_Network
        }

        if (resourceName == "LRdb") {
            dataframe <- .SignalR$BulkSignalR_LRdb
        }
    }

    return(dataframe)
}


####################################################
###     / Parse / Format / Import Ressources     ###
####################################################

#' Import a refernce network of your own
#'
#' Network is a data frame that defines interactions between
#' genes. It's composed of 3 columns named as
#' follows:
#'
#' a.gn: Gene Symbol 1
#' type: controls-expression-of
#' b.gn: Gene Symbol 2
#'
#' When the user provides his own network,
#' 'type' should be set to 'controls-expression-of'.
#'
#' @param network   Network data frame made of 3 columns
#' a.gn, b.gn & type. 'a.gn' & 'b.gn' should be gene symbols
#' of gene interactions. 'type'  should be set as
#' 'controls-expression-of' when a user provides
#' his own file.
#'
#' @return Returns `NULL`, invisibly. 
#'
#' @importFrom cli cli_alert_info
#' @export
#' @examples
#'  BulkSignalR_Network <- getResource(resourceName = "Network",
#'  cache = FALSE)
#' resetNetwork(BulkSignalR_Network)
resetNetwork <- function(network) {
    if (!all(c("a.gn", "type", "b.gn") %in% colnames(network))) {
        stop("Column names of network should be defined as a.gn, type & b.gn.")
    }

    message("")
    cli::cli_alert_info("New resource defined for {.val Network}.\n")

    assign("BulkSignalR_Network", 
        network, envir = .SignalR)

    return(invisible(NULL))
} # resetNetwork


#' Import pathways from a file or data frame
#'
#' \code{resetPathways} is a function
#' we provide to users who want to refresh REACTOME
#' and GO-BP content included in BulkSignalR.
#' 
#' Pathways in `BulkSignalR` (as sets of genes/proteins) are defined
#' after Reactome and GOBP databases.
#' Those can be updated using
#' json files from
#' the Human Molecular Signatures Database (MSigDB)
#' at \url{https://www.gsea-msigdb.org/}
#' Gmt file format also can be imported.
#' A data frame can be used directly also.
#'
#' @param dataframe  Data frame formated as follows.
#' When \code{resourceName} is set to "Reactome",
#' dataframe colnames must be defined as
#' "Reactome ID", "Gene name", and "Reactome name"
#' When \code{resourceName} is set to "GO-BP",
#' #' dataframe colnames must be defined as
#' "GO ID", "Gene name", and "GO name"
#' @param file    Path to file.
#' @param fileType    Default is Json.
#' Other options are gmt or txt files.
#' @param resourceName    Two options "GO-BP" or "Reactome".
#'
#' @return Returns `NULL`, invisibly. 
#'
#' @importFrom cli cli_alert_info
#' @export
#' @examples
#' reactSubset <- getResource(resourceName = "Reactome",
#' cache = TRUE)
#' 
#' subset <- c("REACTOME_BASIGIN_INTERACTIONS",
#' "REACTOME_SYNDECAN_INTERACTIONS",
#' "REACTOME_ECM_PROTEOGLYCANS",
#' "REACTOME_CELL_JUNCTION_ORGANIZATION")
#' 
#' reactSubset <- reactSubset[
#' reactSubset$`Reactome name` %in% subset,]
#' 
#' resetPathways(dataframe = reactSubset,
#' resourceName = "Reactome")
resetPathways <- function(
    dataframe = NULL,
    file = NULL,
    fileType = c("json", "gmt", "txt"),
    resourceName = NULL) {

    if (!resourceName %in% c("GO-BP", "Reactome")) {
        stop("GO-BP and Reactome are the only keywords alllowed.")
    }

    if (!is.null(file)){

        fileType <- match.arg(fileType)

        if (!file.exists(file)) {
            stop("This file doesn't exist.")
        }

        if (fileType == "json") {
            dataframe <- .formatPathwaysFromJson(
                file = file,
                resourceName = resourceName
            )
        } else if (fileType == "gmt") {
            dataframe <- .formatPathwaysFromGmt(
                file = file,
                resourceName = resourceName
            )
        } else if (fileType == "txt") {
            dataframe <- .formatPathwaysFromTxt(
                file = file,
                resourceName = resourceName
            )
        } else {
            stop("File format accepted are `json` , `gmt` or `txt` only.")
        }

        message("")
        cli::cli_alert_info("New resource defined for {.val {resourceName}}.\n")

        if (resourceName == "Reactome") {
            assign("BulkSignalR_Reactome", dataframe,
                envir = .SignalR)
        }
        if (resourceName == "GO-BP") {
            assign("BulkSignalR_Gobp", dataframe,
                envir = .SignalR)
        }

    }
    else {

    if (resourceName == "Reactome") {
        if (!all(
            c("Reactome ID", "Gene name", "Reactome name") %in% 
            names(dataframe))) {
            mess <- "Colnames should be defined with specific names."
            cli::cli_alert_danger(mess)
            stop("Three columns must defined as",
                " 'Reactome ID','Gene name','Reactome name'.")
        }
    } else if (resourceName == "GO-BP") {
        if (!all(c("GO ID", "Gene name", "GO name") %in% names(dataframe))) {
            mess <- "Colnames should be defined with specific names."
            cli::cli_alert_danger(mess)
            stop("Three columns must defined as",
            " 'GO ID','Gene name','GO name'.")
        }
    } else {
        cli::cli_alert_danger("Resource name is not well defined.\n")
        stop("")
    }

    if (resourceName == "Reactome") {
        assign("BulkSignalR_Reactome", dataframe,
            envir = .SignalR)
    }
    if (resourceName == "GO-BP") {
        assign("BulkSignalR_Gobp", dataframe,
            envir = .SignalR)
    }
    }
    
    return(invisible(NULL))
} # resetPathways


#' Read dataframe from txt file
#'
#' @param file    Path to a tabular file.
#' @param resourceName    Two options "GO-BP"  "REACTOME".
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
#'
.formatPathwaysFromTxt <- function(
    file,
    resourceName = NULL) {
    db <- utils::read.csv(file,
        stringsAsFactors = FALSE, sep = "\t", check.names = FALSE
    )

    if (resourceName == "Reactome") {
        if (!all(
            c("Reactome ID", "Gene name", "Reactome name") %in% 
            names(db))) {
            mess <- "Colnames should be defined with specific names."
            cli::cli_alert_danger(mess)
            stop("Three columns must defined as",
                " 'Reactome ID','Gene name','Reactome name'.")
        }
    } else if (resourceName == "GO-BP") {
        if (!all(c("GO ID", "Gene name", "GO name") %in% names(db))) {
            mess <- "Colnames should be defined with specific names."
            cli::cli_alert_danger(mess)
            stop("Three columns must defined as",
            " 'GO ID','Gene name','GO name'.")
        }
    } else {
        cli::cli_alert_danger("Resource name is not well defined.\n")
        stop("")
    }

    return(db)
} # .formatPathwaysFromTxt


#' Format a data frame according to json input
#'
#' @param file    Path to file.
#' @param resourceName    Two options "GO-BP" or "REACTOME".
#'
#' @return Data frame with pathwayID, geneName and pathwayName
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @import jsonlite
.formatPathwaysFromJson <- function(
    file,
    resourceName = NULL) {
    data <- jsonlite::read_json(file, simplifyVector = TRUE)
    indexPathway <- NULL

    if (!all(c("exactSource", "geneSymbols") %in% names(data[[1]]))) {
        cli::cli_alert_danger("Json format is invalid.\n")
        stop("Json must follow the standard",
        " from The Molecular Signatures Database (MSigDB).")
    }

    db <- foreach::foreach(
        indexPathway = seq_len(length(data)),
        .combine = "rbind"
    ) %dopar% {
        data.frame(
            pathwayID = rep(data[[indexPathway]]$exactSource,
                length(data[[indexPathway]]$exactSource)),
            geneName = unlist(data[[indexPathway]]$geneSymbols)[
            seq_len(length(unlist(data[[indexPathway]]$geneSymbols)))],
            pathwayName = rep(names(data)[[indexPathway]], 
                length(names(data)[[indexPathway]]))
        )
    }

    # Due to the fact React and Go are organized differently
    if (resourceName == "Reactome") {
        names(db) <- c("Reactome ID", "Gene name", "Reactome name")
    }

    if (resourceName == "GO-BP") {
        names(db) <- c("GO ID", "Gene name", "GO name")
    }

    return(db)
} # .formatPathwaysFromJson


#' Transform gmt file to data frame
#'
#' We note discrepancies between the formats available
#' from internet sources. Here, we consider a valid gmt file format defined
#' on each lines as follows:
#' First is Pathway name,
#' then comes the ID,
#' finally you will find genes symbols
#' part of the pathway defined on the line.
#'
#' You can find an example here.
#' - For Reactome. (Directly from their website)
#' \url{https://reactome.org/download/current/ReactomePathways.gmt.zip}
#' Note that you need to unzip the file to read the content.

#' The code is inspired from read.gmt function
#' from the gsa R package.
#'
#' @param file    Path to GMT file
#' @param resourceName   Two options "GO-BP" or "REACTOME"
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
.formatPathwaysFromGmt <- function(
    file,
    resourceName = NULL) {
    read1 <- scan(file,
        what = list("", ""), sep = "\t",
        quote = NULL, fill = TRUE, flush = TRUE, multi.line = FALSE
    )

    read1 <- stats::setNames(read1, c("Description", "ID"))
    read1 <- read1[c("ID", "Description")]

    geneset.ids <- read1[1][[1]]
    geneset.descriptions <- read1[2][[1]]

    read2 <- scan(file, what = "", sep = "\t", quote = NULL)
    read2 <- read2[read2 != ""]

    nn <- length(geneset.ids)
    n <- length(read2)
    ox <- rep(NA, nn)

    # Compute indice of pathway
    ii <- 1
    for (i in seq_len(nn)) {
        while ((read2[ii] != geneset.descriptions[i]) | 
            (read2[ii + 1] != geneset.ids[i])) {
            ii <- ii + 1
        }
        ox[i] <- ii
        ii <- ii + 1
    }

    genesets <- vector("list", nn)

    for (i in seq_len(nn - 1)) {
        i1 <- ox[i] + 2
        i2 <- ox[i + 1] - 1
        geneset.descriptions[i] <- read2[ox[i] + 1]
        geneset.ids[i] <- read2[ox[i]]
        genesets[[i]] <- read2[i1:i2]
    }

    geneset.ids[nn] <- read2[ox[nn]]
    geneset.descriptions[nn] <- read2[ox[nn] + 1]
    genesets[[nn]] <- read2[(ox[nn] + 2):n]

    data <- list(
        geneset.ids = geneset.ids,
        geneset.descriptions = geneset.descriptions,
        genesets = genesets
    )

    dataframeFromGmt <- foreach::foreach(
        i = seq_len(length(data$geneset.ids)),
        .combine = "rbind"
    ) %dopar% {
        data.frame(
            a = rep(data$geneset.ids[[i]], length(data$genesets[[i]])),
            b = data$genesets[[i]][seq_len(length(data$genesets[[i]]))],
            c = rep(data$geneset.descriptions[[i]], length(data$genesets[[i]]))
        )
    }

    # Due to the fact React and Go are organized differently
    if (resourceName == "Reactome") {
        names(dataframeFromGmt) <- 
        c("Reactome name", "Gene name", "Reactome ID")
        dataframeFromGmt <- dataframeFromGmt[
        c("Reactome ID", "Gene name", "Reactome name")]
    }

    if (resourceName == "GO-BP") {
        names(dataframeFromGmt) <- c("GO name", "Gene name", "GO ID")
        dataframeFromGmt <- dataframeFromGmt[c("GO ID", "Gene name", "GO name")]
    }

    return(dataframeFromGmt)
} # .formatPathwaysFromGmt

#' Modify LRdb database
#'
#' Users can provide a data frame with 2 columns named
#' ligand and receptor.
#' This can be used to extend or replace the existing
#' LRdb.
#'
#' @param db     A data frame with 2 columns named
#' ligand and receptor.
#' @param switch  A logical indicating whether LRdb should be extended only
#' (FALSE, default) or completely replaced (TRUE).
#'
#' @return Returns `NULL`, invisibly. 
#'
#' @importFrom cli cli_alert_info
#' @export
#' @examples
#' resetLRdb(db = data.frame(ligand = "A2M", receptor = "LRP1"), switch = FALSE)
resetLRdb <- function(db, switch = FALSE) {
    if (all(c("ligand","receptor") %in% names(db))) {
        if (switch) {
            assign("BulkSignalR_LRdb", unique(db[, c("ligand", "receptor")]),
                envir = .SignalR
            )
        } else {
            db <- rbind(
                .SignalR$BulkSignalR_LRdb[, c("ligand", "receptor")],
                db[, c("ligand", "receptor")]
            )
            assign("BulkSignalR_LRdb", unique(db), 
                envir = .SignalR)
        }
    } else {
        stop(
            "db should be a dataframe with",
            "2 columns named 'ligand' and 'receptor'."
        )
    }

    message("")
    cli::cli_alert_info(
        "New database defined for {.val LRdb}."
    )

    return(invisible(NULL))
    
} # resetLRdb
