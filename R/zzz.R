#nameEnv <- ".SignalR-Env"
.SignalR<- new.env(parent = emptyenv()) 

.onLoad <- function(...) {
    
    # handle directory creation over different OS
    cacheDir <- tools::R_user_dir("BulkSignalR", which="cache")
    assign("BulkSignalR_CACHEDIR", cacheDir, envir = .SignalR)

    url <- "https://partage-dev.montp.inserm.fr:9192/CBSB/"

    assign("BulkSignalR_CORE_URL", url, envir = .SignalR)

    .testRemoteServer()
    .testCacheFiles()
    

    ################################
    ##   Resource Cache Files   ###
    ################################
    urlGo <- paste0(url,
        "SignalR/resources/gobp.rds")
    urlReactome <- paste0(url,
        "SignalR/resources/reactome.rds")
    urlNetwork <- paste0(url,
        "SignalR/resources/Network.rds")
    urlLRdb <- paste0(url,
        "SignalR/resources/LRdb.txt")

    assign("BulkSignalR_GO_URL", urlGo, 
    envir = .SignalR)
    assign("BulkSignalR_Reactome_URL", urlReactome, 
    envir = .SignalR)
    assign("BulkSignalR_Network_URL", urlNetwork, 
    envir = .SignalR)
    assign("BulkSignalR_LRdb_URL", urlLRdb, 
    envir = .SignalR)

    createResources(onRequest = FALSE)
   
    BulkSignalR_Reactome <- getResource(resourceName = "Reactome",
        cache = TRUE)
    BulkSignalR_Gobp <- getResource(resourceName = "GO-BP",
        cache = TRUE)
    BulkSignalR_Network <- getResource(resourceName = "Network",
        cache = TRUE)
    BulkSignalR_LRdb <- getResource(resourceName = "LRdb",
        cache = TRUE)  

    assign("BulkSignalR_Reactome", BulkSignalR_Reactome, 
    envir = .SignalR)
    assign("BulkSignalR_Gobp", BulkSignalR_Gobp, 
    envir = .SignalR)
    assign("BulkSignalR_Network", BulkSignalR_Network, 
    envir = .SignalR)
    assign("BulkSignalR_LRdb", BulkSignalR_LRdb, 
    envir = .SignalR)
}
