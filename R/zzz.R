.onLoad <- function(...) {

 # package.BulkSignalR
 assign("LRdb", SingleCellSignalR::LRdb, envir = .GlobalEnv)
 
}