.onLoad <- function(...) {

 # envir = package.BulkSignalR
 assign("LRdb", SingleCellSignalR::LRdb, envir = .GlobalEnv)

}