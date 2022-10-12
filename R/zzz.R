.onLoad <- function(...) {

 # envir = package.BulkSignalR .GlobalEnv
 assign("LRdb", SingleCellSignalR::LRdb, envir = .GlobalEnv)

}