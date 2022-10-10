.onLoad <- function(...) {

 # envir = package.BulkSignalR
 assign("LRdb", SingleCellSignalR::LRdb, envir = .GlobalEnv)

  # envir = package.BulkSignalR
 assign("PwC_ReactomeKEGG", 
    SingleCellSignalR::PwC_ReactomeKEGG, envir = .GlobalEnv)
}