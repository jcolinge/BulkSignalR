.onLoad <- function(...) {

 # Initially I was using envir = .GlobalEnv
 # Here and also in dataPrepare::resetLRdb 
 # But it was causing error when running 
 # examples during R CMD check.

 myEnv <- new.env(parent = globalenv())
 attach(myEnv, name="LRdbEnv")
 assign("LRdb", SingleCellSignalR::LRdb, envir = as.environment("LRdbEnv"))

}