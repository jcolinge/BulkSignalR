#' Salivary duct carcinoma transcriptomes
#'
#' A dataset containing the read counts of salivary duct carcinomas
#'   and adjacent normal tissues.
#'
#' @format A data frame with 19764 rows and 26 variables.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138581}
#' @usage data(sdc)
"sdc"


#' Tumor microenvironment gene signatures
#'
#' A dataset containing gene signatures for some immune and stromal
#'   cell populations that are present in the microenvironment of a tumor.
#'
#' @format A data frame with 209 rows and 2 variables: \describe{
#'   \item{gene}{HUGO gene symbol} \item{signature}{cell population name}git }
#' @source Becht et al., Genome Biol, 2016; Angelova et al., Genome Biol, 2015.
#' @usage data(tme.signatures)
"tme.signatures"


#' Immune cell gene signatures
#'
#' A dataset containing gene signatures for general immune cell populations.
#'
#' @format A data frame with 1541 rows and 2 variables:
#' \describe{
#'   \item{gene}{HUGO gene symbol}
#'   \item{signature}{cell population name}
#' }
#' @source PanglaoDB (Franzén et al., Database, 2019).
"immune.signatures"
