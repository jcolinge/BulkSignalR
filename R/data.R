#' Salivary duct carcinoma transcriptomes
#'
#' A dataset containing the read counts of salivary duct carcinomas and adjacent normal tissues.
#'
#' @format A data frame with 19764 rows and 26 variables.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138581}
"sdc"


#' Tumor microenvironment gene signatures
#'
#' A dataset containing gene signatures for some immune and stromal cell populations
#' that are present in the microenvironment of a tumor.
#'
#' @format A data frame with 209 rows and 2 variables:
#' \describe{
#'   \item{gene}{HUGO gene symbol}
#'   \item{signature}{cell population name}
#' }
#' @source Becht et al., Genome Biol, 2016, and Angelova et al., Genome Biol, 2015.
"tme.signatures"
