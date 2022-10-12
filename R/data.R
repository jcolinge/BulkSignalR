#' Salivary duct carcinoma transcriptomes
#'
#' A dataset containing the read counts of salivary duct carcinomas
#'   and adjacent normal tissues.
#'
#' @format A data frame with 19764 rows and 26 variables.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138581}
#' @usage data(sdc)
"sdc"

#' Mouse transcriptomes across tissues
#'
#' A dataset containing rpkm values of  brain and liver.
#'   
#'
#' @format A data frame with 24543 rows and 8 variables.
#' @source Bin Li & al., Scientific Reports, 2017;
#' @usage data(bodyMap.mouse)
"bodyMap.mouse"


#' Tumor microenvironment gene signatures
#'
#' A dataset containing gene signatures for some immune and stromal
#'   cell populations that are present in the microenvironment of a tumor.
#'
#' @format A data frame with 209 rows and 2 variables: \describe{
#'   \item{gene}{HUGO gene symbol} \item{signature}{cell population name} }
#' @source Becht & al., Genome Biol, 2016; Angelova et al., Genome Biol, 2015.
#' @usage data(tme.signatures)
"tme.signatures"


#' Partial EMT gene signature
#'
#' A dataset containing a partial EMT gene signature.
#'
#' @format A data frame with 100 rows and 1 variables:
#' \describe{
#'   \item{gene}{HUGO gene symbol}
#' }
#' @source Puram, SV & al., Cell, 2017.
#' @usage data(p.EMT)
"p.EMT"


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
#' @usage data(immune.signatures)
"immune.signatures"
