#' Salivary duct carcinoma transcriptoms
#'
#' A data set containing the read counts of salivary duct carcinomas (SDCs)
#'   and adjacent normal tissues.
#'
#' @format A data frame with 19764 rows and 26 variables.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138581}
#' @usage data(sdc)
"sdc"

#' A skinny BSRDataModel object related to sdc.
#'
#' Output from the `learnParameters` function to get
#' BulkSignalR statistical model parameters.
#' 
#' @format An example of an object created by `BSRDataModel`
#' applied to an sdc subset (Patients N20,N22,SDC17,SDC25) and
#' 10 000 genes sampled (seed set to 123)
#' `learnParameters` was also called to get statistical 
#' model parameters.
#' @usage data(bsrdm)
"bsrdm"

#' A skinny BSRInference object related to sdc.
#'
#' From the previous object `bsrdm`, 
#' you can generate inferences by calling its
#' constructor `BSRInference`. 
#' The resulting BSRInference object is `bsrinf`,
#' It contains all the inferred L-R interactions 
#' with their associated pathways and corrected p-values.
#' 
#' @format An example of an object created by inference function 
#' @usage data(bsrinf)
"bsrinf"

#' A skinny BSRDataModelComp object related to sdc.
#'
#' See Vignette BulkSignalR-Differential.
#' 
#' @format An example of an BSRDataModelComp object
#' @usage data(bsrdm.comp)
"bsrdm.comp"

#' A skinny BSRInferenceComp object related to sdc.
#'
#' See Vignette BulkSignalR-Differential.
#' 
#' @format An example of an BSRInferenceComp object
#' @usage data(bsrinf.comp)
"bsrinf.comp"

#' A skinny BSRDataModel object related to a spatial data set
#' 
#' Obtained from STexampleData::Visium_humanDLPFC. 
#' A single sample (sample 151673)
#' of human brain dorsolateral prefrontal cortex (DLPFC)
#' in the human brain, measured using the 10x Genomics 
#' Visium platform. This is a subset of the full dataset 
#' published by Maynard and Collado-Torres et al. (2021).
#' The subset is reproduced in the vignette.
#' name.idx <- c("10x32","3x47","4x50",
#' "17x111","5x59","0x20","8x100",
#' "8x108","14x30","11x39") 
#' 
#' Output from the `learnParameters` function to get
#' BulkSignalR statistical model parameters for a subset
#' of a spatial data set.
#' 
#' @format An example of an object created by `BSRDataModel`
#' applied to a subset of a spatial data set.
#' `learnParameters` was also called to get statistical 
#' model parameters.
#' @source \url{http://spatial.libd.org/spatialLIBD/}
#' @usage data(bsrdm.spa)
"bsrdm.spa"

#' A skinny BSR-inference object related to spatial data set
#'
#' Output from the `learnParameters` function to get
#' BulkSignalR statistical model parameters.
#' 
#' @format From the previous object `bsrdm.spa`, 
#' you can generate inferences by calling its
#' method `BSRInference`. 
#' The resulting BSRInference object is `bsrinf.spa`,
#' It contains all the inferred L-R interactions 
#' with their associated pathways and corrected p-values.
#' `learnParameters` was also called to train the statistical 
#' model parameters.
#' @source \url{http://spatial.libd.org/spatialLIBD/}
#' @usage data(bsrinf.spa)
"bsrinf.spa"

#' A skinny data frame used in the spatial workflow
#'
#' Data frame subset describing the spatial spots
#' 
#' @format Data frame that contains the following columns:
#' barcode_id,sample_id, in_tissue,array_row
#' array_col, ground_truth, reference, cell_count, idSpatial
#' 
#' barcode_id is the id of the spot
#' idSpatial is the spatial id of the spot(array_rowXarray_col)
#' ground_truth is the label (Layer1/2 were only kept)
#' 
#' They are the mandatory data necessary to
#' generate plots for the spatial workflow.
#' 
#' @source \url{http://spatial.libd.org/spatialLIBD/}
#' @usage data(annotation.spa)
"annotation.spa"

#' Mouse transcriptoms across tissues
#'
#' A data set containing RPKM values of brain and liver.
#'
#'
#' @format A data frame with 24543 rows and 8 variables.
#' @source Bin Li & al., Scientific Reports, 2017;
#' @usage data(bodyMap.mouse)
"bodyMap.mouse"

#' A skinny BSRInference object related to bodyMap.mouse
#' 
#' see related workflow for non human organism
#' in the vignette
#'
#' @format An example of an object created by inference function
#' @usage data(bsrinf.mouse)
"bsrinf.mouse"

#' A skinny data frame used in the mouse workflow
#'
#' Synthetic object used during the call to the
#' function `resetToInitialOrganism``
#'
#' @format An example of a data frame created by findOrthoGenes
#' @usage data(ortholog.dict)
"ortholog.dict"

#' Tumor microenvironment gene signatures
#'
#' A data set containing gene signatures for some immune and stromal
#'   cell populations that are present in the microenvironment of a tumor.
#'
#' @format A data frame with 209 rows and 2 variables: \describe{
#'   \item{gene}{HUGO gene symbol} \item{signature}{cell population name} }
#' @source Becht & al., Genome Biol, 2016; Angelova et al., Genome Biol, 2015.
#' @usage data(tme.signatures)
"tme.signatures"


#' Partial EMT gene signature
#'
#' A data set containing a partial EMT gene signature.
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
#' A data set containing gene signatures for general immune cell populations.
#'
#' @format A data frame with 1541 rows and 2 variables:
#' \describe{
#'   \item{gene}{HUGO gene symbol}
#'   \item{signature}{cell population name}
#' }
#' @source PanglaoDB (Franzén et al., Database, 2019).
#' @usage data(immune.signatures)
"immune.signatures"
