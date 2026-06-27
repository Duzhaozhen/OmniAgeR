#' @title Calculate Universal Pan-Mammalian Blood Epigenetic Clocks
#' @description
#' Applies the two universal pan-mammalian blood clocks (Clock 2 and 3)
#' from Lu et al. (2023) to a given dataset of DNA methylation values.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param speciesName A character string or vector specifying the Latin species
#' name(s) (e.g., "Homo sapiens").
#' @param anageData A data.frame containing the AnAge database information.
#'   Must include 'SpeciesLatinName', 'GestationTimeInYears',
#'   'averagedMaturity.yrs', and 'maxAge'.
#' @param minCoverage Numeric (0-1). Minimum required proportion of CpGs present.
#' Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#'
#' @return A data.frame containing the 'Sample', 'SpeciesLatinName',
#'   and the calculated ages: 'DNAmRelativeAge',
#'    'DNAmAgePanMammalianBlood2', 'DNAmRelativeAdultAge' and
#'   'DNAmAgePanMammalianBlood3'. Any additional columns
#'   from the input `sampleInfo` (like 'Age', 'Tissue') will also be returned.
#'
#'
#' @export
#'
#' @references
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#' @examples
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#'   tursiopsExample <- loadOmniAgeRdata(
#'       "omniager_tursiops_example",
#'       verbose = FALSE
#'   )
#'   
#'   # Run the calculation. Note: anageData is automatically loaded if omitted.
#'   clockResults <- panMammalianBlood(
#'       x = tursiopsExample$beta_m,
#'       speciesName = tursiopsExample$PhenoTypes$SpeciesLatinName,
#'       verbose = FALSE
#'   )
#'   print(head(clockResults))
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = tursiopsExample$beta_m),
#'       colData = tursiopsExample$PhenoTypes
#'     )
#'     
#'     se_res <- panMammalianBlood(
#'       x = se_obj,
#'       speciesName = se_obj$SpeciesLatinName,
#'       verbose = FALSE
#'     )
#'   }
#' }
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/shorvath/MammalianMethylationConsortium
# under the MIT License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------

panMammalianBlood <- function(x, speciesName, anageData = NULL, minCoverage = 0, verbose = TRUE) {
  .calculatePanMammalian(
    x = x, speciesName = speciesName, anageData = anageData,
    minCoverage = minCoverage, verbose = verbose,
    coefDataName = "omniager_pan_mammalian_blood_coef",
    yNames = c("Y.pred2", "Y.pred3"),
    outputPrefix = "DNAmAgePanMammalianBlood",
    functionName = "panMammalianBlood"
  )
}