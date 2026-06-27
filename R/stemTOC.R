#' @title  Estimate stemTOC score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' will return the stemTOC score.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' The function will return the 0.95 upper quantile of the 371 stemTOC CpGs.
#' Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for
#' significant DNA hypermethylation with chronological age
#' in large in-vivo datasets
#'
#' @return The stemTOC score of each sample.
#'
#' @references
#' Zhu, T., Tong, H., Du, Z. et al.
#' An improved epigenetic counter to track mitotic age in normal
#' and precancerous tissues.
#' \emph{Nat Commun} 2024
#'
#' @importFrom stats quantile
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- stemTOC(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     predOut <- stemTOC(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
stemTOC <- function(x, minCoverage = 0, verbose = TRUE) {
  .calculateStemTOC(
    x = x, 
    minCoverage = minCoverage, 
    verbose = verbose,
    cpgDataName = "omniager_stemtoc_cpg",
    clockName = "stemTOC"
  )
}


#' @title Internal helper function to calculate stemTOC variants
#'
#' @description
#' A shared internal engine designed to compute the stemTOC and stemTOCvitro 
#' scores. It extracts the 0.95 upper quantile of the specified target CpGs.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'   object containing DNA methylation beta values.
#' @param minCoverage Numeric value between 0 and 1. The minimum required proportion 
#'   of overlapping CpGs. Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#' @param cpgDataName Character string. The internal RData identifier for the target 
#'   CpGs to load (e.g., \code{"omniager_stemtoc_cpg"}).
#' @param clockName Character string. The name of the clock for logging and coverage 
#'   checking (e.g., \code{"stemTOC"}).
#'
#' @return A numeric vector representing the calculated score for each sample. 
#'   Returns \code{NA_real_} for samples failing the coverage check.
#'
#' @importFrom stats quantile setNames
#' @keywords internal
#' @noRd
.calculateStemTOC <- function(x, minCoverage, verbose, cpgDataName, clockName) {
  betaM <- .extractAssayMatrix(x)
  targetCpGsData <- loadOmniAgeRdata(cpgDataName, verbose = verbose)
  
  # Prepare the reference probe weights (all 1s for quantile calculation)
  targetCpGs <- as.character(targetCpGsData)
  clockWeights <- setNames(rep(1, length(targetCpGs)), targetCpGs)
  
  # Perform coverage check
  coverageResult <- .checkCpGCoverage(
    betaM = betaM,
    allWeights = clockWeights,
    clockName = clockName,
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  # Handle coverage failure
  if (!coverageResult$pass) {
    scores <- rep(NA_real_, ncol(betaM))
    names(scores) <- colnames(betaM)
    return(scores)
  }
  
  # Calculate the 0.95 quantile score
  scores <- apply(
    betaM[coverageResult$betaIdx, , drop = FALSE],
    2,
    quantile,
    probs = 0.95,
    na.rm = TRUE
  )
  
  return(scores)
}