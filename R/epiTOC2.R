#' @title
#' Estimate epiTOC2 scores
#'
#' @aliases epiTOC2
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an
#' age vector (optional) and will return the epiTOC2 scores.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param age Optional numeric vector representing chronological ages
#'   of the samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' Building upon a dynamic model of DNA methylation gain in unmethylated
#'  CpG-rich regions, epiTOC2 can directly estimate the cumulative number
#'  of stem cell divisions in a tissue. The details of the algorithm are
#'  described in Teschendorff et al. 2020.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using the full epiTOC2 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using an approximation of epiTOC2
#'   which assumes all epiTOC2 CpGs have beta-values exactly 0 in the
#'   fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the
#'   estimated average lifetime intrinsic rate of stem-cell division per
#'   sample, as derived from epiTOC2
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate
#'   for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#'
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @importFrom stats median
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' epiTOC2Out <- epiTOC2(x = beta_matrix, verbose = FALSE)
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
#'     epiTOC2Out <- epiTOC2(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
#'

epiTOC2 <- function(x, age = NULL, minCoverage = 0, verbose = TRUE) {
  .calculateEpiTOC(
    x = x, 
    age = age, 
    minCoverage = minCoverage, 
    verbose = verbose,
    modelName = "omniager_epitoc2_model",
    clockName = "epiTOC2",
    calcAvETOC3 = FALSE
  )
}


#' @title Internal helper function to calculate epiTOC variants
#'
#' @description
#' A shared internal engine designed to compute the cumulative number of stem-cell 
#' divisions using either the epiTOC2 or epiTOC3 models. This function extracts 
#' shared logic to adhere to DRY (Don't Repeat Yourself) principles.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'   object containing DNA methylation beta values. Rows should be CpG probes and 
#'   columns should be individual samples.
#' @param age Optional numeric vector representing the chronological ages of the samples.
#'   Required for calculating intrinsic rates of stem-cell division.
#' @param minCoverage Numeric value between 0 and 1. The minimum required proportion 
#'   of probe coverage.
#' @param verbose Logical. Whether to print coverage statistics and progress messages.
#' @param modelName Character string. The name of the RData model to load 
#'   (e.g., \code{"omniager_epitoc3_model"}).
#' @param clockName Character string. The specific name of the clock, used primarily 
#'   for coverage check reporting (e.g., \code{"epiTOC3"}).
#' @param calcAvETOC3 Logical. If \code{TRUE}, the function will calculate and append 
#'   the average epiTOC3 score (\code{avETOC3}) to the results list. Default is \code{FALSE}.
#'
#' @return A named list containing the following components:
#' \itemize{
#'   \item{\code{tnsc}: Numeric vector. The estimated cumulative number of stem-cell divisions per stem-cell per year using the full model.}
#'   \item{\code{tnsc2}: Numeric vector. The estimated cumulative number of stem-cell divisions using the approximation model (assuming beta0 = 0).}
#'   \item{\code{irS}: Numeric vector or \code{NULL}. The estimated average lifetime intrinsic rate of stem-cell division per sample. Returned only if \code{age} is provided.}
#'   \item{\code{irS2}: Numeric vector or \code{NULL}. As \code{irS}, but for the approximation model.}
#'   \item{\code{irT}: Numeric. The median estimate over all \code{irS} values representing the tissue-level intrinsic rate.}
#'   \item{\code{irT2}: Numeric. As \code{irT}, but for the approximation model.}
#'   \item{\code{avETOC3}: Numeric vector. Included only if \code{calcAvETOC3} is \code{TRUE}. Represents the average over the epiTOC3 specific sites.}
#' }
#' 
#' @importFrom stats median setNames
#' @keywords internal
#' @noRd

.calculateEpiTOC <- function(x, age, minCoverage, verbose, modelName, clockName, calcAvETOC3 = FALSE) {
  betaM <- .extractAssayMatrix(x)
  estParams <- loadOmniAgeRdata(modelName, verbose = verbose)
  
  dummyWeights <- setNames(seq_len(nrow(estParams)), rownames(estParams))
  
  # Perform coverage check
  coverageResult <- .checkCpGCoverage(
    betaM = betaM,
    allWeights = dummyWeights,
    clockName = clockName,
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  # Handle coverage failure
  if (!coverageResult$pass) {
    retList <- list(
      tnsc = rep(NA_real_, ncol(betaM)), 
      tnsc2 = rep(NA_real_, ncol(betaM)),
      irS = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
      irS2 = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
      irT = NA_real_, 
      irT2 = NA_real_
    )
    # Add avETOC3 placeholder if required by epiTOC3
    if (calcAvETOC3) {
      retList$avETOC3 <- rep(NA_real_, ncol(betaM))
    }
    return(retList)
  }
  
  # Extract the matching data and parameters
  matchedParams <- estParams[names(coverageResult$weightsSubset), , drop = FALSE]
  subBeta <- betaM[coverageResult$betaIdx, , drop = FALSE]
  
  deltaV <- matchedParams[, 1]
  beta0V <- matchedParams[, 2]
  
  # Core algorithm implementation
  # Full Model
  scalingFactor <- 2 / (deltaV * (1 - beta0V))
  tnscV <- colMeans(sweep(subBeta, 1, beta0V, "-") * scalingFactor, na.rm = TRUE)
  
  # Approximation (beta0 = 0)
  scalingFactor2 <- 2 / deltaV
  tnsc2V <- colMeans(subBeta * scalingFactor2, na.rm = TRUE)
  
  # Intrinsic Rate
  irS <- NULL
  irS2 <- NULL
  irT <- NULL
  irT2 <- NULL
  
  if (!is.null(age)) {
    irS <- tnscV / age
    irS2 <- tnsc2V / age
    irT <- median(irS, na.rm = TRUE)
    irT2 <- median(irS2, na.rm = TRUE)
  }
  
  # Prepare return list
  retList <- list(
    tnsc = tnscV,
    tnsc2 = tnsc2V,
    irS = irS,
    irS2 = irS2,
    irT = irT,
    irT2 = irT2
  )
  
  # Add avETOC3 if required by epiTOC3
  if (calcAvETOC3) {
    retList$avETOC3 <- colMeans(subBeta, na.rm = TRUE)
  }
  
  return(retList)
}
