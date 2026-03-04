#' @title Predict Cortical DNA Methylation Clock Age (2020)
#'
#' @description Predicts DNAm age in cortical samples using the elastic net 
#' model by Shireby et al. (2020).
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector of predicted cortical DNAm ages for 
#' each sample.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate 
#' multiple clocks simultaneously.
#'
#' @references
#' Shireby GL, Davies JP, Francis PT, et al.
#' Recalibrating the epigenetic clock: implications for assessing 
#' biological age in the human cortex.
#' \emph{Brain.} 2020
#'
#' @export
#'
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' corticalClockOut <- corticalClock(hannum_bmiq_m)


corticalClock <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
  # --- Step 1: Load and parse coefficients (from package internal data) ---
  data("CorticalClockCoef", envir = environment())
  coefData <- CorticalClock_coef
  refData  <- CorticalClock_ref 
  clockName <- "corticalClock"
  
  # --- 2. Imputation with Ref ---
  clockProbes <- as.character(coefData[-1, 1])
  userProbes  <- rownames(betaM)
  missingProbes <- setdiff(clockProbes, userProbes)
  
  betaComplete <- betaM
  
  if (length(missingProbes) > 0) {
    if (verbose) message(paste0("[", clockName, "] Imputing ", 
                                length(missingProbes), 
                                " missing probes using reference data."))
    
    # Check whether the reference data is complete
    if (!all(missingProbes %in% names(refData))) {
      stop(paste(clockName, 
                 "Error: Reference data is missing for some required probes."))
    }
    
    # Construct the completion matrix: 
    # Rows are missing probes and columns are samples
    
    refVals <- refData[missingProbes]
    refMatrix <- matrix(rep(refVals, ncol(betaM)), 
                        nrow = length(missingProbes), 
                        ncol = ncol(betaM))
    rownames(refMatrix) <- missingProbes
    colnames(refMatrix) <- colnames(betaM)
    
    
    betaComplete <- rbind(betaM, refMatrix)
  }
  
  # --- 3. Prediction ---
  predAgev <- .calLinearClock(
    betaM       = betaComplete,
    coefData    = coefData,
    clockLabel  = clockName,
    minCoverage = minCoverage,
    verbose     = verbose
  )
  
  predAgev <- .antiTrafo(predAgev)
  
  return(predAgev)
}




