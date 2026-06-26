#' @title Predict Cortical DNA Methylation Clock Age (2020)
#'
#' @description Predicts DNAm age in cortical samples using the elastic net
#' model by Shireby et al. (2020).
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector of predicted cortical DNAm ages for
#' each sample.
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
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' corticalClockOut  <- corticalClock(x = beta_matrix, verbose = FALSE)
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
#'     corticalClockOut <- corticalClock(x = se_obj, verbose = FALSE)
#'   }
#' }
#' 

corticalClock <- function(x, minCoverage = 0, verbose = TRUE) {
    betaM <- .extractAssayMatrix(x)
    # --- Step 1: Load and parse coefficients (from package internal data) ---
    CorticalClockList <- loadOmniAgeRdata(
        "omniager_cortical_clock_coef",
        verbose = verbose
    )

    coefData <- CorticalClockList[["CorticalClock_coef"]]
    refData <- CorticalClockList[["CorticalClock_ref"]]
    clockName <- "corticalClock"

    # --- 2. Imputation with Ref ---
    clockProbes <- as.character(coefData[-1, 1])
    userProbes <- rownames(betaM)
    missingProbes <- setdiff(clockProbes, userProbes)

    betaComplete <- betaM

    if (length(missingProbes) > 0) {
        if (verbose) {
            message(
                "[", clockName, "] Imputing ",
                length(missingProbes),
                " missing probes using reference data."
            )
        }

        # Check whether the reference data is complete
        if (!all(missingProbes %in% names(refData))) {
            stop(
                clockName,
                "Reference data is missing for some required probes."
            )
        }

        # Construct the completion matrix:
        # Rows are missing probes and columns are samples

        refVals <- refData[missingProbes]
        refMatrix <- matrix(rep(refVals, ncol(betaM)),
            nrow = length(missingProbes),
            ncol = ncol(betaM)
        )
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
