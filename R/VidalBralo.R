#' @title Calculate the Vidal-Bralo Epigenetic Clock (8-CpG Model)
#'
#' @description
#' Implements the simplified 8-CpG epigenetic clock for estimating
#' chronological age in adults, as described by Vidal-Bralo et al. (2016).
#'
#' @details
#' This function calculates the highly simplified epigenetic clock developed
#' for use in adult whole blood samples. The model is a linear predictor
#' based on 8 specific CpG sites.
#'
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted chronological age (in years)
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#'
#' @export
#'
#' @references
#' Vidal-Bralo L, Lopez-Golan Y, Gonzalez A.
#' Simplified Assay for Epigenetic Age Estimation in Whole Blood of Adults.
#' \emph{Front Genet.} 2016
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_vidalbralo_coef", verbose = FALSE)
#' 
#' # 2. Extract feature names and exclude potential intercept terms
#' allFeatures <- unique(unlist(modelCoef, use.names = FALSE))
#' requiredCpGs <- allFeatures[grep("^cg", allFeatures)]
#' if (length(requiredCpGs) == 0) requiredCpGs <- allFeatures 
#' 
#' # 3. Generate a mock micro-beta matrix for 2 samples in memory
#' mockBetaM <- matrix(
#'     runif(length(requiredCpGs) * 2, min = 0, max = 1),
#'     nrow = length(requiredCpGs),
#'     dimnames = list(requiredCpGs, c("Sample1", "Sample2"))
#' )
#' 
#' # 4. Run the age prediction
#' vidalBraloClockOut <- vidalBraloClock(mockBetaM)

vidalBraloClock <- function(betaM,
                            minCoverage = 0,
                            verbose = TRUE) {
    vidalBraloCoef <- loadOmniAgeRdata(
        "omniager_vidalbralo_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, vidalBraloCoef, "vidalBraloClock",
        minCoverage, verbose
    )
    return(predAgev)
}
