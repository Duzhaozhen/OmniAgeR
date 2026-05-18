#' @title Calculate a DNA Methylation-Based Proxy for IL-6
#'
#' @description
#' Computes a DNAm surrogate score for Interleukin-6 (IL-6)
#' protein levels.
#'
#' @details
#' This function calculates the IL-6 proxy score by applying a pre-defined
#' set of coefficients to the input beta-value matrix.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector containing the calculated IL-6 proxy score for each
#' sample. The vector is named according to the column names (sample IDs) of
#' the input matrix.
#'
#' @references
#' Stevenson AJ et al.
#' Creating and Validating a DNA Methylation-Based Proxy for Interleukin-6.
#' \emph{J Gerontol A Biol Sci Med Sci.} 2021
#'
#' @export
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_il6_coef", verbose = FALSE)
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
#' # 4. Run the age prediction
#' compil6Out <- compIL6(mockBetaM)

compIL6 <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    # --- Step 1: Load and parse coefficients ---
    iL6Coef <- loadOmniAgeRdata(
        "omniager_il6_coef",
        verbose = verbose
    )

    il6Score <- .calLinearClock(
        betaM, iL6Coef, "compIL6",
        minCoverage, verbose
    )

    return(il6Score)
}
