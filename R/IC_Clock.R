#' @title Calculate the Intrinsic Capacity (IC) Epigenetic Clock.
#'
#' @description
#' Implements the DNA methylation-based predictor of Intrinsic Capacity,
#' a biomarker of functional aging and mortality risk.
#'
#' @details
#' This function calculates the Intrinsic Capacity (IC) Clock, a novel
#' DNAm biomarker trained on clinical evaluations of physical and mental
#' capacities
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples. The matrix should not
#' contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @return
#' A **numeric vector** containing the predicted Intrinsic Capacity (IC) score
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#'
#' @export
#'
#' @references
#' Fuentealba M, Rouch L, Guyonnet S, et al.
#' A blood-based epigenetic clock for intrinsic capacity predicts mortality
#' and is associated with clinical, immunological and lifestyle factors.
#' \emph{Nature Aging.} 2025
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_ic_clock_coef", verbose = FALSE)
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
#' icClockO <- icClock(mockBetaM)

icClock <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    icClockCoef <- loadOmniAgeRdata(
        "omniager_ic_clock_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, icClockCoef, "icClock",
        minCoverage, verbose
    )
    return(predAgev)
}
