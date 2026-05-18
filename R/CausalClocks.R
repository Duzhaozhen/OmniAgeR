#' @title Calculate Causal, Damage, and Adaptation Epigenetic Clocks
#'
#' @description
#' Calculates three related epigenetic clocks (Causal, Damage, Adaptation)
#' from a DNA methylation beta value matrix. These clocks were developed to
#' distinguish different aspects of aging.
#' @param betaM DNAm beta value matrix with rows labeling Illumina 450k/EPIC
#' CpGs and columns labeling samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' This function implements the three epigenetic clocks described by
#' Ying et al. (2024) to dissect biological aging into distinct components.
#' The Damage clock is designed to capture the accumulation of molecular damage.
#' The Adaptation clock reflects the body's adaptive responses to this damage.
#' The Causal clock is enriched for CpGs with a causal effect on mortality.
#' Each score is calculated as a weighted linear sum of beta values from
#' its specific set of CpG sites.
#'
#' @return A list containing the predicted scores for the "Causal", "Damage",
#' and "Adaptation" clocks.
#'
#' @references
#' Ying K, Liu H, Tarkhov AE, et al.
#' Causality-enriched epigenetic age uncouples damage and adaptation.
#' \emph{Nat Aging} 2024
#'
#' @export
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_causal_clocks_coef", verbose = FALSE)
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
#' causalClockO <- causalClock(mockBetaM, verbose = FALSE)

causalClock <- function(betaM,
                        minCoverage = 0,
                        verbose = TRUE) {
    # --- Step 1: Load Coefficients ---
    causalClockList <- loadOmniAgeRdata(
        "omniager_causal_clocks_coef",
        verbose = verbose
    )
    # Define the specific names for the three sub-clocks
    clockNames <- c("CausalAge", "DamAge", "AdaptAge")

    # --- Step 2: Calculate Scores for Each Clock ---
    estLv <- list()

    # Loop through the list of coefficients (causalClockList)
    for (i in seq_along(causalClockList)) {
        # Call the internal helper to handle all calculation and logging
        estLv[[i]] <- .calLinearClock(
            betaM = betaM,
            coefData = causalClockList[[i]],
            clockLabel = clockNames[i],
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # Assign names to the result list
    names(estLv) <- clockNames

    return(estLv)
}
