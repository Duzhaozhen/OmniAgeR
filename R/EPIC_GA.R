#' @title The Gestational Age (GA) clock based on 176 Illumina EPIC CpGs
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns
#' conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_epic_ga_coef", verbose = FALSE)
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
#' epicGaOut <- epicGa(mockBetaM, minCoverage = 0)

epicGa <- function(betaM,
                   minCoverage = 0,
                   verbose = TRUE) {
    epicGACoef <- loadOmniAgeRdata(
        "omniager_epic_ga_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, epicGACoef, "epicGa",
        minCoverage, verbose
    )
    ## Convert the number of days into weeks
    predAgev <- predAgev / 7
    return(predAgev)
}
