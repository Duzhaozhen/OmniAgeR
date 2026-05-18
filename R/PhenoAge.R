#' @title Calculate DNAm PhenoAge
#'
#' @description
#' CCalculates the DNA methylation PhenoAge, an epigenetic biomarker
#' of aging, based on a matrix of DNA methylation beta values. The
#' function implements the model originally developed by Levine et al.
#' (2018), which uses a weighted linear combination of 513 specific
#' CpG sites to predict phenotypic age.
#'
#' @details
#' This function calculates DNAm PhenoAge based on the model by
#' Levine et al. (2018), which predicts phenotypic age by
#' calculating a weighted sum of the beta values from 513 specific CpGs.
#' The function automatically loads the required model coefficients,
#' matches them with the CpGs in the input matrix, and computes the age.
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector containing the calculated PhenoAge
#' for each sample.
#'
#' @references
#' Levine ME, Lu AT, Quach A, et al.
#' An epigenetic biomarker of aging for lifespan and healthspan.
#' \emph{Aging} 2018
#
#' @export
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_phenoage_coef", verbose = FALSE)
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
#' phenoAgeOut <- phenoAge(mockBetaM)

phenoAge <- function(betaM,
                     minCoverage = 0,
                     verbose = TRUE) {
    phenoAgeCoef <- loadOmniAgeRdata(
        "omniager_phenoage_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, phenoAgeCoef, "phenoAge",
        minCoverage, verbose
    )
    return(predAgev)
}
