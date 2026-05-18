#' @title The Garagnani ELOVL2-based Epigenetic Age Score
#'
#' @description Calculates the Garagnani epigenetic age score based on the
#' methylation level of the ELOVL2 gene.
#'
#' @details
#' This function implements the ELOVL2-based biomarker described by
#' Garagnani et al. (2012). The study identified ELOVL2 as a specific
#' hypermethylation marker that correlates strongly with chronological age
#' (Spearman's correlation coefficient = 0.92) across the entire human lifespan.
#'
#' Based on the provided coefficients (Intercept = 0, cg16867657 = 1), this
#' function currently returns the methylation beta value of the single most
#' significant CpG site located in the promoter of ELOVL2: \strong{cg16867657}.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes
#' and columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return A numeric vector of the predicted epigenetic score.
#'
#' @export
#'
#' @references
#' Garagnani, P. et al.
#' Methylation of ELOVL2 gene as a new epigenetic marker of age.
#' \emph{Aging Cell} 2012
#'
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_garagnani_coef", verbose = FALSE)
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
#' garagnaniClockOut <- garagnaniClock(mockBetaM)

garagnaniClock <- function(betaM,
                           minCoverage = 0,
                           verbose = TRUE) {
    garagnaniCoef <- loadOmniAgeRdata(
        "omniager_garagnani_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, garagnaniCoef, "garagnaniClock",
        minCoverage, verbose
    )
    return(predAgev)
}
