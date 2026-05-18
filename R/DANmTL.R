#' @title The epigenetic age used for calculating the Leukocyte telomere length.
#'
#' @description A function to calculate the the Leukocyte telomere length (2019)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named vector of predicted Leukocyte telomere length.
#'
#' @export
#'
#' @references
#' Lu AT, Seeboth A, Tsai PC, et al.
#' DNA methylation-based estimator of telomere length
#' \emph{Aging} 2019
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_dnamtl_coef", verbose = FALSE)
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
#' dnamTlO <- dnamTL(mockBetaM, verbose = FALSE)


dnamTL <- function(betaM,
                   minCoverage = 0,
                   verbose = TRUE) {
    dnamTLCoef <- loadOmniAgeRdata(
        "omniager_dnamtl_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, dnamTLCoef, "dnamTL",
        minCoverage, verbose
    )
    return(predAgev)
}
