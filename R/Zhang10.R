#' @title Calculate Zhang10 DNAm Age (2017)
#'
#' @description A function to calculate the Zhang10 epigenetic clock age (2017)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Zhang Y, Wilson R, Heiss J, et al.
#' DNA methylation signatures in peripheral blood strongly predict
#' all-cause mortality.
#' \emph{Nat Commun.} 2017
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_zhang10_coef", verbose = FALSE)
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
#' zhang10Out <- zhang10(mockBetaM)

zhang10 <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    zhang10Coef <- loadOmniAgeRdata(
        "omniager_zhang10_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, zhang10Coef, "zhang10",
        minCoverage, verbose
    )

    return(predAgev)
}
