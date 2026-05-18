#' @title Calculate Weidner Epigenetic Age (3-CpG Blood Clock)
#'
#' @description
#' Estimates biological age using the specific 3-CpG signature described by
#' Weidner et al. (2014). This model is a multivariate linear regression
#' originally designed for pyrosequencing data derived from blood samples,
#' but it can also be applied to microarray data.
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
#' A **numeric vector** containing the predicted DNAm age
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#'
#' @export
#'
#' @references
#' Weidner, C.I., Lin, Q., Koch, C.M. et al.
#' Aging of blood can be tracked by DNA methylation changes
#' at just three CpG sites.
#' \emph{Genome Biol} 2014
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_weidner_coef", verbose = FALSE)
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
#' weidnerClockOut <- weidnerClock(mockBetaM)

weidnerClock <- function(betaM,
                         minCoverage = 0,
                         verbose = TRUE) {
    weidnerCoef <- loadOmniAgeRdata(
        "omniager_weidner_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, weidnerCoef, "weidnerClock",
        minCoverage, verbose
    )
    return(predAgev)
}
