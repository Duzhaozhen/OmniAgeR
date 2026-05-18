#' @title The PedBE (Pediatric Buccal) Clock for DNAm Age in Children
#'
#' @description
#' Implements the Pediatric Buccal Epigenetic (PedBE) clock, specifically
#' developed to estimate DNA methylation age in
#' **pediatric (childhood) samples**, as described by McEwen et al. (2020).
#'
#' @details
#' This clock is specifically trained on and designed for pediatric buccal
#' epithelial (cheek swab) samples. The calculation is a two-step process:
#'
#' 1.  A linear predictor is first calculated from the beta values using the
#' 94-CpG elastic net coefficients . This value represents a *transformed* age.
#' 2.  A non-linear inverse age transformation (via the transformation function
#' from Horvath, 2013) is then applied. This converts the transformed age into
#' a final estimate of chronological age in years.
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted DNAm age
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#' @export
#'
#' @references
#' McEwen LM, O'Donnell KJ, McGill MG, et al.
#' The PedBE clock accurately estimates DNA methylation age in
#' pediatric buccal cells.
#' \emph{Proc Natl Acad Sci U S A.} 2020
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_pedbe_coef", verbose = FALSE)
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
#' pedBEClockOut <- pedBEClock(mockBetaM)

pedBEClock <- function(betaM,
                       minCoverage = 0,
                       verbose = TRUE) {
    pedBECoef <- loadOmniAgeRdata(
        "omniager_pedbe_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, pedBECoef, "pedBEClock",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
