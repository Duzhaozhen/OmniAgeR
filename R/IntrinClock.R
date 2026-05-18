#' @title IntrinClock Age Prediction
#'
#' @description
#' Calculates the "IntrinClock" epigenetic age (intrinsic cellular age) using
#' DNA methylation data(Illumina 450K and EPIC).
#' This clock is designed to be resistant to changes in immune cell composition.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes
#' and columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @details
#' The IntrinClock utilizes an elastic net model trained on 410 CpGs and
#' refined to 381 active predictors. It predicts age by calculating a linear
#' combination of CpG beta values and coefficients, followed by an inverse
#' Horvath transformation to convert the linear predictor to years.
#'
#' The model coefficients are stored in \code{IntrinClockCoef}, which contains:
#' \itemize{
#'   \item Intercept term
#'   \item 381 non-zero CpG coefficients (sparsity optimized)
#' }
#'
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding
#' to the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Tomusiak, A., Floro, A., Tiwari, R. et al.
#' Development of an epigenetic clock resistant to changes in
#' immune cell composition.
#' \emph{Commun Biol} 2024
#'
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_intrin_clock_coef", verbose = FALSE)
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
#' intrinClockO <- intrinClock(mockBetaM)

intrinClock <- function(betaM,
                        minCoverage = 0,
                        verbose = TRUE) {
    intrinClockCoef <- loadOmniAgeRdata(
        "omniager_intrin_clock_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, intrinClockCoef, "intrinClock",
        minCoverage, verbose
    )

    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
