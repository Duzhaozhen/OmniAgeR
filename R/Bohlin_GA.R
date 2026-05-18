#' @title Calculate the Bohlin Gestational Age (Cord Blood)
#'
#' @description
#' Implements the epigenetic clock for predicting gestational age (GA) using
#' newborn cord blood, as described by Bohlin et al. (2016).
#'
#' @param betaM a matrix of methylation beta values.
#' Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return #' A **numeric vector** containing the predicted gestational
#' age (in weeks) for each sample. The vector is named with the sample IDs
#' from the `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially
#' methylated regions. \emph{Genome Biol.} 2016
#'
#' @examples
#' 
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_bohlin_ga_coef", verbose = FALSE)
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
#' bohlinGaOut <- bohlinGa(mockBetaM, verbose = FALSE)

bohlinGa <- function(betaM,
                     minCoverage = 0,
                     verbose = TRUE) {
    bohlinGACoef <- loadOmniAgeRdata(
        "omniager_bohlin_ga_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, bohlinGACoef, "bohlinGa",
        minCoverage, verbose
    )
    ## Convert the number of days into weeks
    predAgev <- predAgev / 7
    return(predAgev)
}
