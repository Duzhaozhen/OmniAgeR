#' @title Calculate the Mayne Placental Gestational Age.
#'
#' @description
#' Implements the 62-CpG placental epigenetic clock for estimating gestational
#' age (GA), as described by Mayne et al. (2017).
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted gestational age (in weeks) for
#'  each sample. The vector is named with the sample IDs from the
#'  `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies
#' identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_mayne_ga_coef", verbose = FALSE)
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
#' mayneGaOut <- mayneGa(mockBetaM)

mayneGa <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    mayneGACoef <- loadOmniAgeRdata(
        "omniager_mayne_ga_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, mayneGACoef, "mayneGa",
        minCoverage, verbose
    )
    return(predAgev)
}
