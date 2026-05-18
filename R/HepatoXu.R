#' @title Calculate HepatoXu ctDNA Methylation Scores
#' for Hepatocellular Carcinoma
#'
#' @description
#' This function implements the diagnostic prediction models for
#' Hepatocellular Carcinoma (HCC) based on circulating tumour DNA (ctDNA)
#' methylation markers as described by Xu et al. (2017)
#'
#' @details
#' The function calculates a composite score using a panel of HCC-specific
#' methylation markers identified through Random Forest and LASSO regression
#' analysis.
#'
#' For diagnosis, the model (cd-score) utilizes 10 genomic markers plus a
#' logistic regression intercept to differentiate HCC patients from healthy
#' controls or those with non-malignant liver diseases.
#'
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector containing the calculated methylation scores
#' (cd-score) for each sample.
#'
#'
#' @export
#'
#' @references
#' Xu, Rh., Wei, W., Krawczyk, M. et al.
#' Circulating tumour DNA methylation markers for diagnosis and prognosis of
#' hepatocellular carcinoma.
#' \emph{Nature Mater} 2017
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_hepato_xu_coef", verbose = FALSE)
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
#' hepatoXuRiskO <- hepatoXuRisk(mockBetaM)

hepatoXuRisk <- function(betaM,
                         minCoverage = 0,
                         verbose = TRUE) {
    hepatoXuCoef <- loadOmniAgeRdata(
        "omniager_hepato_xu_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, hepatoXuCoef, "hepatoXuRisk",
        minCoverage, verbose
    )
    return(predAgev)
}
