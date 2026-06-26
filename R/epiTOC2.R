#' @title
#' Estimate epiTOC2 scores
#'
#' @aliases epiTOC2
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an
#' age vector (optional) and will return the epiTOC2 scores.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param age Optional numeric vector representing chronological ages
#'   of the samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' Building upon a dynamic model of DNA methylation gain in unmethylated
#'  CpG-rich regions, epiTOC2 can directly estimate the cumulative number
#'  of stem cell divisions in a tissue. The details of the algorithm are
#'  described in Teschendorff et al. 2020.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using the full epiTOC2 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using an approximation of epiTOC2
#'   which assumes all epiTOC2 CpGs have beta-values exactly 0 in the
#'   fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the
#'   estimated average lifetime intrinsic rate of stem-cell division per
#'   sample, as derived from epiTOC2
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate
#'   for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#'
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @importFrom stats median
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' epiTOC2Out <- epiTOC2(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     epiTOC2Out <- epiTOC2(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
#'


epiTOC2 <- function(x, age = NULL, minCoverage = 0, verbose = TRUE) {
    betaM <- .extractAssayMatrix(x)
    estParams <- loadOmniAgeRdata(
        "omniager_epitoc2_model",
        verbose = verbose
    )

    dummyWeights <- setNames(seq_len(nrow(estParams)), rownames(estParams))

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = dummyWeights,
        clockName = "epiTOC2",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
        return(list(
            tnsc = rep(NA_real_, ncol(betaM)), tnsc2 = rep(NA_real_, ncol(betaM)),
            irS = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
            irS2 = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
            irT = NA_real_, irT2 = NA_real_
        ))
    }

    # Extract the matching data and parameters
    matchedParams <- estParams[names(coverageResult$weightsSubset), , drop = FALSE]
    subBeta <- betaM[coverageResult$betaIdx, , drop = FALSE]

    deltaV <- matchedParams[, 1]
    beta0V <- matchedParams[, 2]

    # Core algorithm implementation
    # Full Model
    resScale <- 2 / (deltaV * (1 - beta0V))
    tnscV <- colMeans(sweep(subBeta, 1, beta0V, "-") * resScale, na.rm = TRUE)

    # Approximation (Assuming beta0 = 0)
    resScale2 <- 2 / deltaV
    tnsc2V <- colMeans(subBeta * resScale2, na.rm = TRUE)

    # Intrinsic Rate
    irS <- NULL
    irS2 <- NULL
    irT <- NULL
    irT2 <- NULL

    if (!is.null(age)) {
        irS <- tnscV / age
        irS2 <- tnsc2V / age
        irT <- median(irS, na.rm = TRUE)
        irT2 <- median(irS2, na.rm = TRUE)
    }

    return(list(
        tnsc = tnscV,
        tnsc2 = tnsc2V,
        irS = irS,
        irS2 = irS2,
        irT = irT,
        irT2 = irT2
    ))
}
