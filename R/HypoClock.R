#' @title Estimate hypoClock score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will
#' return the HypoClock score.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' The hypoClock score is defined by Teschendorff (2020) and is calculated
#' as 1 minus the average beta value of 678 solo-WCGW sites.
#'
#' @return The HypoClock score of each sample.
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' hypoClockOut  <- hypoClock(x = beta_matrix, verbose = FALSE)
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
#'     hypoClockOut <- hypoClock(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
#'

hypoClock <- function(x, minCoverage = 0, verbose = TRUE) {
    betaM <- .extractAssayMatrix(x)
    hypoClockCpG <- loadOmniAgeRdata(
        "omniager_hypoclock_cpg",
        verbose = verbose
    )

    # Prepare reference sites
    soloCpGs <- as.character(hypoClockCpG)
    clockWeights <- setNames(rep(1, length(soloCpGs)), soloCpGs)

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = clockWeights,
        clockName = "hypoClock",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
        scores <- rep(NA_real_, ncol(betaM))
        names(scores) <- colnames(betaM)
        return(scores)
    }

    # 5. Calculate score(1 - mean beta)
    scores <- 1 - colMeans(betaM[coverageResult$betaIdx, , drop = FALSE],
        na.rm = TRUE
    )

    return(scores)
}
