#' @title
#' Estimate epiTOC1 score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' will return the epiTOC1 score.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' The epiTOC1 score is calculated as the average beta value of 385 CpGs
#' described in Yang et al. 2014. These CpGs are promoter CpG sites that
#' localize to Polycomb group target genes that are unmethylated in 11
#' different fetal tissue types.
#'
#' @return The epiTOC1 score of each sample.
#'
#' @references
#' Yang Z, Wong A, Kuh D, et al.
#' Correlation of an epigenetic mitotic clock with cancer risk.
#' \emph{Genome Biol.} 2016
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' epiTOC1Out  <- epiTOC1(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#' library(SummarizedExperiment)
#' pheno_data <- dnamExample[[2]]
#' rownames(pheno_data) <- colnames(beta_matrix)
#'
#' se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#' epiTOC1Out <- epiTOC1(x = se_obj, verbose = FALSE)
#' }
#' @export
#'


epiTOC1 <- function(x, minCoverage = 0, verbose = TRUE) {
    betaM <- .extractAssayMatrix(x)
    epiTOC1Model <- OmniAgeRData::getOmniAgeRData(
        "omniager_epitoc1_model",
        verbose = verbose
    )
    # Obtain the list of epiTOC1 reference probes
    targetCpGs <- as.character(epiTOC1Model)
    clockWeights <- setNames(rep(1, length(targetCpGs)), targetCpGs)

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = clockWeights,
        clockName = "epiTOC1",
        minCoverage = minCoverage,
        verbose = verbose
    )

    # Calculate the score
    if (!coverageResult$pass) {
        scores <- rep(NA_real_, ncol(betaM))
        names(scores) <- colnames(betaM)
        return(scores)
    }
    scores <- colMeans(betaM[coverageResult$betaIdx, , drop = FALSE], na.rm = TRUE)

    return(scores)
}
