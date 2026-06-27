#' @title Estimate epiTOC3 scores
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' an age vector (optional) and will return the epiTOC3 scores.
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
#' Building upon a dynamic model of DNA methylation gain in 170 unmethylated
#' population doubling associated CpGs, epiTOC3 can directly estimate the
#' cumulative number of stem cell divisions in a tissue.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell
#'   per year and per sample using the full epiTOC3 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell
#'   per year and per sample using an approximation of epiTOC3 which assumes all
#'   epiTOC3 CpGs have beta-values exactly 0 in the fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the
#'   estimated average lifetime intrinsic rate of stem-cell division per
#'   sample, as derived from epiTOC3
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate
#'   for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#' * avETOC3: The simple average over the 170 epiTOC3 sites.
#'
#'
#' @importFrom stats median
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' epiTOC3Out <- epiTOC3(x = beta_matrix, verbose = FALSE)
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
#'     epiTOC3Out <- epiTOC3(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
#'
epiTOC3 <- function(x, age = NULL, minCoverage = 0, verbose = TRUE) {
  .calculateEpiTOC(
    x = x, 
    age = age, 
    minCoverage = minCoverage, 
    verbose = verbose,
    modelName = "omniager_epitoc3_model",
    clockName = "epiTOC3",
    calcAvETOC3 = TRUE
  )
}
