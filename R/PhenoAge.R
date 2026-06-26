#' @title Calculate DNAm PhenoAge
#'
#' @description
#' CCalculates the DNA methylation PhenoAge, an epigenetic biomarker
#' of aging, based on a matrix of DNA methylation beta values. The
#' function implements the model originally developed by Levine et al.
#' (2018), which uses a weighted linear combination of 513 specific
#' CpG sites to predict phenotypic age.
#'
#' @details
#' This function calculates DNAm PhenoAge based on the model by
#' Levine et al. (2018), which predicts phenotypic age by
#' calculating a weighted sum of the beta values from 513 specific CpGs.
#' The function automatically loads the required model coefficients,
#' matches them with the CpGs in the input matrix, and computes the age.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector containing the calculated PhenoAge
#' for each sample.
#'
#' @references
#' Levine ME, Lu AT, Quach A, et al.
#' An epigenetic biomarker of aging for lifespan and healthspan.
#' \emph{Aging} 2018
#
#' @export
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- phenoAge(x = beta_matrix, verbose = FALSE)
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
#'     predOut <- phenoAge(x = se_obj, verbose = FALSE)
#'   }
#' }
#'


phenoAge <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_phenoage_coef", 
    clockName = "phenoAge",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}

