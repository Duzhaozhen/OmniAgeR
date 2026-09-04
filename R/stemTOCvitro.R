#' @title
#' Estimate stemTOCvitro score
#'
#' @aliases stemTOCvitro
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix
#' and will return the stemTOCvitro score.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details The function will return the 0.95 upper quantile of
#' the 629 stemTOCvitro CpGs.
#' @return The stemTOCvitro score of each sample.
#'
#' @references
#' Zhu, T., Tong, H., Du, Z. et al.
#' An improved epigenetic counter to track mitotic age in normal
#' and precancerous tissues.
#' \emph{Nat Commun} 2024
#'
#' @importFrom stats quantile
#' @export
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- stemTOCvitro(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   library(SummarizedExperiment)
#'   pheno_data <- dnamExample[[2]]
#'   rownames(pheno_data) <- colnames(beta_matrix)
#'
#'   se_obj <- SummarizedExperiment(
#'     assays = list(beta = beta_matrix),
#'     colData = pheno_data
#'   )
#'
#'   predOut <- stemTOCvitro(x = se_obj, verbose = FALSE)
#' }
#' @export
stemTOCvitro <- function(x, minCoverage = 0, verbose = TRUE) {
  .calculateStemTOC(
    x = x, 
    minCoverage = minCoverage, 
    verbose = verbose,
    cpgDataName = "omniager_stemtocvitro_cpg",
    clockName = "stemTOCvitro"
  )
}