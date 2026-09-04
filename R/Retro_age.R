#' @title Calculate the Retro-age Epigenetic Clock
#'
#' @description
#' Calculates the "Retro-age," a retroelement-based epigenetic clock for
#' chronological age, based on the models developed by Ndhlovu et al. (2024).
#' This function computes both Version 1 (V1) and Version 2 (V2) of the clock.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#'
#' @param minCoverage A numeric value between 0 and 1 (default is 0).
#' Specifies the minimum proportion of required CpGs that must be present
#' in the input matrix for the clock calculation to proceed.
#' @param verbose A logical value. If TRUE (default), the function will
#' print messages detailing the calculation steps.
#'
#' @return A list containing the predicted age for the "V1" and "V2" clocks.
#'
#' @export
#'
#' @references
#' Ndhlovu LC, Bendall ML, Dwaraka V, et al.
#' Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
#' \emph{Aging Cell.} 2024
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- retroAge(x = beta_matrix, verbose = FALSE)
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
#'   predOut <- retroAge(x = se_obj, verbose = FALSE)
#' }
#'


retroAge <- function(x,
                     minCoverage = 0,
                     verbose = TRUE) {
  
  estLv <- .runMultiEpiClockPipeline(
    x = x,
    coefName = "omniager_retroage_coef",
    clockNames = c("retroAgeV1", "retroAgeV2"),
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  return(estLv)
}

