#' @title Wu's Epigenetic Clock for Pediatric Age Estimation
#'
#' @description
#' Estimates biological age (in years) based on DNA methylation data using the method
#' proposed by Wu et al. (2019). This clock is specifically designed for pediatric
#' cohorts and utilizes a non-linear transformation to capture rapid developmental
#' changes in early life.
#'
#' @details
#' The calculation involves two main steps:
#' \enumerate{
#'   \item Calculation of a linear predictor using weighted CpG beta values.
#'   \item Transformation of the linear predictor into biological age (years) using
#'   a specific "anti-transformation" function with a toddler age offset of 48 months.
#' }
#'
#' \strong{Data Requirements:}
#' The input \code{beta.m} must be a matrix of Beta values (0 to 1). The function
#' expects CpG probes as row names.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'   object containing DNA methylation beta values. Rows should be CpG probes and 
#'   columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. 
#'
#' @export
#'
#' @references
#' Wu, Xiaohui et al.
#' DNA methylation profile is a quantitative measure of biological aging in children
#' \emph{Aging} 2019
#' 
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- wuClock(x = beta_matrix, verbose = FALSE)
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
#'   predOut <- wuClock(x = se_obj, verbose = FALSE)
#' }
#'

wuClock <- function(x,
                    minCoverage = 0,
                    verbose = TRUE) {

  predAgev <- .runEpiClockPipeline(
    x = x,
    coefName = "omniager_wu_clock_coef",
    clockName = "wuClock",
    minCoverage = minCoverage,
    verbose = verbose,
    useHorvathTrafo = FALSE 
  )
  
  predAgev <- .antiTrafo(predAgev, 48)
  
  ## transform to years
  predAgev <- predAgev / 12
  
  return(predAgev)
}