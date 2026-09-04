#' @title Calculate the Bohlin Gestational Age (Cord Blood)
#'
#' @description
#' Implements the epigenetic clock for predicting gestational age (GA) using
#' newborn cord blood, as described by Bohlin et al. (2016).
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return #' A **numeric vector** containing the predicted gestational
#' age (in weeks) for each sample. The vector is named with the sample IDs
#' from the `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially
#' methylated regions. \emph{Genome Biol.} 2016
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' bohlinGaOut <- bohlinGa(x = beta_matrix, verbose = FALSE)
#' head(bohlinGaOut)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'     library(SummarizedExperiment)
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     bohlinGaOut <- bohlinGa(x = se_obj, verbose = FALSE)
#' }

bohlinGa <- function(x,
                     minCoverage = 0,
                     verbose = TRUE) {
  
  predAgeDays <- .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_bohlin_ga_coef", 
    clockName = "bohlinGa",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE 
  )
  ## Convert the number of days into weeks
  predAgeWeeks <- predAgeDays / 7
  
  return(predAgeWeeks)
}
