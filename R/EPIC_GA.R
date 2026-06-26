#' @title The Gestational Age (GA) clock based on 176 Illumina EPIC CpGs
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns
#' conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#' 
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' epicGaOut <- epicGa(x = beta_matrix, verbose = FALSE)
#' head(epicGaOut)
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
#'     epicGaOut <- epicGa(x = se_obj, verbose = FALSE)
#'   }
#' }



epicGa <- function(x,
                  minCoverage = 0,
                  verbose = TRUE) {
  
  predAgeDays <- .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_epic_ga_coef", 
    clockName = "epicGa",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE 
  )
  ## Convert the number of days into weeks
  predAgeWeeks <- predAgeDays / 7
  
  return(predAgeWeeks)
}

