#' @title Calculate Centenarian Epigenetic Clocks (Eric Dec et al.)
#'
#' @description  Calculates the Centenarian epigenetic clocks
#' (ENCen40 and ENCen100) developed by Eric Dec et al. (2023). This function
#' serves as a wrapper that loads the internal clock coefficients and computes
#' the linear predictors for each clock using the helper function.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A list containing the predicted scores for each Centenarian clock.
#' \itemize{
#'   \item `ENCen40`: Elastic net clock trained on individuals aged 40+.
#'   \item `ENCen100`: Elastic net clock specifically trained on centenarians (100+).
#' }
#'
#' @export
#'
#' @references
#' Dec, E., Clement, J., Cheng, K. et al.
#' Centenarian clocks: epigenetic clocks for validating claims of
#' exceptional longevity. \emph{GeroScience} 2023
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' centenarian_res <- centenarianClock(x = beta_matrix, verbose = FALSE)
#' head(centenarian_res$ENCen100)
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
#'     cent_se_res <- centenarianClock(x = se_obj, verbose = FALSE)
#'   }
#' }
centenarianClock <- function(x,
                             minCoverage = 0,
                             verbose = TRUE) {
  estLv <- .runMultiEpiClockPipeline(
    x = x,
    coefName = "omniager_centenarian_coef",
    clockNames = c("ENCen40", "ENCen100"),
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  return(estLv)
}
