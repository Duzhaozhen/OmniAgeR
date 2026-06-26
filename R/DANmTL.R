#' @title The epigenetic age used for calculating the Leukocyte telomere length.
#'
#' @description A function to calculate the the Leukocyte telomere length (2019)
#' from a DNA methylation beta value matrix.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named vector of predicted Leukocyte telomere length.
#'
#' @export
#'
#' @references
#' Lu AT, Seeboth A, Tsai PC, et al.
#' DNA methylation-based estimator of telomere length
#' \emph{Aging} 2019
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' dnamTLOut  <- dnamTL(x = beta_matrix, verbose = FALSE)
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
#'     dnamTLOut <- dnamTL(x = se_obj, verbose = FALSE)
#'   }
#' }
#' 


dnamTL <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_dnamtl_coef", 
    clockName = "dnamTL",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}

