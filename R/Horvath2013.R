#' @title Calculate Horvath's DNAm Age (2013)
#'
#' @description A function to calculate the Horvath epigenetic clock age(2013)
#' from a DNA methylation beta value matrix.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @details Implements the Horvath (2013) pan-tissue clock. The function
#' calculates a weighted linear predictor from 353 CpGs found in the input
#' matrix and then transforms this value using a non-linear function to
#' return the final DNAm age.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Horvath S.
#' DNA methylation age of human tissues and cell types.
#' \emph{Genome Biol.} 2013
#'
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predAge <- horvath2013Clock(x = beta_matrix)
#' head(predAge)
#'
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     
#'     # Extract phenotype and ensure rownames match matrix colnames
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     # Construct the SummarizedExperiment object
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     # The function seamlessly accepts the Bioconductor object
#'     predAge <- horvath2013Clock(x = se_obj)
#'     head(predAge)
#'   }
#' }

horvath2013Clock <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_horvath2013_coef", 
    clockName = "horvath2013Clock",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}
