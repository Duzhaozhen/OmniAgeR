#' @title IntrinClock Age Prediction
#'
#' @description
#' Calculates the "IntrinClock" epigenetic age (intrinsic cellular age) using
#' DNA methylation data(Illumina 450K and EPIC).
#' This clock is designed to be resistant to changes in immune cell composition.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @details
#' The IntrinClock utilizes an elastic net model trained on 410 CpGs and
#' refined to 381 active predictors. It predicts age by calculating a linear
#' combination of CpG beta values and coefficients, followed by an inverse
#' Horvath transformation to convert the linear predictor to years.
#'
#' The model coefficients are stored in \code{IntrinClockCoef}, which contains:
#' \itemize{
#'   \item Intercept term
#'   \item 381 non-zero CpG coefficients (sparsity optimized)
#' }
#'
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding
#' to the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Tomusiak, A., Floro, A., Tiwari, R. et al.
#' Development of an epigenetic clock resistant to changes in
#' immune cell composition.
#' \emph{Commun Biol} 2024
#'
#'
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predRes <- icClock(x = beta_matrix)
#' head(predRes)
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
#'     predRes <- icClock(x = se_obj)
#'     head(predRes)
#'   }
#' }


intrinClock <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_intrin_clock_coef", 
    clockName = "intrinClock",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}


