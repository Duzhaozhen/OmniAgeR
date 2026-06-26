#' @title The PedBE (Pediatric Buccal) Clock for DNAm Age in Children
#'
#' @description
#' Implements the Pediatric Buccal Epigenetic (PedBE) clock, specifically
#' developed to estimate DNA methylation age in
#' **pediatric (childhood) samples**, as described by McEwen et al. (2020).
#'
#' @details
#' This clock is specifically trained on and designed for pediatric buccal
#' epithelial (cheek swab) samples. The calculation is a two-step process:
#'
#' 1.  A linear predictor is first calculated from the beta values using the
#' 94-CpG elastic net coefficients . This value represents a *transformed* age.
#' 2.  A non-linear inverse age transformation (via the transformation function
#' from Horvath, 2013) is then applied. This converts the transformed age into
#' a final estimate of chronological age in years.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted DNAm age
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#' @export
#'
#' @references
#' McEwen LM, O'Donnell KJ, McGill MG, et al.
#' The PedBE clock accurately estimates DNA methylation age in
#' pediatric buccal cells.
#' \emph{Proc Natl Acad Sci U S A.} 2020
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' pedBEClock3Out <- pedBEClock(x = beta_matrix, verbose = FALSE)
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
#'     pedBEClock3Out <- pedBEClock(x = se_obj, verbose = FALSE)
#'   }
#' }
#' @export
#'

pedBEClock <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_pedbe_coef", 
    clockName = "pedBEClock",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}

