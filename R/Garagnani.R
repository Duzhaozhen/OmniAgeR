#' @title The Garagnani ELOVL2-based Epigenetic Age Score
#'
#' @description Calculates the Garagnani epigenetic age score based on the
#' methylation level of the ELOVL2 gene.
#'
#' @details
#' This function implements the ELOVL2-based biomarker described by
#' Garagnani et al. (2012). The study identified ELOVL2 as a specific
#' hypermethylation marker that correlates strongly with chronological age
#' (Spearman's correlation coefficient = 0.92) across the entire human lifespan.
#'
#' Based on the provided coefficients (Intercept = 0, cg16867657 = 1), this
#' function currently returns the methylation beta value of the single most
#' significant CpG site located in the promoter of ELOVL2: \strong{cg16867657}.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return A numeric vector of the predicted epigenetic score.
#'
#' @export
#'
#' @references
#' Garagnani, P. et al.
#' Methylation of ELOVL2 gene as a new epigenetic marker of age.
#' \emph{Aging Cell} 2012
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
#' predAge <- garagnaniClock(x = beta_matrix)
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
#'     predAge <- garagnaniClock(x = se_obj)
#'     head(predAge)
#'   }
#' }


garagnaniClock <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_garagnani_coef", 
    clockName = "garagnaniClock",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}

