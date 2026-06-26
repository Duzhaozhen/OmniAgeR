#' @title Calculate HepatoXu ctDNA Methylation Scores
#' for Hepatocellular Carcinoma
#'
#' @description
#' This function implements the diagnostic prediction models for
#' Hepatocellular Carcinoma (HCC) based on circulating tumour DNA (ctDNA)
#' methylation markers as described by Xu et al. (2017)
#'
#' @details
#' The function calculates a composite score using a panel of HCC-specific
#' methylation markers identified through Random Forest and LASSO regression
#' analysis.
#'
#' For diagnosis, the model (cd-score) utilizes 10 genomic markers plus a
#' logistic regression intercept to differentiate HCC patients from healthy
#' controls or those with non-malignant liver diseases.
#'
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector containing the calculated methylation scores
#' (cd-score) for each sample.
#'
#'
#' @export
#'
#' @references
#' Xu, Rh., Wei, W., Krawczyk, M. et al.
#' Circulating tumour DNA methylation markers for diagnosis and prognosis of
#' hepatocellular carcinoma.
#' \emph{Nature Mater} 2017
#'
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predRes <- hepatoXuRisk(x = beta_matrix)
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
#'     predRes <- hepatoXuRisk(x = se_obj)
#'     head(predRes)
#'   }
#' }


hepatoXuRisk <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_hepato_xu_coef", 
    clockName = "hepatoXuRisk",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}

