#' @title Calculate the Mayne Placental Gestational Age.
#'
#' @description
#' Implements the 62-CpG placental epigenetic clock for estimating gestational
#' age (GA), as described by Mayne et al. (2017).
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted gestational age (in weeks) for
#'  each sample. The vector is named with the sample IDs from the
#'  `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies
#' identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predRes <- mayneGa(x = beta_matrix)
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
#'     predRes <- mayneGa(x = se_obj)
#'     head(predRes)
#'   }
#' }


mayneGa <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_mayne_ga_coef", 
    clockName = "mayneGa",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}
