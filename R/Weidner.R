#' @title Calculate Weidner Epigenetic Age (3-CpG Blood Clock)
#'
#' @description
#' Estimates biological age using the specific 3-CpG signature described by
#' Weidner et al. (2014). This model is a multivariate linear regression
#' originally designed for pyrosequencing data derived from blood samples,
#' but it can also be applied to microarray data.
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
#' for each sample. 
#' @export
#'
#' @references
#' Weidner, C.I., Lin, Q., Koch, C.M. et al.
#' Aging of blood can be tracked by DNA methylation changes
#' at just three CpG sites.
#' \emph{Genome Biol} 2014
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- weidnerClock(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   library(SummarizedExperiment)
#'   pheno_data <- dnamExample[[2]]
#'   rownames(pheno_data) <- colnames(beta_matrix)
#'
#'   se_obj <- SummarizedExperiment(
#'     assays = list(beta = beta_matrix),
#'     colData = pheno_data
#'   )
#'
#'   predOut <- weidnerClock(x = se_obj, verbose = FALSE)
#' }
#'

weidnerClock <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_weidner_coef", 
    clockName = "weidnerClock",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}







