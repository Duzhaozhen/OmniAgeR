#' @title Calculate the Lee gestational age
#'
#' @description
#' Implements the placental epigenetic clocks for estimating gestational
#' age (GA) using DNA methylation data, as described by Lee et al. (2019).
#'
#' @details
#' This function computes three distinct GA clocks derived from the models
#' presented in the Lee et al. (2019) study.
#'
#' The implemented clocks are:
#' \itemize{
#'   \item \strong{`LeeControl`}: The control model (546 CpGs).
#'   \item \strong{`LeeRobust`}: The robust model (558 CpGs).
#'   \item \strong{`LeeRefinedRobust`}: The refined robust model (395 CpGs).
#' }
#' The function iterates through each clock, matches the required
#' CpGs (e.g., 546 for `LeeControl`) with the columns in the input matrix,
#' and calculates a linear prediction of GA.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A `list` containing three named numeric vectors.
#' Each vector represents the predicted gestational age (in weeks)
#' for the corresponding samples.
#' \itemize{
#'   \item \strong{`LeeControl`}: Numeric vector of predicted GAs
#'   from the Control model.
#'   \item \strong{`LeeRobust`}: Numeric vector of predicted GAs
#'   from the Robust model.
#'   \item \strong{`LeeRefinedRobust`}: Numeric vector of predicted GAs
#'   from the Refined Robust model.
#' }
#' Each vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#' @export
#'
#' @references
#' Lee Y, Choufani S, Weksberg R, et al.
#' Placental epigenetic clocks: estimating gestational age using placental
#' DNA methylation levels.
#' \emph{Aging} 2019
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
#' predRes <- LeeGa(x = beta_matrix)
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
#'     predRes <- LeeGa(x = se_obj)
#'     head(predRes)
#'   }
#' }



LeeGa <- function(x,
                  minCoverage = 0,
                  verbose = TRUE) {
  
  # Delegate the entire calculation to the unified multi-model pipeline
  estLv <- .runMultiEpiClockPipeline(
    x = x,
    coefName = "omniager_lee_ga_coef",
    clockNames = c("LeeControl", "LeeRobust", "LeeRefinedRobust"),
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  return(estLv)
}

