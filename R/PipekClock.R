#' @title Pipek's Multi-tissue Elastic Net Epigenetic Clock (239 CpGs)
#'
#' @description
#' Implements the "elasticNet (239)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model was trained using elastic net
#' regression on a large multi-tissue dataset containing methylation
#' data from Illumina 27K, 450K, and EPIC platforms. It is designed to
#' provide improved accuracy on EPIC array data compared to the original
#' Horvath2013 clock.
#'
#' @details
#' Implements the "elasticNet (239)" model using **239 CpGs** shared
#' across Illumina 27K, 450K, and EPIC arrays.
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed prior to input.
#'   \item **Normalization:** Recommended but not strictly required.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. 
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- pipekElasticNet(x = beta_matrix, verbose = FALSE)
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
#'   predOut <- pipekElasticNet(x = se_obj, verbose = FALSE)
#' }
#'

pipekElasticNet <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_pipek_elasticnet_coef", 
    clockName = "pipekElasticNet",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}


#' @title Pipek's Filtered Horvath Epigenetic Clock (272 CpGs)
#'
#' @description
#' Implements the "filtered H (272)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model restricts feature selection to the
#' subset of CpG sites from the original Horvath pan-tissue clock that
#' are also present on the EPIC array.
#'
#' @details
#' Implements the "filtered H (272)" model. Unlike the full ElasticNet
#' model, this clock was trained by limiting candidate features to the
#' intersection of original Horvath probes and the cross-platform
#' (27K/450K/EPIC) probe set.
#'
#' **Key Features:**
#' \itemize{
#'   \item **272 CpGs:** A subset of the original Horvath clock,
#'   re-optimized for better EPIC array compatibility.
#'   \item **Best for:** Datasets pre-filtered to Horvath clock probes
#'   but requiring updated calibration.
#' }
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. 
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- pipekFilteredh(x = beta_matrix, verbose = FALSE)
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
#'   predOut <- pipekFilteredh(x = se_obj, verbose = FALSE)
#' }
#'

pipekFilteredh <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_pipek_filteredh_coef", 
    clockName = "pipekFilteredh",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}


#' @title Pipek's Retrained Horvath Epigenetic Clock (308 CpGs)
#'
#' @description
#' Implements the "retrained H (308)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model is a direct recalibration of the
#' original Horvath clock probes using a large multi-platform training set.
#'
#' @details
#' Implements the "retrained H (308)" model. It includes all **308 CpG sites** #' from the original Horvath pan-tissue clock that are present across 27K,
#' 450K, and EPIC platforms.
#'
#' **Key Features:**
#' \itemize{
#'   \item **No Feature Selection:** Coefficients were simply refitted to the
#'   new data without dropping probes.
#'   \item **Robust Update:** Serves as a direct update to the Horvath clock
#'   to correct for accuracy loss on EPIC arrays.
#' }
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. 
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- pipekRetrainedh(x = beta_matrix, verbose = FALSE)
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
#'   predOut <- pipekRetrainedh(x = se_obj, verbose = FALSE)
#' }
#'


pipekRetrainedh <- function(x, minCoverage = 0, verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_pipek_retrainedh_coef", 
    clockName = "pipekRetrainedh",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = TRUE
  )
}

