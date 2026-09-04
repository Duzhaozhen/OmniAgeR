#' @title Calculate Causal, Damage, and Adaptation Epigenetic Clocks
#'
#' @description
#' Calculates three related epigenetic clocks (Causal, Damage, Adaptation)
#' from a DNA methylation beta value matrix. These clocks were developed to
#' distinguish different aspects of aging.
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' This function implements the three epigenetic clocks described by
#' Ying et al. (2024) to dissect biological aging into distinct components.
#' The Damage clock is designed to capture the accumulation of molecular damage.
#' The Adaptation clock reflects the body's adaptive responses to this damage.
#' The Causal clock is enriched for CpGs with a causal effect on mortality.
#' Each score is calculated as a weighted linear sum of beta values from
#' its specific set of CpG sites.
#'
#' @return A list containing the predicted scores for the "Causal", "Damage",
#' and "Adaptation" clocks.
#'
#' @references
#' Ying K, Liu H, Tarkhov AE, et al.
#' Causality-enriched epigenetic age uncouples damage and adaptation.
#' \emph{Nat Aging} 2024
#'
#' @export
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' causalClockO <- causalClock(x = beta_matrix, verbose = FALSE)
#' head(causalClockO$CausalAge)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'     library(SummarizedExperiment)
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     causal_se_res <- causalClock(x = se_obj, verbose = FALSE)
#' }
causalClock <- function(x,
                        minCoverage = 0,
                        verbose = TRUE) {
  
  estLv <- .runMultiEpiClockPipeline(
    x = x,
    coefName = "omniager_causal_clocks_coef",
    clockNames = c("CausalAge", "DamAge", "AdaptAge"),
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  
  return(estLv)
}

