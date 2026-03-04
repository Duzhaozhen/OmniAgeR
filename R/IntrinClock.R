#' @title IntrinClock Age Prediction
#'
#' @description
#' Calculates the "IntrinClock" epigenetic age (intrinsic cellular age) using 
#' DNA methylation data(Illumina 450K and EPIC).
#' This clock is designed to be resistant to changes in immune cell composition.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes
#' and columns should be individual samples.
#'
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
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' intrinClockO <- intrinClock(hannum_bmiq_m)


intrinClock <- function(betaM,
                    minCoverage = 0.5, 
                    verbose = TRUE) {
  data("IntrinClockCoef", envir = environment())
  
  predAgev <- .calLinearClock(betaM, IntrinClockCoef, "intrinClock", 
                              minCoverage,verbose)
  
  # transformation
  predAgev <- .antiTrafo(predAgev)
  return(predAgev)
  
}

