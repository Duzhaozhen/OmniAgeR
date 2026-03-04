#' @title The epigenetic age used for calculating the Leukocyte telomere length.
#'
#' @description A function to calculate the the Leukocyte telomere length (2019) 
#' from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and 
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' 
#' @return A named vector of predicted Leukocyte telomere length.
#'
#' @export
#'
#' @references
#' Lu AT, Seeboth A, Tsai PC, et al.
#' DNA methylation-based estimator of telomere length
#' \emph{Aging} 2019
#'
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' dnamTlO <- dnamTL(hannum_bmiq_m)


dnamTL <- function(betaM,
                   minCoverage = 0.5, 
                   verbose = TRUE) {
  data("DNAmTLCoef", envir = environment()) 
  predAgev <- .calLinearClock(betaM, DNAmTLCoef, "dnamTL", 
                              minCoverage, verbose)
  return(predAgev)
  
}

