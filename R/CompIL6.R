#' @title Calculate a DNA Methylation-Based Proxy for IL-6
#'
#' @description
#' Computes a DNAm surrogate score for Interleukin-6 (IL-6) 
#' protein levels. 
#'
#' @details
#' This function calculates the IL-6 proxy score by applying a pre-defined
#' set of coefficients to the input beta-value matrix. The coefficients
#' are loaded from the internal `IL6Coef` data object.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and 
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector containing the calculated IL-6 proxy score for each 
#' sample. The vector is named according to the column names (sample IDs) of 
#' the input matrix.
#'
#' @references
#' Stevenson AJ et al.
#' Creating and Validating a DNA Methylation-Based Proxy for Interleukin-6. 
#' \emph{J Gerontol A Biol Sci Med Sci.} 2021
#'
#' @export
#'
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' compil6Out <- compIL6(hannum_bmiq_m)

compIL6 <- function(betaM,
                    minCoverage = 0.5, 
                    verbose = TRUE) {
  # --- Step 1: Load and parse coefficients ---
  data("IL6Coef", envir = environment())
  
  il6Score <- .calLinearClock(betaM, IL6Coef, "compIL6", 
                              minCoverage,verbose)
 
  return(il6Score)
}
