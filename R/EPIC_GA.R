#' @title The Gestational Age (GA) clock based on 176 Illumina EPIC CpGs
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).  
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' 
#' @return A named vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns 
#' conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' epicGaOut <- epicGA(hannum_bmiq_m)
#' 

epicGa <- function(betaM,
                   minCoverage = 0.5, 
                   verbose = TRUE) {
  data("EPICGACoef", envir = environment())
  
  predAgev <- .calLinearClock(betaM, EPICGACoef, "epicGa", 
                              minCoverage, verbose)
  ## Convert the number of days into weeks
  predAgev <- predAgev/7
  return(predAgev)
  
}

