#' @title Calculate the Bohlin Gestational Age (Cord Blood)
#'
#' @description
#' Implements the epigenetic clock for predicting gestational age (GA) using
#' newborn cord blood, as described by Bohlin et al. (2016).
#'
#' @param betaM a matrix of methylation beta values. 
#' Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return #' A **numeric vector** containing the predicted gestational 
#' age (in weeks) for each sample. The vector is named with the sample IDs 
#' from the `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially 
#' methylated regions. \emph{Genome Biol.} 2016
#'
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' bohlinGaOut <- bohlinGa(hannum_bmiq_m)

bohlinGa <- function(betaM,
                     minCoverage = 0.5, 
                     verbose = TRUE) {
  data("BohlinGACoef", envir = environment())
  
  predAgev <- .calLinearClock(betaM, BohlinGACoef, "bohlinGa", 
                              minCoverage, verbose)
  ## Convert the number of days into weeks
  predAgev <- predAgev/7
  return(predAgev)
  
}
