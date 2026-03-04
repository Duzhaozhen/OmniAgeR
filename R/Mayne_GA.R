#' @title Calculate the Mayne Placental Gestational Age.
#'
#' @description
#' Implements the 62-CpG placental epigenetic clock for estimating gestational 
#' age (GA), as described by Mayne et al. (2017).
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).  
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of 
#'   required CpGs that must be present. Default is 0.5.
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
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' Mayne_GA.out <- Mayne_GA(hannum_bmiq_m)


mayneGa <- function(betaM,
                     minCoverage = 0.5, 
                     verbose = TRUE) {
  data("MayneGACoef", envir = environment())
  
  predAgev <- .calLinearClock(betaM, MayneGACoef, "mayneGa", 
                              minCoverage, verbose)
  return(predAgev)
  
}
