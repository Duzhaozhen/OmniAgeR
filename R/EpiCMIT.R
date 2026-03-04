#' @title Estimate EpiCMIT scores
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will 
#' return EpiCMIT scores.
#' 
#' @details
#' EpiCMIT is based on two groups of CpGs, a 184 age-associated hypermethylated 
#' CpG list and a 1164 hypomethylated CpG list (Duran-Ferrer et al. 2020) . 
#' The EpiCMIT function calculates the average DNAm of the two groups separately 
#' and returns the EpiCMIT_hyper and EpiCMIT_hypo scores.
#' 
#' @param betaM A numeric matrix of DNAm beta values (probes as rows).
#' @param minCoverage Numeric (0-1). Minimum required probe coverage. 
#'   Default is 0.5. 
#' @param verbose Logical. Whether to print coverage statistics.
#' 
#' 
#' @return A list containing the following entries.
#'
#' * hyperScore: EpiCMIT-hyper score of each sample.
#'
#' * hypoScore: EpiCMIT-hypo score of each sample.
#'
#' @references
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors 
#' and predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#'
#' @examples
#' downloadOmniAgeRExample("LungInv")
#' loadOmniAgeRExample("LungInv")
#' epicmitOut <- epiCMIT(betaM = bmiq.m)
#' @export
#'

epiCMIT <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
  data('EpiCMITModel', envir = environment())
  
  # 1. Prepare the list of probes
  hyperProbes <- epiCMIT.df[grep("hyper", epiCMIT.df$epiCMIT.class), 1]
  hypoProbes  <- epiCMIT.df[grep("hypo", epiCMIT.df$epiCMIT.class), 1]
  
  hyperWeights <- setNames(rep(1, length(hyperProbes)), hyperProbes)
  hypoWeights  <- setNames(rep(1, length(hypoProbes)), hypoProbes)
  
  # 2. Call the helper function
  resHyper <- .checkCpGCoverage(betaM, hyperWeights, "EpiCMIT-Hyper", 
                                minCoverage, verbose)
  resHypo  <- .checkCpGCoverage(betaM, hypoWeights, "EpiCMIT-Hypo", 
                                minCoverage, verbose)
  
  # 3. Predcition
  hyperScore <- if (resHyper$pass) {
    colMeans(betaM[resHyper$betaIdx, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(betaM))
  }
  
  hypoScore <- if (resHypo$pass) {
    1 - colMeans(betaM[resHypo$betaIdx, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(betaM))
  }
  
  return(list(
    hyperScore = hyperScore,
    hypoScore  = hypoScore
  ))
}

