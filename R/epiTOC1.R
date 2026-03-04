#' @title
#' Estimate epiTOC1 score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and 
#' will return the epiTOC1 score.
#'
#' @param betaM A numeric matrix of DNAm beta values (probes as rows).
#' @param minCoverage Numeric (0-1). Minimum required probe coverage. 
#'   Default is 0.5. 
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' The epiTOC1 score is calculated as the average beta value of 385 CpGs 
#' described in Yang et al. 2014. These CpGs are promoter CpG sites that 
#' localize to Polycomb group target genes that are unmethylated in 11 
#' different fetal tissue types.
#'
#' @return The epiTOC1 score of each sample.
#'
#' @references
#' Yang Z, Wong A, Kuh D, et al.
#' Correlation of an epigenetic mitotic clock with cancer risk.
#' \emph{Genome Biol.} 2016
#' 
#' @examples
#' downloadOmniAgeRExample("LungInv")
#' loadOmniAgeRExample("LungInv")
#' eppitoc1Out <- epiTOC1(betaM = bmiq.m)
#'
#' @export
#'

epiTOC1 <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
  
  data('epiTOC1Model', envir = environment())
  # Obtain the list of epiTOC1 reference probes
  targetCpGs <- as.character(epiTOC1Model)
  clockWeights <- setNames(rep(1, length(targetCpGs)), targetCpGs)
  
  # Perform coverage check
  coverageResult <- .checkCpGCoverage(
    betaM = betaM, 
    allWeights = clockWeights, 
    clockName = "epiTOC1", 
    minCoverage = minCoverage, 
    verbose = verbose
  )
  
  # Calculate the score
  if (!coverageResult$pass) {
    scores <- rep(NA_real_, ncol(betaM))
    names(scores) <- colnames(betaM)
    return(scores)
  }
  scores <- colMeans(betaM[coverageResult$betaIdx, , drop = FALSE], na.rm = TRUE)
  
  return(scores)
}




