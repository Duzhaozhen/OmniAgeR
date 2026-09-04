#' @title Calculate Universal Pan-Mammalian Epigenetic Clocks
#' @description
#' Applies the three universal pan-mammalian epigenetic clocks (Clock 1, 2, and 3)
#' from Lu et al. (2023) to a given dataset of DNA methylation values.
#'
#' @details
#' This function encapsulates the complete logic from the original publication's
#' script. It merges sample metadata with AnAge data, processes methylation
#' values, calculates the three clocks by matrix multiplication with pre-trained
#' models, and performs the necessary inverse mathematical transformations
#' to report age in years.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param speciesName A character string or vector specifying the Latin species
#' name(s) (e.g., "Homo sapiens").
#' @param anageData A data.frame containing the AnAge database information.
#'   Must include 'SpeciesLatinName', 'GestationTimeInYears',
#'   'averagedMaturity.yrs', and 'maxAge'.
#' @param minCoverage Numeric (0-1). Minimum required proportion of CpGs present.
#' Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#' @return A data.frame containing the 'Sample', 'SpeciesLatinName',
#'   and the calculated ages: 'DNAmAgePanMammalianClock1','DNAmRelativeAge',
#'    'DNAmAgePanMammalianClock2', 'DNAmRelativeAdultAge' and
#'   'DNAmAgePanMammalianClock3'. Any additional columns
#'   from the input `sampleInfo` (like 'Age', 'Tissue') will also be returned.
#'
#' @export
#'
#' @references
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#'
#'
#' @examples
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#'   tursiopsExample <- OmniAgeRData::getOmniAgeRData(
#'       "omniager_tursiops_example",
#'       verbose = FALSE
#'   )
#'   
#'   # Run the calculation. Note: anageData is automatically loaded if omitted.
#'   clockResults <- universalPanMammalianClocks(
#'       x = tursiopsExample$beta_m,
#'       speciesName = tursiopsExample$PhenoTypes$SpeciesLatinName,
#'       verbose = FALSE
#'   )
#'   print(head(clockResults))
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'   library(SummarizedExperiment)
#'
#'   se_obj <- SummarizedExperiment(
#'     assays = list(beta = tursiopsExample$beta_m),
#'     colData = tursiopsExample$PhenoTypes
#'   )
#'
#'   se_res <- universalPanMammalianClocks(
#'     x = se_obj,
#'     speciesName = se_obj$SpeciesLatinName,
#'     verbose = FALSE
#'   )
#' }
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/shorvath/MammalianMethylationConsortium
# under the MIT License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------
universalPanMammalianClocks <- function(x, speciesName, anageData = NULL, minCoverage = 0, verbose = TRUE) {
  .calculatePanMammalian(
    x = x, speciesName = speciesName, anageData = anageData,
    minCoverage = minCoverage, verbose = verbose,
    coefDataName = "omniager_pan_mammalian_clock_coef",
    yNames = c("Y.pred1", "Y.pred2", "Y.pred3"),
    outputPrefix = "DNAmAgePanMammalianClock",
    functionName = "UniversalPanMammalianClocks"
  )
}



#' @title Internal unified engine to calculate Pan-Mammalian epigenetic clocks
#'
#' @description
#' A shared internal unified function that encapsulates the common pipeline for 
#' universal pan-mammalian clocks (Skin, Blood, and Universal). It standardizes 
#' data extraction, AnAge database cross-referencing, linear score predictions, 
#' and inverse mathematical age transformations.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'   object containing DNA methylation beta values. Rows must represent CpG probes 
#'   and columns must represent individual samples.
#' @param speciesName A character string or vector specifying the Latin species 
#'   name(s) corresponding to each sample (e.g., \code{"Homo sapiens"}).
#' @param anageData A \code{data.frame} containing the AnAge database information. 
#'   Must include columns: \code{'SpeciesLatinName'}, \code{'GestationTimeInYears'}, 
#'   \code{'averagedMaturity.yrs'}, and \code{'maxAge'}. If \code{NULL}, the default 
#'   database is loaded internally via \code{OmniAgeRData::getOmniAgeRData}.
#' @param minCoverage Numeric value between 0 and 1. The minimum required proportion 
#'   of overlapping CpGs present for clock estimation. Default is 0.
#' @param verbose Logical. Whether to print diagnostic messages and progress status.
#' @param coefDataName Character string. The internal RData identifier of the 
#'   trained model coefficient matrix to load (e.g., \code{"omniager_pan_mammalian_skin_coef"}).
#' @param yNames Character vector. The target intermediate linear prediction columns 
#'   to compute (e.g., \code{c("Y.pred2", "Y.pred3")}).
#' @param outputPrefix Character string. The prefix string used to dynamically name 
#'   the final absolute DNAm age columns (e.g., \code{"DNAmAgePanMammalianSkin"}).
#' @param functionName Character string. The name of the parent user-facing function, 
#'   used exclusively for formatted console logging and error reporting.
#'
#' @return A \code{data.frame} combining sample identifiers, matched AnAge life-history 
#'   traits, and the calculated epigenetic ages. Depending on the active components in 
#'   \code{yNames}, the output columns include:
#' \itemize{
#'   \item{\code{Sample}: Character vector of sample column names.}
#'   \item{\code{SpeciesLatinName}: Character vector of the provided Latin species names.}
#'   \item{\code{GestationTimeInYears}: Numeric vector of gestation period in years.}
#'   \item{\code{averagedMaturity.yrs}: Numeric vector of age at sexual maturity in years.}
#'   \item{\code{maxAge}: Numeric vector of maximum lifespan in years.}
#'   \item{\code{DNAmRelativeAge}: Numeric vector of relative age (included if \code{"Y.pred2"} is active).}
#'   \item{\code{DNAmRelativeAdultAge}: Numeric vector of relative adult age (included if \code{"Y.pred3"} is active).}
#'   \item{\code{Calculated Absolute Ages}: Dynamic absolute DNAm age columns appended and named using the \code{outputPrefix} suffix (e.g., \code{Prefix2}, \code{Prefix3}).}
#' }
#'
#' @keywords internal
#' @noRd
.calculatePanMammalian <- function(x, speciesName, anageData, minCoverage, verbose, 
                                   coefDataName, yNames, outputPrefix, functionName) {
  if (verbose) message(sprintf("[%s] Initializing calculation...", functionName))
  
  # --- Step 0 & 1: Extraction & Load Internal Data ---
  betaM <- .extractAssayMatrix(x)
  modelCoefs <- OmniAgeRData::getOmniAgeRData(coefDataName, verbose = verbose)
  
  if (is.null(anageData)) {
    if (verbose) message(sprintf("[%s] Loading default AnAge database...", functionName))
    anageData <- OmniAgeRData::getOmniAgeRData("omniager_anage_data", verbose = FALSE)
  }
  
  # --- Step 2: Data Preparation & Merging ---
  sampleInfo <- data.frame(
    Sample = colnames(betaM),
    SpeciesLatinName = speciesName,
    stringsAsFactors = FALSE
  )
  
  requiredCols <- c("SpeciesLatinName", "GestationTimeInYears", "averagedMaturity.yrs", "maxAge")
  anageSubset <- anageData[, requiredCols, drop = FALSE]
  
  if (anyDuplicated(anageSubset$SpeciesLatinName)) {
    duplicatedSpecies <- unique(anageSubset$SpeciesLatinName[duplicated(anageSubset$SpeciesLatinName)])
    stop(sprintf("[%s] 'anageData' contains duplicated SpeciesLatinName entries: %s", 
                 functionName, paste(duplicatedSpecies, collapse = ", ")))
  }
  
  idx <- match(sampleInfo$SpeciesLatinName, anageSubset$SpeciesLatinName)
  info <- cbind(
    sampleInfo,
    anageSubset[idx, c("GestationTimeInYears", "averagedMaturity.yrs", "maxAge"), drop = FALSE]
  )
  
  if (any(is.na(info$maxAge))) {
    missingSp <- unique(info$SpeciesLatinName[is.na(info$maxAge)])
    warning(sprintf("[%s] Missing AnAge data for species: %s", 
                    functionName, paste(missingSp, collapse = ", ")))
  }
  
  # Process maximum lifespan limits for Clock 2
  info$HighmaxAge <- info$maxAge * 1.3
  specialSpecies <- c("Homo sapiens", "Mus musculus")
  info$HighmaxAge[info$SpeciesLatinName %in% specialSpecies] <- 
    info$maxAge[info$SpeciesLatinName %in% specialSpecies]
  
  # --- Step 3: Prediction (Linear Scores) ---
  for (k in seq_along(yNames)) {
    modelDf <- modelCoefs[[k]]
    clockLabel <- paste0(outputPrefix, if (length(yNames) == 3) k else k + 1)
    
    info[[yNames[k]]] <- .calLinearClock(
      betaM = betaM[, info$Sample, drop = FALSE],
      coefData = modelDf,
      clockLabel = clockLabel,
      minCoverage = minCoverage,
      verbose = verbose
    )
  }
  
  # --- Step 4: Post-Processing (Inverse Transformations) ---
  if (verbose) message(sprintf("[%s] Applying age transformations...", functionName))
  
  # Clock 1: Simple log-linear (Only applicable for universal multi-tissue clock)
  if ("Y.pred1" %in% yNames) {
    info[[paste0(outputPrefix, "1")]] <- exp(info$Y.pred1) - 2
  }
  
  # Clock 2: Relative Age conversion
  if ("Y.pred2" %in% yNames) {
    info$DNAmRelativeAge <- exp(-exp(-info$Y.pred2))
    info[[paste0(outputPrefix, "2")]] <- info$DNAmRelativeAge * (info$HighmaxAge + info$GestationTimeInYears) - info$GestationTimeInYears
  }
  
  # Clock 3: Relative Adult Age conversion
  if ("Y.pred3" %in% yNames) {
    a2 <- info$GestationTimeInYears / info$averagedMaturity.yrs
    m1 <- 5 * (a2^0.38)
    y3 <- info$Y.pred3
    relAdultAge <- ifelse(y3 < 0, (exp(y3) - 1) * m1 + m1, y3 * m1 + m1)
    
    info$DNAmRelativeAdultAge <- relAdultAge
    info[[paste0(outputPrefix, "3")]] <- relAdultAge * (info$averagedMaturity.yrs + info$GestationTimeInYears) - info$GestationTimeInYears
  }
  
  # --- Step 5: Clean and Return Final Output ---
  finalResults <- info[, setdiff(names(info), c(yNames, "HighmaxAge"))]
  return(finalResults)
}