#' @title Calculate Universal Pan-Mammalian Blood Epigenetic Clocks
#' @description
#' Applies the two universal pan-mammalian blood clocks (Clock 2 and 3)
#' from Lu et al. (2023) to a given dataset of DNA methylation values.
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
#'
#' @return A data.frame containing the 'Sample', 'SpeciesLatinName',
#'   and the calculated ages: 'DNAmRelativeAge',
#'    'DNAmAgePanMammalianBlood2', 'DNAmRelativeAdultAge' and
#'   'DNAmAgePanMammalianBlood3'. Any additional columns
#'   from the input `sampleInfo` (like 'Age', 'Tissue') will also be returned.
#'
#'
#' @export
#'
#' @references
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#' @examples
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#'   tursiopsExample <- loadOmniAgeRdata(
#'       "omniager_tursiops_example",
#'       verbose = FALSE
#'   )
#'   
#'   # Run the calculation. Note: anageData is automatically loaded if omitted.
#'   clockResults <- panMammalianBlood(
#'       x = tursiopsExample$beta_m,
#'       speciesName = tursiopsExample$PhenoTypes$SpeciesLatinName,
#'       verbose = FALSE
#'   )
#'   print(head(clockResults))
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = tursiopsExample$beta_m),
#'       colData = tursiopsExample$PhenoTypes
#'     )
#'     
#'     se_res <- panMammalianBlood(
#'       x = se_obj,
#'       speciesName = se_obj$SpeciesLatinName,
#'       verbose = FALSE
#'     )
#'   }
#' }
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/shorvath/MammalianMethylationConsortium
# under the MIT License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------
panMammalianBlood <- function(x,
                              speciesName,
                              anageData = NULL,
                              minCoverage = 0,
                              verbose = TRUE) {
    if (verbose) message("[panMammalianBlood] Initializing calculation...")
    # --- Step 0: Universal Matrix Extraction ---
    betaM <- .extractAssayMatrix(x)
    # --- 1. Load Internal Data ---
    panMammalianBloodCoef <- loadOmniAgeRdata(
        "omniager_pan_mammalian_blood_coef",
        verbose = verbose
    )
    
    if (is.null(anageData)) {
      if (verbose) message("[panMammalianBlood] Loading default AnAge database...")
      anageData <- loadOmniAgeRdata("omniager_anage_data", verbose = FALSE)
    }
    # --- 2. Data Preparation & Merging ---
    sampleInfo <- data.frame(
        Sample = colnames(betaM),
        SpeciesLatinName = speciesName,
        stringsAsFactors = FALSE
    )

    requiredCols <- c(
      "SpeciesLatinName",
      "GestationTimeInYears",
      "averagedMaturity.yrs",
      "maxAge"
    )
    
    anageSubset <- anageData[, requiredCols, drop = FALSE]
    
    if (anyDuplicated(anageSubset$SpeciesLatinName)) {
      duplicatedSpecies <- unique(
        anageSubset$SpeciesLatinName[
          duplicated(anageSubset$SpeciesLatinName)
        ]
      )
      stop(
        "[panMammalianBlood]'anageData' contains duplicated SpeciesLatinName entries: ",
        paste(duplicatedSpecies, collapse = ", ")
      )
    }
    
    idx <- match(sampleInfo$SpeciesLatinName, anageSubset$SpeciesLatinName)
    
    info <- cbind(
      sampleInfo,
      anageSubset[
        idx,
        c("GestationTimeInYears", "averagedMaturity.yrs", "maxAge"),
        drop = FALSE
      ]
    )
    
    
    if (any(is.na(info$maxAge))) {
        missingSp <- unique(info$SpeciesLatinName[is.na(info$maxAge)])
        warning("[PanMammalianBlood] Missing AnAge data for species: ", paste(missingSp, collapse = ", "))
    }

    # The HighmaxAge required for processing Clock 2
    info$HighmaxAge <- info$maxAge * 1.3
    info$HighmaxAge[info$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus")] <-
        info$maxAge[info$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus")]

    # --- 3. Prediction (Step 1: Linear Scores) ---
    yNames <- c("Y.pred2", "Y.pred3")

    for (k in seq_along(yNames)) {
        modelDf <- panMammalianBloodCoef[[k]]

        info[[yNames[k]]] <- .calLinearClock(
            betaM = betaM[, info$Sample, drop = FALSE],
            coefData = modelDf,
            clockLabel = paste0("PanMammalian_Clock", k+1),
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # --- 4. Post-Processing (Step 2: Inverse Transformations) ---
    if (verbose) message("[PanMammalianBlood] Applying age transformations...")

    # Clock 2: Relative Age
    info$DNAmRelativeAge <- exp(-exp(-info$Y.pred2))

    info$DNAmAgePanMammalianBlood2 <- info$DNAmRelativeAge * (info$HighmaxAge + info$GestationTimeInYears) - info$GestationTimeInYears

    # Clock 3: Relative Adult Age

    a2 <- info$GestationTimeInYears / info$averagedMaturity.yrs
    m1 <- 5 * (a2^0.38)

    # Reversal conversion Clock 3
    y3 <- info$Y.pred3
    relAdultAge <- ifelse(y3 < 0, (exp(y3) - 1) * m1 + m1, y3 * m1 + m1)

    info$DNAmRelativeAdultAge <- relAdultAge
    info$DNAmAgePanMammalianBlood3 <- relAdultAge * (info$averagedMaturity.yrs + info$GestationTimeInYears) - info$GestationTimeInYears

    # --- 5. Return Output ---
    finalResults <- info[, setdiff(names(info), c("Y.pred2", "Y.pred3", "HighmaxAge"))]

    return(finalResults)
}
