#' @title Calculate GrimAge1
#'
#' @description
#' Calculates DNA methylation GrimAge1, a composite biomarker of mortality risk
#' and biological aging.
#'
#' @details
#' This function calculates DNAm GrimAge1 in a multi-step process. First, it
#' predicts DNAm-based surrogate biomarkers for several plasma proteins from
#' the input beta values. These predicted biomarkers, along with chronological
#' age and sex, are then used to calculate a composite mortality risk score.
#' This score is calibrated to the scale of chronological age to produce the
#' final `DNAmGrimAge1`.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param age A numeric vector of chronological ages for the samples
#'  corresponding to the columns in `betaM`.
#' @param sex A character vector of sample sexes. Must contain "Male" or
#'  "Female" for each sample.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A data.frame containing the following columns:
#' \itemize{
#'   \item `SampleID`: Identifier for each sample.
#'   \item `DNAm...`: Columns for each of the predicted surrogate biomarkers
#'   (e.g., `DNAmADM`, `DNAmGDF15`).
#'   \item `DNAmGrimAge1`: The final calibrated GrimAge1 score.
#' }
#'
#'
#'
#' @references
#' Lu AT, Quach A, Wilson JG, et al.
#' DNA methylation GrimAge strongly predicts lifespan and healthspan
#' \emph{Aging} 2019
#'
#' @export
#'
#' @examples
#' # ====================================================================
#' # Example 1: Direct Matrix Input 
#' # ====================================================================
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' phenoTypes_df <- dnamExample[[2]]
#' 
#' age <- phenoTypes_df$Age
#' sex <- ifelse(phenoTypes_df$Sex == "F", "Female", "Male")
#' 
#' # Calculate GrimAge first
#' GrimAge1O <- grimAge1(x = beta_matrix, age = age, sex = sex, verbose = FALSE)
#' 
#' \dontrun{
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     GrimAge1O <- grimAge1(
#'       x = se_obj,
#'       age = se_obj$Age,
#'       sex = ifelse(se_obj$Sex == "F", "Female", "Male"),
#'       verbose = FALSE
#'     )
#'   }
#' }

grimAge1 <- function(x, age, sex, minCoverage = 0, verbose = TRUE) {
  .calculateGrimAge(
    x = x, 
    age = age, 
    sex = sex, 
    minCoverage = minCoverage, 
    verbose = verbose,
    modelName = "omniager_grimage1_model",
    clockName = "GrimAge1",
    outputColName = "DNAmGrimAge1",
    applyRenameMap = FALSE
  )
}




#' @title Internal helper function to calculate GrimAge variants
#'
#' @description
#' A shared internal engine designed to compute the GrimAge1 and GrimAge2 
#' biological aging clocks. It executes a three-phase pipeline: predicting surrogate 
#' biomarkers, calculating a composite COX mortality risk score, and calibrating 
#' the score to chronological age.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'   object containing DNA methylation beta values.
#' @param age A numeric vector of chronological ages for the samples.
#' @param sex A character vector of sample sexes (\code{"Male"} or \code{"Female"}).
#' @param minCoverage Numeric value between 0 and 1. The minimum required proportion 
#'   of required CpGs.
#' @param verbose Logical. Whether to print status messages.
#' @param modelName Character string. The internal RData identifier for the GrimAge 
#'   model to load (e.g., \code{"omniager_grimage1_model"}).
#' @param clockName Character string. Used for coverage check reporting 
#'   (e.g., \code{"GrimAge1"}).
#' @param outputColName Character string. The specific name for the final age 
#'   column (e.g., \code{"DNAmGrimAge1"}).
#' @param applyRenameMap Logical. If \code{TRUE}, applies a specific column renaming 
#'   map for surrogate biomarkers (used specifically for GrimAge 2).
#'
#' @return A \code{data.frame} containing the \code{SampleID}, predicted surrogate 
#'   biomarkers, and the final calibrated GrimAge score.
#'
#' @importFrom stats setNames
#' @keywords internal
#' @noRd
.calculateGrimAge <- function(x, age, sex, minCoverage, verbose, 
                              modelName, clockName, outputColName, applyRenameMap = FALSE) {
  # --- Step 0: Extraction & Load Data ---
  betaM <- .extractAssayMatrix(x)
  grimageModel <- loadOmniAgeRdata(modelName, verbose = verbose)
  
  protCoefs <- grimageModel[[1]]    # Phase 1 weights
  finalModel <- grimageModel[[2]]   # Phase 2 weights (COX)
  calibParams <- grimageModel[[3]]  # Phase 3 parameters
  
  uqCpgs <- unique(protCoefs$var[startsWith(protCoefs$var, "cg")])
  fakeWeights <- setNames(rep(1, length(uqCpgs)), uqCpgs)
  
  # Coverage check
  coverage <- .checkCpGCoverage(betaM, fakeWeights, clockName, minCoverage, verbose)
  
  # --- Step 1: Handle Covariates ---
  femaleVec <- ifelse(sex == "Female", 1, 0)
  availableCpGs <- intersect(protCoefs$var, rownames(betaM))
  
  if (verbose) {
    nTotalCpGs <- length(setdiff(unique(protCoefs$var), c("Intercept", "Age")))
    message(sprintf(
      "[%s] Found %d / %d required CpGs (%.1f%%).",
      clockName, length(availableCpGs), nTotalCpGs, (length(availableCpGs) / nTotalCpGs) * 100
    ))
  }
  
  # --- Phase 1: Predict Surrogate Biomarkers ---
  proteinNames <- unique(protCoefs$Y.pred)
  protPredList <- list()
  
  for (pName in proteinNames) {
    pSub <- protCoefs[protCoefs$Y.pred == pName, ]
    presentVars <- intersect(pSub$var, c(availableCpGs, "Intercept", "Age"))
    pSubValid <- pSub[pSub$var %in% presentVars, ]
    
    # Intercept
    score <- if ("Intercept" %in% pSubValid$var) pSubValid$beta[pSubValid$var == "Intercept"] else 0
    
    # Age effect
    if ("Age" %in% pSubValid$var) {
      score <- score + (pSubValid$beta[pSubValid$var == "Age"] * age)
    }
    
    # CpG effect (Optimized Matrix Multiplication)
    cpgVars <- intersect(pSubValid$var, availableCpGs)
    if (length(cpgVars) > 0) {
      score <- score + as.vector(t(betaM[cpgVars, , drop = FALSE]) %*% 
                                   pSubValid$beta[match(cpgVars, pSubValid$var)])
    }
    protPredList[[pName]] <- score
  }
  
  protDf <- as.data.frame(protPredList)
  
  # --- Phase 2: Calculate Combined Risk Score (COX) ---
  finalInput <- cbind(Age = age, Female = femaleVec, protDf)
  finalInput$Intercept <- 1
  
  availableFinalVars <- intersect(finalModel$var, colnames(finalInput))
  finalWeightsSub <- finalModel[match(availableFinalVars, finalModel$var), ]
  
  coxScore <- as.numeric(as.matrix(finalInput[, availableFinalVars]) %*% finalWeightsSub$beta)
  
  # --- Phase 3: Calibration to Chronological Age ---
  coxParams <- calibParams[calibParams$var == "COX", ]
  ageParams <- calibParams[calibParams$var == "Age", ]
  
  zCox <- (coxScore - coxParams$mean) / coxParams$sd
  finalAgeScore <- (zCox * ageParams$sd) + ageParams$mean
  
  # --- Final Formatting ---
  res <- data.frame(
    SampleID = colnames(betaM),
    protDf,
    stringsAsFactors = FALSE
  )
  res[[outputColName]] <- finalAgeScore
  
  # Rename columns if specified (specifically for GrimAge2 backward compatibility)
  if (applyRenameMap) {
    renameMap <- c(
      "DNAmadm" = "DNAmADM", "DNAmCystatin_C" = "DNAmCystatinC",
      "DNAmGDF_15" = "DNAmGDF15", "DNAmleptin" = "DNAmLeptin",
      "DNAmpai_1" = "DNAmPAI1", "DNAmTIMP_1" = "DNAmTIMP1",
      "DNAmlog.CRP" = "DNAmlogCRP", "DNAmlog.A1C" = "DNAmlogA1C"
    )
    names(res) <- ifelse(names(res) %in% names(renameMap), renameMap[names(res)], names(res))
  }
  
  return(res)
}