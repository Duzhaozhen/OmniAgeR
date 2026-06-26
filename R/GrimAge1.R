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




grimAge1 <- function(x, age, sex,
                     minCoverage = 0, verbose = TRUE) {
    # --- Step 0: Universal Matrix Extraction ---
    betaM <- .extractAssayMatrix(x)
    # 1. Load model weights
    grimage1 <- loadOmniAgeRdata(
        "omniager_grimage1_model",
        verbose = verbose
    )
    # 2. Extract coefficients from model object (grimage1)
    protCoefs <- grimage1[[1]] # CpG weights for proteins
    finalModel <- grimage1[[2]] # Final weights for COX
    calibParams <- grimage1[[3]] # Calibration means/sds
    uqCpgs <- unique(protCoefs$var[startsWith(protCoefs$var, "cg")])
    fakeWeights <- setNames(rep(1, length(uqCpgs)), uqCpgs)
    coverage <- .checkCpGCoverage(betaM, fakeWeights, "GrimAge1", minCoverage, verbose)
    # 3. Handle Covariates (Age and Sex)
    femaleVec <- ifelse(sex == "Female", 1, 0)
    # 4. Phase 1: Predict DNAm Protein Biomarkers
    availableCpGs <- intersect(protCoefs$var, rownames(betaM))
    proteinNames <- unique(protCoefs$Y.pred)
    protPredList <- list()
    for (pName in proteinNames) {
        # Subset coefficients for this specific protein
        pSub <- protCoefs[protCoefs$Y.pred == pName, ]
        # Intersection with available data
        presentVars <- intersect(pSub$var, c(availableCpGs, "Intercept", "Age"))
        pSubValid <- pSub[pSub$var %in% presentVars, ]
        # Calculate score using only present features
        score <- 0
        if ("Intercept" %in% pSubValid$var) {
            score <- pSubValid$beta[pSubValid$var == "Intercept"]
        }
        # Age component (if required)
        if ("Age" %in% pSubValid$var) {
            score <- score + (pSubValid$beta[pSubValid$var == "Age"] * age)
        }
        # CpG component: matrix multiplication of present sites
        cpgVars <- intersect(pSubValid$var, availableCpGs)
        if (length(cpgVars) > 0) {
            # t(betaM) ensures samples are rows for the multiplication
            score <- score + as.vector(t(betaM[cpgVars, , drop = FALSE]) %*%
                pSubValid$beta[match(cpgVars, pSubValid$var)])
        }
        protPredList[[pName]] <- score
    }
    protDf <- as.data.frame(protPredList)
    # 5. Phase 2: Calculate Mortality Risk Score (COX)
    # Feature set: Age, Female, and predicted DNAm Proteins
    finalInput <- cbind(Age = age, Female = femaleVec, protDf)
    finalInput$Intercept <- 1
    # Match variables for the final COX model
    availableFinalVars <- intersect(finalModel$var, colnames(finalInput))
    finalWeightsSub <- finalModel[match(availableFinalVars, finalModel$var), ]
    coxScore <- as.numeric(as.matrix(finalInput[, availableFinalVars]) %*% finalWeightsSub$beta)
    # 6. Phase 3: Calibration to Chronological Age
    coxParams <- calibParams[calibParams$var == "COX", ]
    ageParams <- calibParams[calibParams$var == "Age", ]

    # Formula: DNAmGrimAge = ((COX - MeanCOX)/SdCOX * SdAge) + MeanAge
    zCox <- (coxScore - coxParams$mean) / coxParams$sd
    grimAgeScore <- (zCox * ageParams$sd) + ageParams$mean

    # 7. Final Output
    res <- data.frame(
        SampleID = colnames(betaM),
        protDf,
        DNAmGrimAge1 = grimAgeScore,
        stringsAsFactors = FALSE
    )

    return(res)
}
