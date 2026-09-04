#' @title Predict DNA methylation age using CTS clocks
#' @description
#' This is the function for computing DNA methylation age using CTS (cell type
#' specific) clocks. The inputs include a DNAm matrix (or multi-omics object), 
#' the CTS clocks you want to use, the data type (bulk vs sorted), cell type 
#' fraction matrix (for Neu-In/Glia-In/Brain clocks), tissue type, and parallel 
#' computing options.
#'
#' @details
#' This function supports a variety of Cell-Type-Specific (CTS) clocks.
#'
#' **Available `CTSclocks` include:**
#' * `'Neu-In'`
#' * `'Glia-In'`
#' * `'Brain'`
#' * `'Neu-Sin'`
#' * `'Glia-Sin'`
#' * `'Hep'`
#' * `'Liver'`
#'
#'  The clocks are grouped below based on their biological target:
#'
#' **1. Cell-Type Specific Clocks**
#'
#' These clocks are trained to measure aging in specific cell populations.
#'
#' * `'Neu-In'` (Intrinsic): Measures cell-intrinsic aging of **Neurons**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Glia-In'` (Intrinsic): Measures cell-intrinsic aging of **Glial cells**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Neu-Sin'` (Semi-intrinsic): Measures aging of **Neurons** using raw data.
#'     (Reflects both intrinsic aging and composition changes).
#' * `'Glia-Sin'` (Semi-intrinsic): Measures aging of **Glial cells** using raw data.
#'     (Reflects both intrinsic aging and composition changes).
#' * `'Hep'` (Semi-intrinsic): Measures aging of **Hepatocytes** using raw data.
#'
#' **2. Non Cell-Type Specific Clocks**
#'
#' * `'Brain'` (Intrinsic): A intrinsic clock for **whole brain tissue**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Liver'` (Semi-intrinsic): A intrinsic clock for **whole liver tissue**.
#'     (Uses raw data).
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment}, 
#'   object containing DNA methylation beta values. Rows should be CpGs and 
#'   columns should be samples.
#' @param compClocks A character vector of one or more clocks to apply.
#'                  (e.g., 'Neu-In', 'Hep', c('Neu-In', 'Neu-Sin', 'Brain')).
#' @param dataType Type of the samples ('bulk' or 'sorted').
#' @param ctfM Optional cell type fraction matrix
#' (rows: samples, columns: cell types).
#'             Required for 'Intrinsic' bulk clocks if tissue is not 'brain'.
#' @param tissue What tissue are your samples from ('brain' or 'otherTissue').
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose Logical. If `TRUE` (default), the function will print
#'   progress messages to the console.
#'
#' @return A data.frame of predicted DNAm ages (samples x clocks).
#'
#' @importFrom stats coef lm sd
#'
#' @export
#'
#' @references
#' Tong H, Guo X, Jacques M, Luo Q, Eynon N, Teschendorff AE.
#' Cell-type specific epigenetic clocks to quantify biological age
#' at cell-type resolution.
#' \emph{Aging} 2024
#'
#' @examples
#' # ====================================================================
#' # Example 1: Direct Matrix Input (Bulk Tissue)
#' # ====================================================================
#' murphyBetaM <- OmniAgeRData::getOmniAgeRData(
#'     "omniager_cts_murphy_gse88890",
#'     verbose = FALSE
#' )[[1]]
#'
#' agePred_df <- ctsClocks(
#'     x = murphyBetaM,
#'     compClocks = c("Neu-In", "Neu-Sin"),
#'     dataType = "bulk",
#'     ctfM = NULL,
#'     tissue = "brain",
#'     verbose = FALSE
#' )
#' # ====================================================================
#' # Example 2: Direct Matrix Input (Sorted Cells)
#' # ====================================================================
#' \dontrun{
#' paiBetaM <- OmniAgeRData::getOmniAgeRData(
#'     "omniager_cts_pai_gse112179",
#'     verbose = FALSE
#' )[[1]]
#' 
#' agePred_df_sorted <- ctsClocks(
#'     x = paiBetaM,
#'     compClocks = c("Neu-In", "Neu-Sin"),
#'     dataType = "sorted",
#'     ctfM = NULL,
#'     tissue = "brain",
#'     verbose = FALSE
#' )
#' 
#' # ====================================================================
#' # Example 3: SummarizedExperiment Input
#' # ====================================================================
#' library(SummarizedExperiment)
#' se_obj <- SummarizedExperiment(
#'   assays = list(beta = murphyBetaM)
#' )
#'     
#' agePred_se <- ctsClocks(
#'    x = se_obj,
#'    compClocks = c("Neu-In", "Neu-Sin"),
#'    dataType = "bulk",
#'    ctfM = NULL,
#'    tissue = "brain",
#'    verbose = FALSE
#' )
#' }
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/HGT-UwU/CTSclocks
# under the GPL-3.0 License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------
ctsClocks <- function(x,
                      compClocks = c("Neu-In"),
                      dataType = c("bulk", "sorted"),
                      ctfM = NULL,
                      tissue = c("brain", "otherTissue"),
                      minCoverage = 0,
                      verbose = TRUE) {
  
  dataType <- match.arg(dataType)
  tissue <- match.arg(tissue)
  
  # --- Step 0: Universal Matrix Extraction ---
  betaM <- .extractAssayMatrix(x)
  
  ctsClocksCoef <- OmniAgeRData::getOmniAgeRData(
    "omniager_cts_clocks_coef",
    verbose = verbose
  )
  
  # --- 1. Perform full data processing (deconvolution) ---
  needsIntrinsic <- any(c("Neu-In", "Glia-In", "Brain") %in% compClocks)
  
  if (needsIntrinsic && dataType == "bulk" && is.null(ctfM)) {
    if (tissue == "brain") {
      if (verbose) message("[CTS] Deconvolving brain tissue fractions using full matrix...")
      estF <- HiBED::HiBED_deconvolution(betaM, h = 1) / 100
      ctfM <- as.matrix(estF[, c(3, 2, 1)]) # Neu, Glia, EndoStrom
    } else {
      stop("ctfM is required for non-brain bulk tissue.")
    }
  }
  
  # --- 2. Probe Selection & Coverage Check ---
  allClockCpGs <- unique(unlist(lapply(ctsClocksCoef[compClocks], function(modelObj) {
    if (inherits(modelObj, "glmnet")) {
      coefMat <- as.matrix(stats::coef(modelObj))
      probes <- rownames(coefMat)[rowSums(coefMat != 0) > 0]
      return(probes[probes != "(Intercept)"])
    }
    probes <- modelObj$probe[modelObj$coef != 0]
    return(probes[probes != "(Intercept)"])
  })))
  
  dummyWeights <- setNames(rep(1, length(allClockCpGs)), allClockCpGs)
  coverage <- .checkCpGCoverage(
    betaM = betaM,
    allWeights = dummyWeights,
    clockName = "CTS_Global",
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  if (!coverage$pass) {
    if (verbose) warning("[CTS] Calculation aborted: Data integrity check failed.")
    return(as.data.frame(matrix(NA, nrow = ncol(betaM), ncol = length(compClocks)),
                         row.names = colnames(betaM), col.names = compClocks
    ))
  }
  
  # --- 3. Prepare the data subset ---
  presentBetaMat <- betaM[coverage$betaIdx, , drop = FALSE]
  
  # --- 4. Perform residual calculation and standardization ---
  processedMat <- NULL
  if (needsIntrinsic) {
    processedMat <- .processCtsData(
      betaM = presentBetaMat,
      dataType = dataType,
      tissue = tissue,
      ctfM = ctfM,
      verbose = verbose
    )
  }
  
  # --- 5. Iterative calculation of predicted values ---
  resultsList <- list()
  for (clockLabel in compClocks) {
    modelObj <- ctsClocksCoef[[clockLabel]]
    
    targetMat <- if (clockLabel %in% c("Neu-In", "Glia-In", "Brain")) processedMat else presentBetaMat
    
    if (inherits(modelObj, "glmnet")) {
      coefMat <- as.matrix(stats::coef(modelObj))
      nonzero_idx <- rowSums(coefMat != 0) > 0
      nonzero_coefs <- coefMat[nonzero_idx, , drop = FALSE]
      
      intercept <- nonzero_coefs["(Intercept)", 1]
      probe_weights <- nonzero_coefs[rownames(nonzero_coefs) != "(Intercept)", 1]
      weights <- setNames(as.numeric(probe_weights), names(probe_weights))
    } else {
      intercept <- modelObj$coef[1]
      weights <- setNames(modelObj$coef[-1], modelObj$probe[-1])
    }
    
    resultsList[[clockLabel]] <- .calculateLinearPredictor(
      betaM = targetMat,
      coefLv = list(intercept, weights),
      clockName = clockLabel,
      minCoverage = 0, 
      verbose = FALSE
    )
  }
  
  return(as.data.frame(resultsList, row.names = colnames(betaM)))
}


#' Internal Data Preprocessing for CTS Intrinsic Clocks
#'
#' @description
#' Prepares methylation data for 'Intrinsic' clocks (Neu-In, Glia-In, Brain)
#' by removing extrinsic aging factors (cell composition) and standardizing
#' the data.
#'
#' @param betaM A numeric matrix of methylation beta values.
#' @param dataType Character, either "bulk" or "sorted".
#' @param tissue Character, the tissue source of the data.
#' @param ctfM A numeric matrix of cell type fractions.
#' @param verbose Logical.
#'
#' @return A numeric matrix of processed (residuals/Z-scores) methylation values.
#'
#' @importFrom stats lm sd
#' @keywords internal
#' @noRd
.processCtsData <- function(betaM, dataType, tissue, ctfM, verbose) {
  if (dataType == "sorted") {
    # Z-score standardization across samples
    rowMeansV <- rowMeans(betaM, na.rm = TRUE)
    rowSdsV <- apply(betaM, 1, stats::sd, na.rm = TRUE)
    rowSdsV[rowSdsV == 0] <- 1 # [极其关键的安全阀] 防止方差为 0 导致出现 NaN
    return((betaM - rowMeansV) / rowSdsV)
  } else if (dataType == "bulk") {
    if (verbose) message("[CTS] Regressing out cell type effects from subset matrix...")
    
    # Linear model: methylation ~ cell type fractions
    fit <- stats::lm(t(betaM) ~ ctfM)
    resM <- t(fit$residuals)
    
    # Standardize Residuals
    rowMeansR <- rowMeans(resM, na.rm = TRUE)
    rowSdsR <- apply(resM, 1, stats::sd, na.rm = TRUE)
    rowSdsR[rowSdsR == 0] <- 1
    return((resM - rowMeansR) / rowSdsR)
  }
}
