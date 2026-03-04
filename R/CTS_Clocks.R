#' @title Predict DNA methylation age using CTS clocks
#' @description
#' This is the function for computing DNA methylation age using CTS (cell type
#' specific) clocks. The inputs including DNAm matrix, the CTS clock you want to
#' use, is your DNAm data from bulk tissue sample or sorted cell sample, cell
#' type fraction matrix if you want to use Neu-In/Glia-In/Brain clock, tissue of
#' your DNAm data samples and the number of cores if you want to do parallel computing.
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
#''
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
#' @param betaM  A DNAm matrix (row: CpGs, column: samples) of the samples you 
#' want to get  DNAm age predicted by a CTS clock.
#' @param CTSclocks A character vector of one or more clocks to apply.
#'                  (e.g., 'Neu-In', 'Hep', c('Neu-In', 'Neu-Sin', 'Brain')).
#' @param dataType Type of the samples ('bulk' or 'sorted').
#' @param ctfM Optional cell type fraction matrix 
#' (rows: samples, columns: cell types).
#'              Required for 'Intrinsic' bulk clocks if tissue is not 'brain'.
#' @param tissue What tissue are your samples from ('brain' or 'otherTissue').
#'
#' @param verbose Logical. If `TRUE` (default), the function will print
#'   progress messages to the console.
#'
#' @return A data.frame of predicted DNAm ages (samples x clocks).
#'
#' @import glmnet
#' @import parallel
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
#' downloadOmniAgeRExample("CTS_MurphyGSE88890")
#'loadOmniAgeRExample("CTS_MurphyGSE88890")
#'agePred_df <- ctsClocks(
#'  Murphy_beta_m,
#'  ctsClocks = c("Neu-In", "Neu-Sin"),
#'  dataType = "bulk",
#'  ctfM = NULL,
#'  tissue = "brain"
#'  )
#' downloadOmniAgeRExample("CTS_PaiGSE112179")
#' loadOmniAgeRExample("CTS_PaiGSE112179")
#' agePred_df <- ctsClocks(
#'   Pai_beta_m,
#'   ctsClocks = c("Neu-In", "Neu-Sin"),
#'   dataType = "sorted",
#'   ctfM = NULL,
#'   tissue = "brain"
#' )
#'
#' downloadOmniAgeRExample("CTS_ExampleData_Liver")
#' loadOmniAgeRExample("CTS_ExampleData_Liver")
#' agePred_df <- ctsClocks(
#'   Liver_beta_m,
#'   ctsClocks = c("Hep", "Liver"),
#'   dataType = "bulk",
#'   ctfM = NULL,
#'   tissue = "otherTissue"
#' )
#'
#' @export
ctsClocks <- function(betaM, 
                      ctsClocks = c("Neu-In"), 
                      dataType = c("bulk", "sorted"),
                      ctfM = NULL, 
                      tissue = c("brain", "otherTissue"),
                      minCoverage = 0.5, 
                      verbose = TRUE) {
  
  dataType <- match.arg(dataType)
  tissue <- match.arg(tissue)
  
  data("CTS_Clocks_Coef", envir = environment())
  # --- 1. Perform full data processing (deconvolution) ---
  needsIntrinsic <- any(c("Neu-In", "Glia-In", "Brain") %in% ctsClocks)
  
  if (needsIntrinsic && dataType == "bulk" && is.null(ctfM)) {
    if (tissue == "brain") {
      if (verbose) message("[CTS] Deconvolving brain tissue fractions using full matrix...")
      estF <- HiBED::HiBED_deconvolution(betaM, h = 1) / 100
      ctfM <- as.matrix(estF[, c(3, 2, 1)]) # Neu, Glia, EndoStrom
    } else {
      stop("ctfM is required for non-brain bulk tissue.")
    }
  }
  
  # --- 2. Deconvolution ---
  allClockCpGs <- unique(unlist(lapply(CTS_Clocks_Coef[ctsClocks], function(x) {
    if (methods::is(x, "glmnet")) return(rownames(x$beta))
    return(x$probe)
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
    return(as.data.frame(matrix(NA, nrow = ncol(betaM), ncol = length(ctsClocks)),
                         row.names = colnames(betaM), col.names = ctsClocks))
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
  for (clockLabel in ctsClocks) {
    modelObj <- CTS_Clocks_Coef[[clockLabel]]
    
    if (clockLabel %in% c("Neu-In", "Glia-In", "Brain", "Neu-Sin", "Glia-Sin")) {
      targetMat <- if (clockLabel %in% c("Neu-In", "Glia-In", "Brain")) processedMat else presentBetaMat
      intercept <- modelObj$coef[1]
      weights <- setNames(modelObj$coef[-1], modelObj$probe[-1])
      
      resultsList[[clockLabel]] <- .calculateLinearPredictor(
        betaM = targetMat, 
        coefLv = list(intercept, weights), 
        clockName = clockLabel, 
        minCoverage = 0, 
        verbose = FALSE
      )
    } else {
      commonCpGs <- intersect(rownames(modelObj$beta), rownames(presentBetaMat))
      if (length(commonCpGs) < (length(rownames(modelObj$beta)) * minCoverage)) {
        resultsList[[clockLabel]] <- rep(NA_real_, ncol(betaM))
      } else {
        trimmedBeta <- presentBetaMat[commonCpGs, , drop = FALSE]
        resultsList[[clockLabel]] <- as.vector(stats::predict(modelObj, newx = t(trimmedBeta)))
      }
    }
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
#' @details
#' The processing steps depend on the `dataType`:
#' \itemize{
#'   \item \strong{sorted}: Performs row-wise Z-score standardization. This 
#'   assumes the samples are already pure cell types.
#'   \item \strong{bulk}: Fits a linear model (\code{DNAm ~ cell_type_fractions}) 
#'   and extracts the residuals. These residuals represent methylation 
#'   changes independent of cell composition. The residuals are then 
#'   Z-score standardized.
#' }
#'
#' @param betaM A numeric matrix of methylation beta values (probes x samples).
#' @param dataType Character, either "bulk" or "sorted".
#' @param tissue Character, the tissue source of the data.
#' @param ctfM A numeric matrix of cell type fractions (samples x cell types). 
#' Required if \code{dataType} is "bulk".
#' @param verbose Logical, whether to display progress messages.
#'
#' @return A numeric matrix of processed (residuals/Z-scores) methylation values.
#'
#' @importFrom stats lm sd
#' @keywords internal
#' @noRd

.processCtsData <- function(betaM, dataType, tissue, ctfM, verbose) {
  
  if (dataType == "sorted") {
    # Z-score standardization across samples
    rowMeansV <- rowMeans(betaM)
    rowSdsV <- apply(betaM, 1, stats::sd)
    return((betaM - rowMeansV) / rowSdsV)
    
  } else if (dataType == "bulk") {
    # Here, ctfM is guaranteed to be present (either from user or calculated in main function)
    if (verbose) message("[CTS] Regressing out cell type effects from subset matrix...")
    
    # Linear model: methylation ~ cell type fractions
    # Use residuals to isolate cell-intrinsic aging
    fit <- stats::lm(t(betaM) ~ ctfM)
    resM <- t(fit$residuals)
    
    # Standardize Residuals
    rowMeansR <- rowMeans(resM)
    rowSdsR <- apply(resM, 1, stats::sd)
    return((resM - rowMeansR) / rowSdsR)
  }
}

