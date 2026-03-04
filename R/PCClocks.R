#' @title Calculate Principal Component (PC) Epigenetic Clocks
#'
#' @description
#' Calculates a suite of Principal Component (PC)-based epigenetic clocks
#' based on the methodology from Higgins-Chen et al. (2022).
#'
#' This function computes PC-based versions of Horvath2013, Horvath2018, Hannum,
#' PhenoAge, and GrimAge1, along with their principal components.
#'
#' @param betaM A numeric matrix of DNA methylation (beta) values.
#'   **Samples must be in columns** and CpG probes in rows.
#' @param age A numeric vector of chronological ages for each sample,
#'   in the same order as the columns of `DNAm`.
#' @param sex A character vector of biological sex for each sample, in the
#'   same order as the columns of `DNAm`. Values of "Female" are
#'   encoded as 1; all other values are encoded as 0.
#' @param clockData The pre-loaded data object from 
#'  \code{loadOmniAgeRData("PCClocks_data")}.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage. 
#'  Default is 0.5.
#' @param verbose Logical. Whether to print status messages.
#'
#'
#' @return
#' A data.frame containing the original `Sample_ID`, `Age`, and `Female`
#' columns, appended with 14 new columns for the calculated PC clock values
#' (e.g., `PCHorvath2013`, `PCHannum`, `PCGrimAge1`, etc.).
#'
#' @references
#' Higgins-Chen AT, Thrush KL, Wang Y, et al.
#' A computational solution for bolstering reliability of epigenetic clocks: 
#' Implications for clinical trials and longitudinal tracking.
#' \emph{Nat Aging.} (2022).
#'
#'
#' @export
#'
#' @examples
#'
#' # Download the external data
#' download_OmniAgeR_data(clocks = "PCClocks") #  ZENODO_DOI: "10.5281/zenodo.17162604"
#'
#' # Either path to the data
#' RData <- getOmniAgeRPath()
#' # OR
#' RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
#'
#'
#' downloadOmniAgeRExample("Hannum_example")
#' loadOmniAgeRExample("Hannum_example")
#' age <- PhenoTypesHannum_lv$Age
#' sex <- ifelse(PhenoTypesHannum_lv$Sex == "F", "Female", "Male")
#' pcClocksOut <- pcClocks(hannum_bmiq_m, age, sex, RData)
#'


pcClocks <- function(betaM, age, sex, clockData, minCoverage = 0.5, verbose = TRUE) {
  
  if (verbose) message("[PCClocks] Initializing PC-based clock pipeline...")
  
  # --- 1. Input Validation and Conversion ---
  if (!is.matrix(betaM)) stop("Input 'betaM' must be a matrix.")
  
  # Validate clockData Hash (Security & Integrity Check)
  if (rlang::hash(clockData) != "46386ec4be2b2a5239cf67b242d7dc24") {
    stop("[PCClocks] Invalid or corrupted clockData object. Please re-download.")
  }
  
  pheno <- data.frame(
    SampleID = colnames(betaM), 
    Age = age,
    Female = ifelse(sex == "Female", 1, 0),
    stringsAsFactors = FALSE
  )
  
  # --- 2. Standardized Preprocessing & Coverage Check ---
  # Transpose for PC operations (Rows = Samples)
  betaTrans <- t(betaM)
  
  # Use optimized internal helper for detection and mean imputation
  betaProcessed <- .preprocessPcData(
    betaM = betaTrans,
    requiredCpGs = clockData$imputeMissingCpGs,
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  if (is.null(betaProcessed)) {
    res <- pheno
    res[, 4:17] <- NA_real_ # Fill with NAs if coverage fails
    return(res)
  }
  
  # --- 3. PC Projections and Clock Estimation ---
  # 
  if (verbose) message("[PCClocks] Projecting data onto principal components...")
  
  # Helper to calculate individual PC Clocks
  calcPc <- function(dat, mod, transform = FALSE) {
    # Formula: anti.trafo( (Beta - Center) %*% Rotation %*% Weights + Intercept )
    val <- (sweep(dat, 2, mod$center) %*% mod$rotation %*% mod$model) + mod$intercept
    if (transform) return(as.numeric(.antiTrafo(val)))
    return(as.numeric(val))
  }
  
  pheno$PCHorvath2013 <- calcPc(betaProcessed, clockData$CalcPCHorvath1, transform = TRUE)
  pheno$PCHorvath2018 <- calcPc(betaProcessed, clockData$CalcPCHorvath2, transform = TRUE)
  pheno$PCHannum      <- calcPc(betaProcessed, clockData$CalcPCHannum)
  pheno$PCPhenoAge    <- calcPc(betaProcessed, clockData$CalcPCPhenoAge)
  pheno$PCDNAmTL      <- calcPc(betaProcessed, clockData$CalcPCDNAmTL)
  
  # --- 4. Complex PCGrimAge Logic ---
  if (verbose) message("[PCClocks] Estimating PCGrimAge components...")
  
  # Project beta into PC space for GrimAge
  grimPcSpace <- sweep(betaProcessed, 2, clockData$CalcPCGrimAge$center) %*% clockData$CalcPCGrimAge$rotation
  grimFeatures <- cbind(grimPcSpace, Female = pheno$Female, Age = pheno$Age)
  
  # Internal function for GrimAge Sub-biomarkers
  calcGrimSub <- function(feat, subMod, subInt) {
    as.numeric(feat[, names(subMod)] %*% subMod + subInt)
  }
  
  pheno$PCPACKYRS   <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPACKYRS.model, clockData$CalcPCGrimAge$PCPACKYRS.intercept)
  pheno$PCADM       <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCADM.model, clockData$CalcPCGrimAge$PCADM.intercept)
  pheno$PCB2M       <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCB2M.model, clockData$CalcPCGrimAge$PCB2M.intercept)
  pheno$PCCystatinC <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCCystatinC.model, clockData$CalcPCGrimAge$PCCystatinC.intercept)
  pheno$PCGDF15     <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCGDF15.model, clockData$CalcPCGrimAge$PCGDF15.intercept)
  pheno$PCLeptin    <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCLeptin.model, clockData$CalcPCGrimAge$PCLeptin.intercept)
  pheno$PCPAI1      <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPAI1.model, clockData$CalcPCGrimAge$PCPAI1.intercept)
  pheno$PCTIMP1     <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCTIMP1.model, clockData$CalcPCGrimAge$PCTIMP1.intercept)
  
  # Final integrated PCGrimAge1
  grimComp <- pheno[, clockData$CalcPCGrimAge$components]
  pheno$PCGrimAge1 <- as.numeric(as.matrix(grimComp) %*% clockData$CalcPCGrimAge$PCGrimAge.model + clockData$CalcPCGrimAge$PCGrimAge.intercept)
  
  return(pheno)
}





