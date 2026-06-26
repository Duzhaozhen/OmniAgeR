#' @title Calculate Principal Component (PC) Epigenetic Clocks
#'
#' @description
#' Calculates a suite of Principal Component (PC)-based epigenetic clocks
#' based on the methodology from Higgins-Chen et al. (2022).
#'
#' This function computes PC-based versions of Horvath2013, Horvath2018, Hannum,
#' PhenoAge, and GrimAge1, along with their principal components.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param age A numeric vector of chronological ages for each sample,
#'   in the same order as the columns of `DNAm`.
#' @param sex A character vector of biological sex for each sample, in the
#'   same order as the columns of `DNAm`. Values of "Female" are
#'   encoded as 1; all other values are encoded as 0.
#' @param clockData The pre-loaded data object from
#'  \code{loadOmniAgeRData("PCClocks_data")}.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'  Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#'
#' @return
#' A data.frame containing the original `SampleID` columns, 
#' appended with 14 new columns for the calculated PC clock values
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
#' \dontrun{
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#'   hannumExample <- loadOmniAgeRdata(
#'       "omniager_hannum_example",
#'       verbose = FALSE
#'   )
#'   
#'   pcClockData <- loadOmniAgeRdata(
#'       "PCClocks_data",
#'       verbose = FALSE
#'   )
#'   
#'   hannumBmiqM <- hannumExample[[1]]
#'   phenoTypesHannum <- hannumExample[[2]]
#'   age <- phenoTypesHannum$Age
#'   sex <- ifelse(phenoTypesHannum$Sex == "F", "Female", "Male")
#'   pcClocksOut <- pcClocks(hannumBmiqM, age, sex, pcClockData)
#'   
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = hannumBmiqM),
#'       colData = hannumExample[[2]]
#'     )
#'     
#'     pcClocksOut <-  pcClocks(se_obj,age, sex, pcClockData)
#'   }
#' }
#'
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/HigginsChenLab/methylCIPHER
# under the BSD 3-Clause License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------
pcClocks <- function(x, age, sex, clockData, minCoverage = 0, verbose = TRUE) {
    if (verbose) message("[PCClocks] Initializing PC-based clock pipeline...")
    betaM <- .extractAssayMatrix(x)
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
        if (transform) {
            return(as.numeric(.antiTrafo(val)))
        }
        return(as.numeric(val))
    }

    pheno$PCHorvath2013 <- calcPc(betaProcessed, clockData$CalcPCHorvath1, transform = TRUE)
    pheno$PCHorvath2018 <- calcPc(betaProcessed, clockData$CalcPCHorvath2, transform = TRUE)
    pheno$PCHannum <- calcPc(betaProcessed, clockData$CalcPCHannum)
    pheno$PCPhenoAge <- calcPc(betaProcessed, clockData$CalcPCPhenoAge)
    pheno$PCDNAmTL <- calcPc(betaProcessed, clockData$CalcPCDNAmTL)

    # --- 4. Complex PCGrimAge Logic ---
    if (verbose) message("[PCClocks] Estimating PCGrimAge components...")
    # Project beta into PC space for GrimAge
    grimPcSpace <- sweep(betaProcessed, 2, clockData$CalcPCGrimAge$center) %*% clockData$CalcPCGrimAge$rotation
    grimFeatures <- cbind(grimPcSpace, Female = pheno$Female, Age = pheno$Age)
    # Internal function for GrimAge Sub-biomarkers
    calcGrimSub <- function(feat, subMod, subInt) {
        as.numeric(feat[, names(subMod)] %*% subMod + subInt)
    }
    pheno$PCPACKYRS <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPACKYRS.model, clockData$CalcPCGrimAge$PCPACKYRS.intercept)
    pheno$PCADM <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCADM.model, clockData$CalcPCGrimAge$PCADM.intercept)
    pheno$PCB2M <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCB2M.model, clockData$CalcPCGrimAge$PCB2M.intercept)
    pheno$PCCystatinC <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCCystatinC.model, clockData$CalcPCGrimAge$PCCystatinC.intercept)
    pheno$PCGDF15 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCGDF15.model, clockData$CalcPCGrimAge$PCGDF15.intercept)
    pheno$PCLeptin <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCLeptin.model, clockData$CalcPCGrimAge$PCLeptin.intercept)
    pheno$PCPAI1 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPAI1.model, clockData$CalcPCGrimAge$PCPAI1.intercept)
    pheno$PCTIMP1 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCTIMP1.model, clockData$CalcPCGrimAge$PCTIMP1.intercept)
    # Final integrated PCGrimAge1
    grimComp <- pheno[, clockData$CalcPCGrimAge$components]
    pheno$PCGrimAge1 <- as.numeric(as.matrix(grimComp) %*% clockData$CalcPCGrimAge$PCGrimAge.model + clockData$CalcPCGrimAge$PCGrimAge.intercept)
    
    pheno$Age <-NULL
    pheno$Female <-NULL
    return(pheno)
}
