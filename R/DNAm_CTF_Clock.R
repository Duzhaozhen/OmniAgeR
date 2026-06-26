#' DNA Methylation Cell-Type Fraction (CTF) Aging Clock
#'
#' @description
#' Predicts biological age based on immune cell type fractions derived from
#' DNA methylation data.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment}, 
#'   object containing celltype fraction matrix. Rows should be cell types and 
#'   columns should be samples. Must contain the specific cell types required 
#'   by the model (e.g., predicted by EpiDISH).
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector of predicted ages.
#' @export
#' @importFrom randomForest randomForest importance
#'
#'
#' @examples
#' 
#' ctfExample <- loadOmniAgeRdata(
#'     "omniager_tzh_example_ctf",
#'     verbose = FALSE
#' )
#' dnamCTFClockOut <- dnamCTFClock(t(ctfExample[[2]]))
#'
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     pheno_data <- ctfExample[[1]]
#'     rownames(pheno_data) <- pheno_data[["Sample"]]
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(ctf=t(ctfExample[[2]])),
#'       colData = pheno_data
#'     )
#'     
#'     dnamCTFClockOut <- dnamCTFClock(x = se_obj, verbose = FALSE)
#'   }
#' }
#' 

dnamCTFClock <- function(x, verbose = TRUE) {
    # --- 1. Load the internal model ---
    ctfM <- t(.extractAssayMatrix(x))
    dnamCtfModel <- loadOmniAgeRdata(
        "omniager_dnam_ctf_model",
        verbose = verbose
    )
    # --- 2. Verify feature integrity ---
    requiredFeatures <- rownames(dnamCtfModel$importance)

    missingCols <- setdiff(requiredFeatures, colnames(ctfM))
    if (length(missingCols) > 0) {
        stop(
            "[dnamCTFClock] Missing required cell types: ",
            paste(missingCols, collapse = ", "),
            ". Ensure you provide estimated fractions for all required types."
        )
    }

    # --- 3. Data alignment and NA inspection ---
    dataForPred <- ctfM[, requiredFeatures, drop = FALSE]

    if (any(is.na(dataForPred))) {
        stop(
            "[dnamCTFClock] Input contains NA values. ",
            "Random Forest model requires complete data."
        )
    }

    # --- 5. Prediction ---
    predAge <- stats::predict(dnamCtfModel, newdata = dataForPred)

    if (!is.null(rownames(ctfM))) {
        names(predAge) <- rownames(ctfM)
    }

    return(predAge)
}
