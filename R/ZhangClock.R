#' @title Predicts DNA methylation age using the Zhang clock(Elastic Net model)
#' #'
#' @description
#' This function takes a matrix of DNA methylation beta values and calculates
#' the epigenetic age for each sample based on the elastic net model developed
#' by Zhang et al. (2019).
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' This function implements the elastic net epigenetic clock from
#' Zhang et al. (2019). A unique feature of this clock is the
#' pre-processing step. Instead of using raw beta values, it
#' first performs a per-sample standardization. Specifically, it calculates a
#' Z-score for each CpG based on the mean and standard deviation of all
#' measured CpGs within that same sample. The final age is then predicted
#' by taking the weighted sum of these standardized values using the model's
#' pre-defined coefficients.
#'
#' @return A numeric vector of predicted ages, with sample names preserved.
#'
#'
#' @references
#' Zhang Q, Vallerga CL, Walker RM, et al.
#' Improved precision of epigenetic clock estimates across tissues and
#' its implication for biological ageing.
#' \emph{Genome Med.} 2019
#'
#' @export
#'
#' @examples
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#' 
#' # Example 1: Direct Matrix Input
#' predOut <- zhangClock(x = beta_matrix, verbose = FALSE)
#' 
#' # Example 2: SummarizedExperiment Input
#' \dontrun{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     library(SummarizedExperiment)
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     predOut <- zhangClock(x = se_obj, verbose = FALSE)
#'   }
#' }
#'

zhangClock <- function(x,
                       minCoverage = 0,
                       verbose = TRUE) {
    
    zhangClockCoef <- loadOmniAgeRdata(
        "omniager_zhang_clock_coef",
        verbose = verbose
    )
    betaM <- .extractAssayMatrix(x)
    # --- 1. Per-Sample Standardization ---
    if (verbose) message("[zhangClock] Performing per-sample standardization...")

    sample_means <- colMeans(betaM, na.rm = TRUE)
    sample_sds <- apply(betaM, 2, sd, na.rm = TRUE)

    # Safe handling: Prevent division by zero errors caused by a standard deviation of 0
    sample_sds[sample_sds == 0] <- 1

    # Z-score
    beta_scaled <- sweep(betaM, 2, sample_means, "-")
    beta_scaled <- sweep(beta_scaled, 2, sample_sds, "/")

    # --- 2. Predict Age ---
    predicted_age <- .calLinearClock(
        betaM       = beta_scaled,
        coefData    = zhangClockCoef,
        clockLabel  = "zhangClock",
        minCoverage = minCoverage,
        verbose     = verbose
    )

    return(predicted_age)
}
