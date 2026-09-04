#' @title Calculate GrimAge2
#'
#' @description
#' Calculates DNA methylation GrimAge2, a composite biomarker of mortality risk
#' and biological aging. This function implements the updated GrimAge2 model.
#'
#' @details
#' This function calculates DNAm GrimAge2 in a multi-step process. First, it
#' predicts DNAm-based surrogate biomarkers for several plasma proteins from
#' the input beta values. These predicted biomarkers, along with chronological
#' age and sex, are then used to calculate a composite mortality risk score.
#' This score is calibrated to the scale of chronological age to produce the
#' final `DNAmGrimAge2`.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing DNA methylation beta values. Rows should be CpG probes and 
#'  columns individual samples.
#' @param age A numeric vector of chronological ages for the samples corresponding
#'   to the columns in `betaM`.
#' @param sex A character vector of sample sexes. Must contain "Male" or "Female"
#'   for each sample.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @return
#' A data.frame containing the following columns:
#' \itemize{
#'   \item `SampleID`: Identifier for each sample.
#'   \item `DNAm...`: Columns for each of the predicted surrogate biomarkers
#'   (e.g., `DNAmADM`, `DNAmGDF15`).
#'   \item `DNAmGrimAge2`: The final calibrated GrimAge2 score.
#' }
#'
#'
#'
#' @references
#' Lu AT, Binder AM, Zhang J, et al.
#' DNA methylation GrimAge version 2.
#' \emph{Aging} 2022
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
#' GrimAge2O <- grimAge2(x = beta_matrix, age = age, sex = sex, verbose = FALSE)
#' 
#' \dontrun{
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#'   library(SummarizedExperiment)
#'
#'   pheno_data <- dnamExample[[2]]
#'   rownames(pheno_data) <- colnames(beta_matrix)
#'
#'   se_obj <- SummarizedExperiment(
#'     assays = list(beta = beta_matrix),
#'     colData = pheno_data
#'   )
#'
#'   GrimAge2O <- grimAge2(
#'     x = se_obj,
#'     age = se_obj$Age,
#'     sex = ifelse(se_obj$Sex == "F", "Female", "Male"),
#'     verbose = FALSE
#'   )
#' }

grimAge2 <- function(x, age, sex, minCoverage = 0, verbose = TRUE) {
  .calculateGrimAge(
    x = x, 
    age = age, 
    sex = sex, 
    minCoverage = minCoverage, 
    verbose = verbose,
    modelName = "omniager_grimage2_model",
    clockName = "GrimAge2",
    outputColName = "DNAmGrimAge2",
    applyRenameMap = TRUE
  )
}
