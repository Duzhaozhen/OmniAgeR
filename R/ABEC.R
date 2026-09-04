#' @title Adult Blood-based EPIC Clock (ABEC)
#'
#' @description
#' Predicts biological age using the Adult Blood-based EPIC Clock (ABEC).
#' Developed by Lee et al., this model was trained on DNA methylation (DNAm)
#' data from the Norwegian Mother, Father and Child Cohort Study (MoBa)
#' (n = 1,592, age range: 19–59 years) using the Illumina EPIC platform.
#'
#' @details
#' The function extracts the necessary CpG coefficients from the internal
#' \code{ABEC_Coef} dataset and applies them to the provided beta value matrix.
#' It uses an internal helper to handle missing probes and compute the
#' final age estimates.
#'
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#' object containing DNA methylation beta values. Rows should be CpG probes and 
#' columns individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#'
#' @return A numeric vector of predicted biological ages.
#'
#' @export
#'
#' @references
#' Lee, Y., Haftorn, K.L., Denault, W.R.P. et al.
#' Blood-based epigenetic estimators of chronological age in human adults
#' using DNA methylation data from the Illumina MethylationEPIC array.
#' \emph{BMC Genomics} 2020
#'
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predAge <- leeABEC(x = beta_matrix)
#' head(predAge)
#'
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'     library(SummarizedExperiment)
#'     
#'     # Extract phenotype and ensure rownames match matrix colnames
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     # Construct the SummarizedExperiment object
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     # The function seamlessly accepts the Bioconductor object
#'     predAge <- leeABEC(x = se_obj)
#'     head(predAge)
#' }

leeABEC <- function(x,
                    minCoverage = 0,
                    verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_abec_coef", 
    clockName = "leeABEC",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}


#' @title Extended Adult Blood-based EPIC Clock (eABEC)
#'
#' @description
#' Predicts biological age using the Extended Adult Blood-based EPIC Clock
#' (eABEC). This model extends the training set of ABEC by incorporating
#' public data from the Gene Expression Omnibus (GEO), resulting in a
#' broader age-span (n = 2,227, age range: 18–88 years).
#'
#' @inheritParams leeABEC
#' @inherit leeABEC return
#'
#' @details
#' Similar to \code{leeABEC}, this function utilizes the \code{eABEC_Coef}
#' dataset. It is designed for applications where a wider range of adult
#' ages is expected.
#'
#' @export
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predAge <- leeExtendedABEC(x = beta_matrix)
#' head(predAge)
#'
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'     library(SummarizedExperiment)
#'     
#'     # Extract phenotype and ensure rownames match matrix colnames
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     # Construct the SummarizedExperiment object
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     # The function seamlessly accepts the Bioconductor object
#'     predAge <- leeExtendedABEC(x = se_obj)
#'     head(predAge)
#' }
leeExtendedABEC <- function(x,
                            minCoverage = 0,
                            verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_eabec_coef", 
    clockName = "leeExtendedABEC",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}


#' @title Common Adult Blood-based EPIC Clock (cABEC)
#'
#' @description
#' Predicts biological age using the Common Adult Blood-based EPIC Clock
#' (cABEC). This model uses the same extended training set as \code{eABEC}
#' but is restricted to CpGs common to both Illumina 450K and EPIC arrays,
#' ensuring backward compatibility and robustness across platforms.
#'
#' @inheritParams leeABEC
#' @inherit leeABEC return
#'
#' @details
#' The function uses coefficients from the \code{cABEC_Coef} dataset.
#'
#' @export
#' @examples
#' # Load the included example data (a list containing a matrix and phenotypes)
#' data(dnamExample)
#' beta_matrix <- dnamExample[[1]]
#'
#' # ====================================================================
#' # Example 1: Direct Matrix Input
#' # ====================================================================
#' predAge <- leeCommonABEC(x = beta_matrix)
#' head(predAge)
#'
#' # ====================================================================
#' # Example 2: SummarizedExperiment Input
#' # ====================================================================
#' \dontrun{
#'     library(SummarizedExperiment)
#'     
#'     # Extract phenotype and ensure rownames match matrix colnames
#'     pheno_data <- dnamExample[[2]]
#'     rownames(pheno_data) <- colnames(beta_matrix)
#'     
#'     # Construct the SummarizedExperiment object
#'     se_obj <- SummarizedExperiment(
#'       assays = list(beta = beta_matrix),
#'       colData = pheno_data
#'     )
#'     
#'     # The function seamlessly accepts the Bioconductor object
#'     predAge <- leeCommonABEC(x = se_obj)
#'     head(predAge)
#' }


leeCommonABEC <- function(x,
                    minCoverage = 0,
                    verbose = TRUE) {
  .runEpiClockPipeline(
    x = x, 
    coefName = "omniager_cabec_coef", 
    clockName = "leeCommonABEC",
    minCoverage = minCoverage, 
    verbose = verbose,
    useHorvathTrafo = FALSE
  )
}
