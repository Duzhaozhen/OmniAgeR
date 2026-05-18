#' @title Compute a DNAm-based Smoking Index
#'
#' @description
#' Computes a DNA methylation-based smoking index (S-index) using 1,501 smoking-associated
#' CpG sites. This index has been shown to correlate with cancer risk across multiple
#' tissues and provides a robust epigenetic signature of smoking exposure.
#'
#' @details
#' The function calculates the score by first standardizing the methylation data (z-score
#' transformation) across samples for each CpG site, and then computing a weighted
#' average based on the provided signature coefficients.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named \code{numeric} vector of smoking index for each sample
#' in \code{data.m}.
#'
#'
#' @importFrom stats sd
#' @export
#'
#' @references
#' Teschendorff, Andrew E et al.
#' Correlation of Smoking-Associated DNA Methylation Changes in Buccal Cells
#' With DNA Methylation Changes in Epithelial Cancer
#' \emph{JAMA oncology} 2015
#'
#'
#' @examples
#' # 1. Load the lightweight clock coefficient table
#' modelCoef <- loadOmniAgeRdata("omniager_coeff_smk_idx", verbose = FALSE)
#' 
#' # 2. Extract feature names and exclude potential intercept terms
#' allFeatures <- names(modelCoef)
#' requiredCpGs <- allFeatures[grep("^cg", allFeatures)]
#' if (length(requiredCpGs) == 0) requiredCpGs <- allFeatures 
#' 
#' # 3. Generate a mock micro-beta matrix for 2 samples in memory
#' mockBetaM <- matrix(
#'     runif(length(requiredCpGs) * 2, min = 0, max = 1),
#'     nrow = length(requiredCpGs),
#'     dimnames = list(requiredCpGs, c("Sample1", "Sample2"))
#' )
#' # 4. Run the age prediction
#' SmokeIndexOut <- compSmokeIndex(mockBetaM)

compSmokeIndex <- function(betaM, minCoverage = 0, verbose = TRUE) {
    coeffSmkIdx <- loadOmniAgeRdata(
        "omniager_coeff_smk_idx",
        verbose = verbose
    )

    # Perform coverage check
    coverage <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = coeffSmkIdx,
        clockName = "SmokeIndex",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverage$pass) {
        scores <- rep(NA_real_, ncol(betaM))
        names(scores) <- colnames(betaM)
        return(scores)
    }


    subBeta <- betaM[coverage$betaIdx, , drop = FALSE]
    weights <- coverage$weightsSubset

    # Z-score
    rowMeansV <- rowMeans(subBeta, na.rm = TRUE)
    rowSdsV <- apply(subBeta, 1, sd, na.rm = TRUE)

    keepIdx <- which(rowSdsV > 0)
    if (length(keepIdx) < length(rowSdsV)) {
        if (verbose) {
            message(sprintf(
                "[SmokeIndex] Excluding %d constant probes.",
                length(rowSdsV) - length(keepIdx)
            ))
        }
        subBeta <- subBeta[keepIdx, , drop = FALSE]
        rowMeansV <- rowMeansV[keepIdx]
        rowSdsV <- rowSdsV[keepIdx]
        weights <- weights[keepIdx]
    }

    zMatrix <- (subBeta - rowMeansV) / rowSdsV

    # Calculate the final score
    scores <- colMeans(weights * zMatrix, na.rm = TRUE)

    return(scores)
}
