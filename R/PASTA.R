#' @title Calculate PASTA-based transcriptomic age scores
#'
#' @description
#' This function computes transcriptomic age acceleration scores based on the
#' models described by Salignon et al. (2025), including the primary PASTA
#' score, a standard regression (REG) score, and the CT46 score.
#' @param x A numeric matrix, \code{data.frame}, or \code{SummarizedExperiment} 
#'  object containing an expression matrix (numeric) with **genes as rows** and
#'   **samples as columns**
#' @param filterGenes Logical. If \code{TRUE} (default), the matrix is
#'   subsetted to retain only the genes utilized by the pre-trained models.
#' @param rankNorm Logical. If \code{TRUE} (default), applies a rank-based
#'   inverse normal transformation (rank-normalization) to the expression data.
#' @param reg Logical. If \code{TRUE}, computes the REG (regression)
#'   age score.
#' @param pasta Logical. If \code{TRUE} (default), computes the PASTA age score.
#' @param ct46 Logical. If \code{TRUE}, computes the CT46 age score.
#' @param verbose Logical. Whether to print status messages.
#'   Default is \code{TRUE}.
#'
#' @return
#' A \code{list} where each element is a numeric vector of predicted age scores
#' for the requested model(s) (e.g., \code{res_list$PASTA}, \code{res_list$REG}).
#'
#' @export
#'
#' @references
#' Salignon, J., Tsiokou, M., Marqués, P. et al.
#' Pasta, an age-shift transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation.
#' \emph{bioRxiv.} 2025
#'
#' @examples
#' library(glmnet)
#' data("ScPbmcExample")
#' expr_mat <- as.matrix(SummarizedExperiment::assay(ScPbmcExample, "logcounts"))
#' cell_meta <- as.data.frame(SummarizedExperiment::colData(ScPbmcExample))
#' groups <- paste(cell_meta$donor_id, cell_meta$celltype, sep = "_")
#' unique_groups <- unique(groups)
#'  pseudobulk_data_list <- lapply(unique_groups, function(grp) {
#'  rowSums(expr_mat[, groups == grp, drop = FALSE])
#' })
#' 
#' pseudobulk_data_matrix <- do.call(cbind, pseudobulk_data_list)
#' colnames(pseudobulk_data_matrix) <- unique_groups
#' 
#' pastaRes <- pastaScores(pseudobulk_data_matrix, filterGenes = TRUE, 
#'                         rankNorm = TRUE)
#'
#'
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/jsalignon/pasta
# under the MIT License.
# Modifications: Added generic object support (SummarizedExperiment), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------

pastaScores <- function(x, filterGenes = TRUE, rankNorm = TRUE,
                        reg = FALSE, pasta = TRUE, ct46 = FALSE, verbose = TRUE) {
    mat <- .extractAssayMatrix(x)
    # 1. Load model data
    pastaGenesModel <- loadOmniAgeRdata(
        "omniager_pasta_gene",
        verbose = verbose
    )

    resList <- list()

    # 2. Preprocessing
    if (filterGenes) {
        mat <- filterAgeModelGenes(mat, pastaGenesModel)
    }

    if (rankNorm) {
        mat <- applyRankNormalization(mat)
    }

    matT <- t(mat)

    # 3. Prediction
    if (reg) resList[["REG"]] <- .predictAgeScore(matT, modelType = "REG")
    if (pasta) resList[["PASTA"]] <- .predictAgeScore(matT, modelType = "PASTA")
    if (ct46) resList[["CT46"]] <- .predictAgeScore(matT, modelType = "CT46")

    return(resList)
}


#' Filter Age Model Genes from Count Matrix
#'
#' Subsets the count matrix to include only genes used in the age prediction model.
#'
#' @param mat Matrix. Count matrix.
#' @param pastaGenesModel The genes used for building the model.
#' @return Matrix. Filtered count matrix with median imputation.
#' @export
#' @examples
#' # 1. Create a mock count matrix (5 genes, 3 samples)
#' mock_mat <- matrix(1:15, nrow = 5, ncol = 3)
#' rownames(mock_mat) <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
#' colnames(mock_mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # 2. Define the genes required by the model
#' # Note: "GeneF" is intentionally missing from the matrix to test imputation
#' model_genes <- c("GeneA", "GeneC", "GeneF")
#'
#' # 3. Run the filter function
#' filtered_mat <- filterAgeModelGenes(mat = mock_mat, pastaGenesModel = model_genes)
#'
#' # 4. View the result
#' print(filtered_mat)
filterAgeModelGenes <- function(mat, pastaGenesModel) {
    if (is.null(rownames(mat))) {
        stop("Matrix must have row names (gene symbols).")
    }

    idx <- match(pastaGenesModel, rownames(mat))
    matFiltered <- mat[idx, , drop = FALSE]

    medianVal <- stats::median(matFiltered, na.rm = TRUE)
    matFiltered[is.na(matFiltered)] <- medianVal

    rownames(matFiltered) <- pastaGenesModel
    return(matFiltered)
}

#' Apply Rank Normalization to Matrix
#'
#' Applies rank normalization across each column of the matrix.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Rank-normalized matrix.
#' @export
#' @examples
#' # 1. Create a mock count matrix (3 genes, 3 samples)
#' # Notice that Sample1 has a tie (10 and 10), and Sample3 has identical values
#' mock_mat <- matrix(
#'     c(
#'         10, 20, 10, # Sample 1
#'         5, 50, 15, # Sample 2
#'         8, 8, 8
#'     ), # Sample 3
#'     nrow = 3, ncol = 3
#' )
#' rownames(mock_mat) <- c("GeneA", "GeneB", "GeneC")
#' colnames(mock_mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # 2. View the original matrix
#' print(mock_mat)
#'
#' # 3. Apply rank normalization
#' norm_mat <- applyRankNormalization(mat = mock_mat)
#'
#' # 4. View the rank-normalized matrix
#' # In Sample1, the two 10s will share the average rank of 1.5
#' print(norm_mat)
applyRankNormalization <- function(mat) {
    matNormalized <- apply(mat, 2, function(x) {
        r <- rank(x, ties.method = "average")
        return(r)
    })
    return(matNormalized)
}

#' Predict Age Score from Gene Expression Matrix
#'
#' Uses pre-trained models to predict age scores based on gene expression.
#'
#' @param mat Matrix. Processed count matrix.
#' @param modelType Character. Model type ('PASTA', 'REG', or 'CT46').
#' @param verbose Logical. Whether to print status messages during data loading. Default is FALSE.
#' @return Numeric vector. Predicted age scores.
#' @noRd

.predictAgeScore <- function(mat, modelType = "PASTA", verbose = FALSE) {
    # 1. Load data according to the type
    dataName <- switch(modelType,
        "REG" = "omniager_cvfit_reg",
        "PASTA" = "omniager_cvfit_pasta",
        "CT46" = "omniager_cvfit_c46",
        stop("Invalid modelType. Choose PASTA, REG, or CT46.")
    )
    # Obtain the current model object
    curModel <- loadOmniAgeRdata(
        dataName,
        verbose = verbose
    )

    # 2. Prediction
    vAgeScores <- as.numeric(stats::predict(curModel, mat,
        s = "lambda.min",
        type = "link"
    )[, 1])

    # 3. Scaled
    if (modelType == "PASTA") {
        betaPASTA <- loadOmniAgeRdata(
            "omniager_beta_pasta",
            verbose = verbose
        )
        vAgeScores <- vAgeScores * betaPASTA
    } else if (modelType == "CT46") {
        betaC46 <- loadOmniAgeRdata(
            "omniager_beta_c46",
            verbose = verbose
        )
        vAgeScores <- vAgeScores * betaC46
    }

    return(vAgeScores)
}

