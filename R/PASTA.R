#' @title Calculate PASTA-based transcriptomic age scores
#'
#' @description
#' This function computes transcriptomic age acceleration scores based on the
#' models described by Salignon et al. (2025), including the primary PASTA
#' score, a standard regression (REG) score, and the CT46 score.
#'
#' @param mat An expression matrix (numeric) with **genes as rows** and
#'   **samples as columns**.
#' @param filterGenes Logical. If \code{TRUE} (default), the matrix is
#'   subsetted to retain only the genes utilized by the pre-trained models.
#' @param rankNorm Logical. If \code{TRUE} (default), applies a rank-based
#'   inverse normal transformation (rank-normalization) to the expression data.
#' @param reg Logical. If \code{TRUE} (default), computes the REG (regression)
#'   age score.
#' @param pasta Logical. If \code{TRUE} (default), computes the PASTA age score.
#' @param ct46 Logical. If \code{TRUE} (default), computes the CT46 age score.
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
#' \dontrun{
#' library(magrittr)
#' library(Seurat)
#' library(glmnet)
#' downloadOmniAgeRExample("seu_gabitto_2024_filtered")
#' loadOmniAgeRExample("seu_gabitto_2024_filtered")
#'
#' seu$age <- seu$development_stage %>%
#'   gsub("-year.*", "", .) %>%
#'   gsub("-", " ", .) %>%
#'   gsub("80 year old and over stage", "85", .)
#'
#' set.seed(42)
#' seuBulk <- makePseudobulksPasta(
#'   seu,
#'   poolBy = c("cell_type", "age"),
#'   chunkSize = 512,
#'   verbose = FALSE
#' )
#'
#' lognormMatrix <- GetAssayData(seuBulk, assay = "RNA", layer = "data")
#' lognormMatrix <- as.matrix(lognormMatrix)
#'
#' seuBulkMeta <- seuBulk[[c("chunk_size", "cell_type", "age")]]
#' seuBulkMeta$age <- as.numeric(seuBulkMeta$age)
#' pastaRes <- pastaScores(lognormMatrix, filterGenes = TRUE, rankNorm = TRUE)
#' }
#' 
#' 
makePseudobulksPasta <- function(seu, poolBy = c("cell_type", "age"), 
                                 chunkSize = 1000, verbose = TRUE) 


pastaScores <- function(mat, filterGenes = TRUE, rankNorm = TRUE, 
                        reg = TRUE, pasta = TRUE, ct46 = TRUE) {
  
  # 1. Load model data
  data("PASTA_Gene", envir = environment())
  pastaGenesModel <- PASTA_genes_model
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
  if (reg) resList[["REG"]] <- predictAgeScore(matT, modelType = "REG")
  if (pasta) resList[["PASTA"]] <- predictAgeScore(matT, modelType = "PASTA")
  if (ct46) resList[["CT46"]] <- predictAgeScore(matT, modelType = "CT46")
  
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
#' @noRd
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
#' @noRd
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
#' @return Numeric vector. Predicted age scores.
#' @export

predictAgeScore <- function(mat, modelType = "PASTA") {
  cvfit_REG <- cvfit_PASTA <- cvfit_C46 <- NULL
  beta_PASTA <- beta_C46 <- NULL
  
  # 1. Load data according to the type
  dataName <- switch(modelType,
                     "REG" = "cvfit_REG",
                     "PASTA" = "cvfit_PASTA",
                     "CT46" = "cvfit_C46",
                     stop("Invalid modelType. Choose PASTA, REG, or CT46.")
  )
  data(list = dataName, envir = environment())
  
  # Obtain the current model object
  curModel <- get(dataName)
  
  # 2. Prediction
  vAgeScores <- as.numeric(stats::predict(curModel, mat, 
                                          s = "lambda.min", 
                                          type = "link")[, 1])
  
  # 3. Scaled
  if (modelType == "PASTA") {
    utils::data("beta_PASTA", envir = environment(), package = "OmniAgeR")
    vAgeScores <- vAgeScores * beta_PASTA
  } else if (modelType == "CT46") {
    utils::data("beta_C46", envir = environment(), package = "OmniAgeR")
    vAgeScores <- vAgeScores * beta_C46
  }
  
  return(vAgeScores)
}


#' Create Pseudobulk Samples from Seurat Object
#'
#' Aggregates single-cell expression data into pseudobulk samples based on
#' user-defined metadata variables and a specified chunk size.
#'
#' @param seu A Seurat object.
#' @param poolBy A character vector of column names from `@meta.data`. Cells are
#'   grouped by unique combinations of these variables prior to chunking.
#'   Default: `c("cell_type", "age")`.
#' @param chunkSize A numeric value. The target number of cells per pseudobulk
#'   sample. If set to 1, no aggregation is performed and the original
#'   object is returned (with metadata updated). Default: `1000`.
#' @param verbose Logical. If `TRUE`, prints a summary table detailing the
#'   number of pseudobulk samples generated per group. Default: `TRUE`.
#'
#' @return A new Seurat object where columns represent pseudobulk samples.
#'   The meta.data slot includes the original `poolBy` variables and the
#'   `chunkSize` used for aggregation.
#' @export
#
makePseudobulksPasta <- function(seu, poolBy = c("cell_type", "age"), 
                                 chunkSize = 1000, verbose = TRUE) {
  
  if (!all(poolBy %in% colnames(seu@meta.data))) {
    stop("Required metadata columns not found in Seurat object.")
  }
  
  if (chunkSize == 1) return(seu)
  
  # 1. Generate grouping factors
  groupFactor <- interaction(seu@meta.data[, poolBy], drop = TRUE)
  
  # 2. Block indexing calculation
  indices <- seq_len(nrow(seu@meta.data))
  chunkIds <- character(length(indices))
  
  splitIndices <- split(indices, groupFactor)
  
  for (groupName in names(splitIndices)) {
    idx <- splitIndices[[groupName]]
    n <- length(idx)
    numChunks <- ceiling(n / chunkSize)
    groupChunks <- rep(seq_len(numChunks), each = chunkSize, length.out = n)
    chunkIds[idx] <- paste(groupName, sample(groupChunks), sep = "-")
  }
  
  seu$tempChunkId <- make.names(chunkIds)
  
  # 3. Aggregation
  seuBulk <- Seurat::AggregateExpression(
    seu, 
    group.by = "tempChunkId", 
    return.seurat = TRUE, 
    verbose = FALSE
  )
  
  # 4. Metadata recovery 
  firstMatchIdx <- match(colnames(seuBulk), seu$tempChunkId)
  originalMeta <- seu@meta.data[firstMatchIdx, poolBy, drop = FALSE]
  
  for (col in poolBy) {
    seuBulk[[col]] <- originalMeta[[col]]
  }
  
  seuBulk$chunkSize <- chunkSize
  
  if (verbose) {
    print(table(seuBulk[[poolBy[1]]]))
  }
  
  return(seuBulk)
}