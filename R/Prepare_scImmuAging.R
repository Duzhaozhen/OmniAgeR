#' Prepare single-cell expression data for statistical modeling
#'
#' @description
#' This internal function transforms standard single-cell expression matrices
#' (genes as rows) into a flattened, cell-based data frame (cells as rows).
#' It ensures proper alignment between expression data and metadata, 
#' effectively creating a "wide-format" data structure suitable for 
#' statistical modeling (e.g., linear regression analysis).
#'
#' @param expr A numeric \code{matrix} of gene expression values, where 
#'   rows represent genes and columns represent cells.
#' @param metadata A \code{data.frame} containing cell-level metadata. Must 
#'   contain at least the columns specified in \code{donorCol} and \code{ageCol}.
#' @param donorCol A character string specifying the column name in 
#'   \code{metadata} representing the donor identity. Defaults to \code{"donor_id"}.
#' @param ageCol A character string specifying the column name in 
#'   \code{metadata} representing the subject age. Defaults to \code{"age"}.
#'
#' @return A \code{data.frame} where each row represents a unique cell, 
#'   containing the metadata columns (donorId, age) followed by the 
#'   gene expression values as columns.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that the expression matrix is numeric.
#'   \item Aligns the rows of \code{metadata} with the columns of \code{expr} 
#'     using column names (cell identifiers).
#'   \item Transposes the expression matrix to switch from (genes x cells) 
#'     to (cells x genes).
#'   \item Merges metadata and expression data into a single \code{data.frame}.
#' }
#'
#' @keywords internal
#' @noRd
.scImmuAgingPreprocessingCore <- function(expr,
                                          metadata,
                                          donorCol = "donor_id",
                                          ageCol = "age") {
  
  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }
  
  if (!is.numeric(expr)) {
    stop("'expr' must contain numeric expression values.")
  }
  
  if (is.null(colnames(expr))) {
    stop("'expr' must have cell names as column names.")
  }
  
  # Ensure metadata rows are aligned with expression matrix columns
  if (!is.null(rownames(metadata)) &&
      all(colnames(expr) %in% rownames(metadata))) {
    metadata <- metadata[colnames(expr), , drop = FALSE]
  }
  
  if (!donorCol %in% colnames(metadata)) {
    stop("Metadata must contain donor column: ", donorCol)
  }
  
  if (!ageCol %in% colnames(metadata)) {
    stop("Metadata must contain age column: ", ageCol)
  }
  
  metaDataSubset <- data.frame(
    donorId = metadata[[donorCol]],
    age = metadata[[ageCol]],
    stringsAsFactors = FALSE
  )
  
  inputMtx <- t(expr)
  
  combinedInput <- data.frame(
    metaDataSubset,
    as.data.frame(inputMtx, check.names = FALSE),
    check.names = FALSE
  )
  
  return(combinedInput)
}

#' Generate pseudocells for downstream statistical analysis
#'
#' @description
#' This function creates "pseudocells" (also known as pseudo-bulk samples) 
#' by aggregating single-cell expression data of a specific cell type. 
#' It reduces technical noise and dropout effects by combining multiple 
#' single cells into averaged expression profiles, which are more suitable 
#' for linear modeling and differential aging analysis.
#'
#' @param expr A numeric \code{matrix} of gene expression (rows: genes, cols: cells).
#' @param metadata A \code{data.frame} of cell-level metadata.
#' @param cellType A character string identifying the cell type to be processed.
#' @param markerGenes A character vector of genes to include in the pseudocell generation.
#' @param donorCol Column name in metadata for donor IDs. Defaults to \code{"donor_id"}.
#' @param ageCol Column name in metadata for subject age. Defaults to \code{"age"}.
#' @param cellTypeCol Column name in metadata for cell types. Defaults to \code{"celltype"}.
#' @param pseudocellSize Integer. The number of single cells to aggregate per pseudocell.
#' @param pseudocellN Integer. The target number of pseudocells to generate.
#' @param replace Logical or "dynamic". Whether to sample cells with replacement.
#' @param verbose Logical. Whether to print progress and missing gene messages.
#'
#' @return A \code{data.frame} of aggregated expression values (pseudocells), 
#'   with metadata columns (donorId, age) and gene expression columns.
#'
#' @details
#' The function follows these steps:
#' \enumerate{
#'   \item Filters expression data for the specific \code{cellType} and \code{markerGenes}.
#'   \item Formats data using \code{.scImmuAgingPreprocessingCore}.
#'   \item Groups data by donor and age.
#'   \item Performs aggregation into pseudocells using the specified \code{pseudocellSize}.
#' }
#'
#' @importFrom stats aggregate
#' @export

.scImmuAgingMakePseudocells <- function(expr,
                                        metadata,
                                        cellType,
                                        markerGenes,
                                        donorCol = "donor_id",
                                        ageCol = "age",
                                        cellTypeCol = "celltype",
                                        pseudocellSize = 15,
                                        pseudocellN = 100,
                                        replace = "dynamic",
                                        verbose = TRUE) {
  if (!cellTypeCol %in% colnames(metadata)) {
    stop("Metadata must contain cell type column: ", cellTypeCol)
  }
  
  keepCells <- metadata[[cellTypeCol]] == cellType
  keepCells[is.na(keepCells)] <- FALSE
  
  if (!any(keepCells)) {
    stop("No cells found for cell type: ", cellType)
  }
  
  markerGenes <- unique(markerGenes)
  presentGenes <- markerGenes[markerGenes %in% rownames(expr)]
  
  if (length(presentGenes) == 0L) {
    stop("None of the marker genes were found in the input object.")
  }
  
  missingGenes <- setdiff(markerGenes, presentGenes)
  
  if (verbose && length(missingGenes) > 0L) {
    message(
      "[scImmuAging] ",
      length(missingGenes),
      " marker genes were not found for cell type ",
      cellType,
      "."
    )
  }
  
  exprSub <- expr[presentGenes, keepCells, drop = FALSE]
  metadataSub <- metadata[keepCells, , drop = FALSE]
  
  cellDf <- .scImmuAgingPreprocessingCore(
    expr = exprSub,
    metadata = metadataSub,
    donorCol = donorCol,
    ageCol = ageCol
  )
  
  groupKey <- paste(cellDf$donorId, cellDf$age, sep = "\r")
  groupKey <- factor(groupKey, levels = unique(groupKey))
  groupIdx <- split(seq_len(nrow(cellDf)), groupKey)
  
  outList <- lapply(groupIdx, function(idx) {
    donorId <- cellDf$donorId[idx[1]]
    age <- cellDf$age[idx[1]]
    
    geneDf <- cellDf[
      idx,
      setdiff(colnames(cellDf), c("donorId", "age")),
      drop = FALSE
    ]
    
    pseudo <- .pseudocellScImmuAging(
      inputData = geneDf,
      size = pseudocellSize,
      n = pseudocellN,
      replace = replace
    )
    
    data.frame(
      donorId = rep(donorId, nrow(pseudo)),
      age = rep(age, nrow(pseudo)),
      pseudo,
      check.names = FALSE
    )
  })
  
  preprocessedData <- do.call(rbind, outList)
  rownames(preprocessedData) <- NULL
  
  return(preprocessedData)
}



## -------------------------------------------------------------------------
## Exported function: pseudocell generation
## -------------------------------------------------------------------------

#' Generate pseudocells for scImmuAging
#'
#' @param inputData A numeric matrix or data.frame with cells in rows and genes
#' in columns.
#' @param size Number of cells sampled to generate each pseudocell.
#' @param n Number of pseudocells to generate.
#' @param replace If \code{"dynamic"}, sampling with replacement is used when
#' the number of available cells is less than or equal to \code{size}. Otherwise
#' use \code{TRUE} or \code{FALSE}.
#'
#' @return A data.frame of pseudocell expression values.
#' @noRd
#' @keywords internal
.pseudocellScImmuAging <- function(inputData,
                                   size = 15,
                                   n = 100,
                                   replace = "dynamic") {
  mat <- as.matrix(inputData)
  
  if (!is.numeric(mat)) {
    stop("'inputData' must contain numeric values.")
  }
  
  numRows <- nrow(mat)
  
  if (numRows == 0L) {
    stop("'inputData' must contain at least one cell.")
  }
  
  if (!is.numeric(size) || length(size) != 1L || size <= 0L) {
    stop("'size' must be a positive numeric value.")
  }
  
  if (!is.numeric(n) || length(n) != 1L || n <= 0L) {
    stop("'n' must be a positive numeric value.")
  }
  
  size <- as.integer(size)
  n <- as.integer(n)
  
  if (identical(replace, "dynamic")) {
    replace <- numRows <= size
  }
  
  if (!is.logical(replace) || length(replace) != 1L) {
    stop("'replace' must be TRUE, FALSE, or 'dynamic'.")
  }
  
  if (!replace && numRows < size) {
    stop(
      "Cannot sample ", size, " cells without replacement from only ",
      numRows, " cells. Use replace = TRUE or replace = 'dynamic'."
    )
  }
  
  indices <- replicate(
    n,
    sample(seq_len(numRows), size = size, replace = replace)
  )
  
  pseudoMat <- vapply(seq_len(n), function(i) {
    colMeans(mat[indices[, i], , drop = FALSE])
  }, numeric(ncol(mat)))
  
  pseudoDf <- as.data.frame(t(pseudoMat), check.names = FALSE)
  colnames(pseudoDf) <- colnames(mat)
  
  return(pseudoDf)
}


## -------------------------------------------------------------------------
## Exported function: prediction calculator
## -------------------------------------------------------------------------

#' Predict age for each pseudocell
#'
#' @param preprocessedData Output from \code{scImmuAgingPreProcess()}.
#' @param model The cell type-specific aging clock model.
#' @param markerGenes Character vector of marker genes used by the model.
#' @param minCoverage Minimum required proportion of marker genes present.
#' Default is 0.5.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A data.frame with columns \code{donorId}, \code{age}, and
#' \code{prediction}.
#' @noRd
#' @keywords internal
.scImmuAgingCalculator <- function(preprocessedData,
                                  model,
                                  markerGenes,
                                  minCoverage = 0.5,
                                  verbose = TRUE) {
  requiredCols <- c("donorId", "age")
  missingCols <- setdiff(requiredCols, colnames(preprocessedData))
  
  if (length(missingCols) > 0L) {
    stop(
      "'preprocessedData' is missing required columns: ",
      paste(missingCols, collapse = ", ")
    )
  }
  
  geneCols <- setdiff(colnames(preprocessedData), requiredCols)
  
  if (length(geneCols) == 0L) {
    stop("'preprocessedData' does not contain gene expression columns.")
  }
  
  testMat <- t(as.matrix(preprocessedData[, geneCols, drop = FALSE]))
  
  if (!is.numeric(testMat)) {
    stop("Gene expression columns in 'preprocessedData' must be numeric.")
  }
  
  fakeWeights <- rep(0, length(markerGenes))
  names(fakeWeights) <- markerGenes
  
  coverage <- .checkCpGCoverage(
    betaM = testMat,
    allWeights = fakeWeights,
    clockName = "scImmuAging",
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  if (!coverage$pass) {
    return(data.frame(
      donorId = preprocessedData$donorId,
      age = preprocessedData$age,
      prediction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  finalMat <- matrix(
    0,
    nrow = length(markerGenes),
    ncol = ncol(testMat),
    dimnames = list(markerGenes, colnames(testMat))
  )
  
  presentGenes <- rownames(testMat)[coverage$betaIdx]
  
  finalMat[presentGenes, ] <- testMat[
    coverage$betaIdx,
    ,
    drop = FALSE
  ]
  
  testPredictions <- stats::predict(
    model,
    newx = t(finalMat),
    s = "lambda.min"
  )
  
  return(data.frame(
    donorId = preprocessedData$donorId,
    age = preprocessedData$age,
    prediction = as.numeric(testPredictions),
    stringsAsFactors = FALSE
  ))
}




## -------------------------------------------------------------------------
## Exported function: aggregate prediction to donor level
## -------------------------------------------------------------------------

#' Aggregate scImmuAging predictions by donor
#'
#' @param predictRes Output from \code{.scImmuAgingCalculator()}.
#'
#' @return A data.frame with one row per donor and columns \code{donorId},
#' \code{age}, and \code{predicted}.
#' @noRd
#' @keywords internal
.ageDonor <- function(predictRes) {
  if (!"donorId" %in% colnames(predictRes)) {
    stop("'predictRes' must contain a 'donorId' column.")
  }
  
  if (!"prediction" %in% colnames(predictRes)) {
    stop("'predictRes' must contain a 'prediction' column.")
  }
  
  donorIds <- unique(predictRes$donorId)
  
  out <- lapply(donorIds, function(id) {
    tmp <- predictRes[predictRes$donorId == id, , drop = FALSE]
    
    predicted <- if (all(is.na(tmp$prediction))) {
      NA_real_
    } else {
      mean(tmp$prediction, na.rm = TRUE)
    }
    
    data.frame(
      donorId = id,
      age = if ("age" %in% colnames(tmp)) tmp$age[1] else NA,
      predicted = predicted,
      stringsAsFactors = FALSE
    )
  })
  
  donorDf <- do.call(rbind, out)
  rownames(donorDf) <- NULL
  
  return(donorDf)
}





