#' Pre-Processing of input dataset for scImmuAging
#'
#' @param inputObj A Seurat object. Meta data must contain "donorId" and 
#' "age" columns.
#' @param cellType A string specifying the cell type.
#' @param model The pre-trained model object.
#' @param markerGenes Character vector of selected marker genes.
#' @return A nested and unnested data frame with pseudocells.
#' @import Seurat dplyr purrr tidyr
#' @export
#' 

scImmuAgingPreProcess <- function(inputObj, cellType, model, markerGenes) {

  newSeurat <- subset(inputObj, features = markerGenes)
  
  preprocessedData <- scImmuAgingPreprocessing(newSeurat) %>%
    dplyr::group_by(donorId, age) %>%
    tidyr::nest()
  
  preprocessedData <- preprocessedData %>%
    dplyr::mutate(pseudocells = purrr::map(data, pseudocellScImmuAging))
  
  preprocessedData$data <- NULL
  preprocessedData <- tidyr::unnest(preprocessedData, pseudocells)
  
  return(preprocessedData)
}





#' Pre-Processing of input dataset
#'
#' input is a seurat object. meta.data should have "donor_id" and "age" columns.
#' @param seuratObj seurat object with metadata "donor_id", "age"
#' @return Return pre-processed input
#' @import Seurat dplyr glmnet purrr
#' @export

scImmuAgingPreprocessing <- function(seuratObj) {
  DefaultAssay(seuratObj) <- "RNA"
  metaData <- seuratObj@meta.data
  
  if (!"donorId" %in% colnames(metaData)) {
    if ("donor_id" %in% colnames(metaData)) {
      metaData$donorId <- metaData$donor_id
    } else {
      stop("Error: please add 'donorId' (or 'donor_id') in your metadata!")
    }
  }
  

  if (!"age" %in% colnames(metaData)) {
    stop("Error: please add 'age' in your metadata!")
  }
  
  metaDataSubset <- metaData[, c("donorId", "age"), drop = FALSE]
  
  inputMtx <- t(as.matrix(GetAssayData(seuratObj, assay = "RNA", layer = "data")))
  
  combinedInput <- dplyr::as_tibble(cbind(metaDataSubset, inputMtx))
  return(combinedInput)
}

#' Pre-Processing of input dataset
#'
#' input could be interested gene list. Each row is one gene and columns are 
#' "gene" and "value". "value" could be log2FC or specific value.
#' colnames of other metadata should be "subtype".
#' @param inputData the output from scImmuAgingPreProcess() function.
#' @param size how many cells would be randomly selected to generate pseudocells
#' @param n how many times to repeat random selection
#' @param replace Character or Logical. Strategy for sampling.
#'   If "dynamic" (default), it automatically sets replace to TRUE 
#'   if nrow(input) <= size. Otherwise, can be set to TRUE or FALSE explicitly.
#' @return Return ranked data.frame based on "value"
#' @import Seurat dplyr glmnet purrr
#' @export

pseudocellScImmuAging <- function(inputData, size = 15, n = 100, replace = "dynamic") {
  if (replace == "dynamic") {
    replace <- nrow(inputData) <= size
  }
  
  mat <- as.matrix(inputData)
  numRows <- nrow(mat)
  
  indices <- replicate(n, sample(seq_len(numRows), size = size, replace = replace))
  
  pseudoMat <- vapply(seq_len(n), function(i) {
    colMeans(mat[indices[, i], , drop = FALSE])
  }, numeric(ncol(mat)))
  
  return(dplyr::as_tibble(t(pseudoMat)))
}





#' predict age for each individual
#'
#' @param preprocessedData the output from scImmuAging_PreProcess() function.
#' @param model the cell type aging clock model
#' @param markerGenes selected marker genes from the corresponding model
#' @return Return dataframe of predicted age
#' @import Seurat dplyr glmnet purrr
#' @export
#' @noRd


scImmuAgingCalculator <- function(preprocessedData, model, markerGenes, 
                                  minCoverage = 0.5, verbose = TRUE) {

  testMat <- t(as.matrix(preprocessedData[, -c(1, 2)])) 
  

  fakeWeights <- rep(0, length(markerGenes))
  names(fakeWeights) <- markerGenes
  
  # 1. Check the coverage
  coverage <- .checkCpGCoverage(
    betaM = testMat, 
    allWeights = fakeWeights, 
    clockName = "scImmuAging", 
    minCoverage = minCoverage, 
    verbose = verbose
  )
  
  if (!coverage$pass) {
    return(data.frame(donorId = preprocessedData[[1]], 
                      age = preprocessedData[[2]], 
                      prediction = NA_real_))
  }
  
  # 2. Prepare the final matrix
  finalMat <- matrix(0, nrow = length(markerGenes), ncol = ncol(testMat))
  rownames(finalMat) <- markerGenes
  finalMat[rownames(testMat)[coverage$betaIdx], ] <- testMat[coverage$betaIdx, ]
  
  # 3. Prediction
  testPredictions <- predict(model, newx = t(finalMat), s = "lambda.min")
  
  return(data.frame(
    donorId = preprocessedData[[1]],
    age = preprocessedData[[2]],
    prediction = as.numeric(testPredictions)
  ))
}





#' predicted age for each individual
#'
#' @param predictRes the output from scImmuAgingCalculator() function
#' @return Return predicted age for each individual
#' @import Seurat dplyr glmnet purrr
#' @export


ageDonor <- function(predictRes) {
  donorDf <- predictRes %>%
    dplyr::group_by(donorId) %>%
    dplyr::mutate(predicted = round(mean(prediction, na.rm = TRUE))) %>%
    dplyr::select(-prediction) %>%
    dplyr::distinct()
  
  return(as.data.frame(donorDf))
}
