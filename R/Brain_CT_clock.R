#' Run Brain Cell Type Specific Clocks
#'
#' @description
#' Run brain cell type-specific aging clock prediction for one or more cell
#' types and one or more sample types.
#'
#' @param x A \code{SingleCellExperiment}, \code{SummarizedExperiment}, Seurat
#' object, matrix, sparse Matrix, or data.frame containing log-normalized
#' expression values. Rows should be genes and columns should be cells.
#' For Bioconductor workflows, \code{SingleCellExperiment} input is recommended.
#' @param cellTypes Character vector of cell types to analyze. Available cell
#' types include \code{"Oligodendrocytes"}, \code{"Astrocytes"},
#' \code{"Microglia"}, \code{"OPCs"}, \code{"Excitatory Neurons"}, and
#' \code{"Inhibitory Neurons"} when corresponding models are available.
#' @param modelNames Character vector specifying which sample types to run.
#' Use \code{"all"} to run \code{"SC"}, \code{"Pseudobulk"}, and
#' \code{"Bootstrap"}.
#' @param metadata Optional metadata. Required when \code{x} is matrix-like.
#' @param assayName Assay name for \code{SingleCellExperiment} or
#' \code{SummarizedExperiment} input.
#' @param donorCol Column in metadata containing donor IDs.
#' @param ageCol Column in metadata containing donor ages.
#' @param cellTypeCol Column in metadata containing cell type labels.
#' @param verbose Logical. Whether to print status messages.
#' @param seuratAssay Assay name for Seurat input.
#' @param seuratLayer Layer or slot name for Seurat input.
#'
#' @return A named list where each element corresponds to a sample type that
#' was run. Each element contains a data.frame of predictions.
#'
#' @details
#' This function supports Bioconductor-native workflows by accepting
#' \code{SingleCellExperiment} and \code{SummarizedExperiment} objects directly.
#' Expression values are extracted from an assay, and cell-level metadata are
#' extracted from \code{colData()}.
#'
#' Seurat input is also supported when the Seurat package is installed, but
#' Seurat is not required for Bioconductor workflows.
#'
#' @references
#' Muralidharan C, Zakar-Polyák E, Adami A, et al.
#' Human Brain Cell-Type-Specific Aging Clocks
#' Based on Single-Nuclei Transcriptomics.
#' \emph{Adv Sci(Weinh).} 2025
#' @export
#' @examples
#' data("ScBrainExample")
#'
#' # Recommended Bioconductor workflow: SingleCellExperiment input
#' brain_res <- brainCtClock(
#'     x = ScBrainExample,
#'     cellTypes = "Oligodendrocytes",
#'     modelNames = "SC",
#'     assayName = "logcounts",
#'     donorCol = "donor_id",
#'     ageCol = "age",
#'     cellTypeCol = "celltype",
#'     verbose = FALSE
#' )
#'
#' names(brain_res)
#' head(brain_res$SC)
#' 
#' \dontrun{
#' # Matrix input is also supported when metadata are supplied separately
#' expr_mat <- SummarizedExperiment::assay(ScBrainExample, "logcounts")
#' cell_meta <- as.data.frame(SummarizedExperiment::colData(ScBrainExample))
#'
#' matrix_res <- brainCtClock(
#'     x = expr_mat,
#'     metadata = cell_meta,
#'     cellTypes = "Oligodendrocytes",
#'     modelNames = "SC",
#'     donorCol = "donor_id",
#'     ageCol = "age",
#'     cellTypeCol = "celltype",
#'     verbose = FALSE
#' )
#'
#' head(matrix_res$SC)
#' # Optional Seurat input
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     expr_mat <- SummarizedExperiment::assay(ScBrainExample, "logcounts")
#'     cell_meta <- as.data.frame(SummarizedExperiment::colData(ScBrainExample))
#'     cell_meta <- cell_meta[colnames(expr_mat), , drop = FALSE]
#'
#'     seurat_obj <- Seurat::CreateSeuratObject(
#'         counts = expr_mat,
#'         meta.data = cell_meta
#'     )
#'
#'     seurat_obj <- tryCatch(
#'         Seurat::SetAssayData(
#'             object = seurat_obj,
#'             assay = "RNA",
#'             layer = "data",
#'             new.data = expr_mat
#'         ),
#'         error = function(e) {
#'             Seurat::SetAssayData(
#'                 object = seurat_obj,
#'                 assay = "RNA",
#'                 slot = "data",
#'                 new.data = expr_mat
#'             )
#'         }
#'     )
#'
#'     seurat_res <- brainCtClock(
#'         x = seurat_obj,
#'         cellTypes = "Oligodendrocytes",
#'         modelNames = "SC",
#'         seuratAssay = "RNA",
#'         seuratLayer = "data",
#'         donorCol = "donor_id",
#'         ageCol = "age",
#'         cellTypeCol = "celltype",
#'         verbose = FALSE
#'     )
#'
#'     head(seurat_res$SC)
#' }
#' }
#' 

brainCtClock <- function(x,
                         cellTypes,
                         modelNames = "all",
                         metadata = NULL,
                         assayName = "logcounts",
                         donorCol = "donor_id",
                         ageCol = "age",
                         cellTypeCol = "celltype",
                         verbose = TRUE,
                         seuratAssay = "RNA",
                         seuratLayer = "data") {
  validModels <- c("SC", "Pseudobulk", "Bootstrap")
  
  if (length(modelNames) == 1L && identical(modelNames, "all")) {
    modelsToRun <- validModels
  } else {
    modelsToRun <- modelNames
    
    if (!all(modelsToRun %in% validModels)) {
      stop(
        "Invalid 'modelNames'. Valid values are 'SC', ",
        "'Pseudobulk', 'Bootstrap', or 'all'."
      )
    }
  }
  
  if (verbose) {
    message("[brainCtClock] Starting brain cell type clock prediction...")
  }
  
  allResults <- list()
  
  for (currentModelType in modelsToRun) {
    predResult <- .runPredictionPipelineBrainCt(
      sampleType = currentModelType,
      x = x,
      cellTypes = cellTypes,
      metadata = metadata,
      assayName = assayName,
      donorCol = donorCol,
      ageCol = ageCol,
      cellTypeCol = cellTypeCol,
      verbose = verbose,
      seuratAssay = seuratAssay,
      seuratLayer = seuratLayer
    )
    
    allResults[[currentModelType]] <- predResult
  }
  
  if (verbose) {
    message("[brainCtClock] All brain clock predictions complete.")
  }
  
  allResults
}



#' Extract and pre-process data for brain cell type clocks
#'
#' @description
#' Subsets a single-cell expression object by cell type and processes the
#' expression data into one of three formats: single-cell, pseudobulk, or
#' bootstrap-resampled pseudobulk.
#'
#' @param x A \code{SingleCellExperiment}, \code{SummarizedExperiment}, Seurat
#' object, matrix, sparse Matrix, or data.frame containing log-normalized
#' expression values. Rows should be genes and columns should be cells.
#' For Bioconductor workflows, \code{SingleCellExperiment} input is recommended.
#' @param cellType Character. Cell type to analyze.
#' @param sampleType Character. One of \code{"SC"}, \code{"Pseudobulk"}, or
#' \code{"Bootstrap"}.
#' @param metadata Optional cell-level metadata. Required when \code{x} is a
#' matrix, sparse Matrix, or data.frame.
#' @param assayName Assay name for \code{SingleCellExperiment} or
#' \code{SummarizedExperiment} input. Default is \code{"logcounts"}.
#' @param donorCol Column in metadata containing donor IDs.
#' @param ageCol Column in metadata containing donor ages.
#' @param cellTypeCol Column in metadata containing cell type labels.
#' @param seuratAssay Assay name for Seurat input.
#' @param seuratLayer Layer or slot name for Seurat input.
#' @param bootstrapReps Number of bootstrap replicates per donor.
#' @param featureNames Optional character vector of model features to retain.
#' If provided, only these features are extracted before converting expression
#' values to a data.frame.
#' 
#' @return A data.frame containing processed expression data and metadata.
#' @keywords internal
#' @noRd

.getDfBrainCt <- function(x,
                          cellType,
                          sampleType,
                          metadata = NULL,
                          assayName = "logcounts",
                          donorCol = "donor_id",
                          ageCol = "age",
                          cellTypeCol = "celltype",
                          seuratAssay = "RNA",
                          seuratLayer = "data",
                          bootstrapReps = 100,
                          featureNames = NULL) {
  validSampleTypes <- c("SC", "Pseudobulk", "Bootstrap")
  
  if (!sampleType %in% validSampleTypes) {
    stop(
      "'sampleType' must be one of: ",
      paste(validSampleTypes, collapse = ", ")
    )
  }
  
  if (!is.numeric(bootstrapReps) ||
      length(bootstrapReps) != 1L ||
      bootstrapReps <= 0L) {
    stop("'bootstrapReps' must be a positive numeric value.")
  }
  
  bootstrapReps <- as.integer(bootstrapReps)
  
  input <- .extractSingleCellAssay(
    x = x,
    metadata = metadata,
    assayName = assayName,
    donorCol = donorCol,
    ageCol = ageCol,
    cellTypeCol = cellTypeCol,
    seuratAssay = seuratAssay,
    seuratLayer = seuratLayer
  )
  expr <- input$expr
  metadata <- input$metadata
  
  ## Keep only model features before dense conversion.
  ## This is important for sparse SingleCellExperiment/Seurat assays.
  if (!is.null(featureNames)) {
    featureNames <- unique(featureNames)
    presentFeatures <- intersect(featureNames, rownames(expr))
    
    if (length(presentFeatures) == 0L) {
      warning("No model features found in input for cell type: ", cellType)
      return(data.frame())
    }
    
    expr <- expr[presentFeatures, , drop = FALSE]
  }
  
  numCellsMap <- list(
    "Oligodendrocytes" = 200,
    "Astrocytes" = 50,
    "Microglia" = 50,
    "OPCs" = 50,
    "Excitatory Neurons" = 100,
    "Inhibitory Neurons" = 100
  )
  
  keepCells <- metadata[[cellTypeCol]] == cellType
  keepCells[is.na(keepCells)] <- FALSE
  
  if (!any(keepCells)) {
    warning("No cells found for: ", cellType)
    return(data.frame())
  }
  
  exprSub <- expr[, keepCells, drop = FALSE]
  metaSub <- metadata[keepCells, , drop = FALSE]
  
  metaOut <- data.frame(
    donorId = metaSub[[donorCol]],
    age = metaSub[[ageCol]],
    celltype = metaSub[[cellTypeCol]],
    stringsAsFactors = FALSE
  )
  
  ## Convert only the retained feature-by-cell matrix to dense matrix,
  ## then transpose to cell-by-gene format.
  exprSub <- as.matrix(exprSub)
  
  if (!is.numeric(exprSub)) {
    stop("The selected expression values must be numeric.")
  }
  
  exprCellByGene <- t(exprSub)
  
  exprDf <- as.data.frame(
    exprCellByGene,
    check.names = FALSE
  )
  
  exprDf[] <- lapply(exprDf, as.numeric)
  
  combinedData <- data.frame(
    metaOut,
    exprDf,
    check.names = FALSE
  )
  
  geneCols <- colnames(exprDf)
  
  if (length(geneCols) == 0L) {
    warning("No expression features remained after preprocessing.")
    return(data.frame())
  }
  
  if (identical(sampleType, "SC")) {
    return(combinedData)
  }
  
  if (identical(sampleType, "Pseudobulk")) {
    groupKey <- paste(
      combinedData$donorId,
      combinedData$age,
      combinedData$celltype,
      sep = "\r"
    )
    
    groupKey <- factor(groupKey, levels = unique(groupKey))
    groupIdx <- split(seq_len(nrow(combinedData)), groupKey)
    
    out <- lapply(groupIdx, function(idx) {
      tmp <- combinedData[idx, , drop = FALSE]
      
      means <- colMeans(
        as.matrix(tmp[, geneCols, drop = FALSE]),
        na.rm = TRUE
      )
      
      data.frame(
        donorId = tmp$donorId[1],
        age = tmp$age[1],
        celltype = tmp$celltype[1],
        as.data.frame(t(means), check.names = FALSE),
        check.names = FALSE
      )
    })
    
    res <- do.call(rbind, out)
    rownames(res) <- NULL
    return(res)
  }
  
  if (identical(sampleType, "Bootstrap")) {
    donors <- unique(combinedData$donorId)
    
    nSample <- numCellsMap[[cellType]]
    if (is.null(nSample)) {
      nSample <- 50
    }
    
    bootstrapList <- lapply(donors, function(d) {
      dfDonor <- combinedData[
        combinedData$donorId == d,
        ,
        drop = FALSE
      ]
      
      numRows <- nrow(dfDonor)
      
      if (numRows == 0L) {
        return(NULL)
      }
      
      indices <- do.call(
        cbind,
        replicate(
          bootstrapReps,
          sample(
            seq_len(numRows),
            size = nSample,
            replace = numRows < nSample
          ),
          simplify = FALSE
        )
      )
      
      bootMat <- vapply(seq_len(bootstrapReps), function(i) {
        colMeans(
          as.matrix(dfDonor[indices[, i], geneCols, drop = FALSE]),
          na.rm = TRUE
        )
      }, numeric(length(geneCols)))
      
      dfBoot <- as.data.frame(t(bootMat), check.names = FALSE)
      colnames(dfBoot) <- geneCols
      
      data.frame(
        donorId = d,
        age = dfDonor$age[1],
        celltype = cellType,
        dfBoot,
        check.names = FALSE
      )
    })
    
    bootstrapList <- Filter(Negate(is.null), bootstrapList)
    
    if (length(bootstrapList) == 0L) {
      return(data.frame())
    }
    
    res <- do.call(rbind, bootstrapList)
    rownames(res) <- NULL
    return(res)
  }
}



# --- Run Prediction Flow ---

#' Apply a Single Clock Model to Expression Data (Vectorized)
#'
#' @description
#' Predicts age using a pre-trained elastic net model (coefficients and intercept)
#' on a given expression data matrix. It handles gene matching, imputation for
#' missing genes, and vectorized prediction.
#'
#' @param inputData A data.frame of expression data (samples in rows, genes in
#'   columns). Must also contain 'age' and 'donor_id' columns.
#' @param imputeData A long-format data.frame with 'feature_name' and
#'   'imputation_value' columns. Used to fill in genes present in the model
#'   but missing from `data`.
#' @param modelObj A long-format data.frame with 'feature_name' (genes + 'intercept')
#'   and 'coefficient' columns, representing one trained clock.
#' @param sampleType A character string (e.g., 'SC', 'Pseudobulk') used to
#'   tag the output data.frame.
#'
#' @return
#' A data.frame with 'predictions' (the predicted age), 'ages' (the true age),
#' 'donors', and 'sample_type'.
#'
#' @details
#' This function is the core prediction engine. It performs a matrix
#' multiplication (`expression_matrix %*% coefficients_vector + intercept`).
#' It ensures that the gene order in the expression matrix exactly matches the
#' coefficient order from the model.
#' @noRd
#' @keywords internal

.predictBrainCtAge <- function(inputData, imputeData, modelObj, sampleType) {
  intercept <- modelObj$coefficient[modelObj$feature_name == "intercept"]
  
  if (length(intercept) == 0L) {
    intercept <- 0
  } else {
    intercept <- intercept[1]
  }
  
  modelGenesDf <- modelObj[modelObj$feature_name != "intercept", , drop = FALSE]
  if (nrow(modelGenesDf) == 0L) {
    return(data.frame(
      prediction = rep(intercept, nrow(inputData)),
      age = inputData$age,
      donorId = inputData$donorId,
      sampleType = sampleType,
      stringsAsFactors = FALSE
    ))
  }
  
  modelGenes <- modelGenesDf$feature_name
  
  presentGenes <- intersect(modelGenes, colnames(inputData))
  missingGenes <- setdiff(modelGenes, presentGenes)
  
  exprData <- inputData[, presentGenes, drop = FALSE]
  exprData <- as.data.frame(exprData, check.names = FALSE)
  
  exprData[] <- lapply(exprData, function(z) {
    if (is.factor(z)) {
      z <- as.character(z)
    }
    as.numeric(z)
  })
  
  exprMat <- as.matrix(exprData)
  
  if (!is.numeric(exprMat)) {
    stop("Expression columns in 'inputData' must be numeric.")
  }
  
  if (anyNA(exprMat)) {
    stop(
      "Some expression columns could not be converted to numeric values. ",
      "Please check whether non-expression columns were included among ",
      "model feature names."
    )
  }
  
  if (length(missingGenes) > 0L) {
    fillValues <- setNames(rep(0, length(missingGenes)), missingGenes)
    
    if (!is.null(imputeData) && nrow(imputeData) > 0L) {
      matchIdx <- match(missingGenes, imputeData$feature_name)
      matched <- !is.na(matchIdx)
      
      fillValues[matched] <- imputeData$imputation_value[
        matchIdx[matched]
      ]
    }
    
    fillMat <- matrix(
      rep(fillValues, each = nrow(exprMat)),
      nrow = nrow(exprMat),
      dimnames = list(NULL, missingGenes)
    )
    
    fullMat <- cbind(exprMat, fillMat)
  } else {
    fullMat <- exprMat
  }
  
  fullMat <- fullMat[, modelGenes, drop = FALSE]
  
  coefVector <- modelGenesDf$coefficient[
    match(modelGenes, modelGenesDf$feature_name)
  ]
  
  predVec <- as.numeric(fullMat %*% coefVector + intercept)
  
  data.frame(
    prediction = predVec,
    age = inputData$age,
    donorId = inputData$donorId,
    sampleType = sampleType,
    stringsAsFactors = FALSE
  )
}




#' Run the Full 5-Fold Averaged Prediction Pipeline
#'
#' @description
#' Orchestrates the brain cell type-specific prediction process for a given
#' sample type and set of cell types.
#'
#' @param sampleType Character. One of \code{"SC"}, \code{"Pseudobulk"}, or
#' \code{"Bootstrap"}.
#' @param x A \code{SingleCellExperiment}, \code{SummarizedExperiment}, Seurat
#' object, matrix, sparse Matrix, or data.frame containing log-normalized
#' expression values.
#' @param cellTypes Character vector of cell types to process.
#' @param metadata Optional metadata. Required when \code{x} is matrix-like.
#' @param assayName Assay name for \code{SingleCellExperiment} or
#' \code{SummarizedExperiment} input.
#' @param donorCol Column in metadata containing donor IDs.
#' @param ageCol Column in metadata containing donor ages.
#' @param cellTypeCol Column in metadata containing cell type labels.
#' @param verbose Logical. Whether to print status messages.
#' @param seuratAssay Assay name for Seurat input.
#' @param seuratLayer Layer or slot name for Seurat input.
#'
#' @return A data.frame containing 5-fold averaged predictions.
#' @keywords internal
#' @noRd


.runPredictionPipelineBrainCt <- function(sampleType,
                                         x,
                                         cellTypes,
                                         metadata = NULL,
                                         assayName = "logcounts",
                                         donorCol = "donor_id",
                                         ageCol = "age",
                                         cellTypeCol = "celltype",
                                         verbose = TRUE,
                                         seuratAssay = "RNA",
                                         seuratLayer = "data") {
  validSampleTypes <- c("SC", "Pseudobulk", "Bootstrap")
  
  if (!sampleType %in% validSampleTypes) {
    stop(
      "'sampleType' must be one of: ",
      paste(validSampleTypes, collapse = ", ")
    )
  }
  
  brainCtResource <- loadOmniAgeRdata(
    "omniager_brain_celltype_specific_clocks_coef",
    verbose = verbose
  )
  
  brainCtClocksCoef <- brainCtResource[["brain_ct_clocks_coef"]]
  brainCtImputationList <- brainCtResource[["Brain_CT_imputation_data_list"]]
  
  if (is.null(brainCtClocksCoef)) {
    stop("The brain clock resource does not contain 'brain_ct_clocks_coef'.")
  }
  
  if (is.null(brainCtImputationList)) {
    stop(
      "The brain clock resource does not contain ",
      "'Brain_CT_imputation_data_list'."
    )
  }
  
  finalResultsList <- list()
  
  for (ct in cellTypes) {
    clockKey <- paste(sampleType, ct, sep = "_")
    
    modelFolds <- brainCtClocksCoef[[clockKey]]
    imputeData <- brainCtImputationList[[clockKey]]
    
    if (is.null(modelFolds)) {
      warning("Skipping ", ct, ": model not found for ", clockKey)
      next
    }
    
    modelGenes <- unique(unlist(lapply(modelFolds, function(m) {
      setdiff(m$feature_name, "intercept")
    }), use.names = FALSE))
    
    dfBase <- .getDfBrainCt(
      x = x,
      cellType = ct,
      sampleType = sampleType,
      metadata = metadata,
      assayName = assayName,
      donorCol = donorCol,
      ageCol = ageCol,
      cellTypeCol = cellTypeCol,
      seuratAssay = seuratAssay,
      seuratLayer = seuratLayer,
      featureNames = modelGenes
    )
    
    if (nrow(dfBase) == 0L) {
      next
    }
    
    foldPreds <- vapply(names(modelFolds), function(f) {
      res <- .predictBrainCtAge(
        inputData = dfBase,
        imputeData = imputeData,
        modelObj = modelFolds[[f]],
        sampleType = sampleType
      )
      
      res$prediction
    }, numeric(nrow(dfBase)))
    
    finalResultsList[[ct]] <- data.frame(
      prediction = rowMeans(foldPreds, na.rm = TRUE),
      age = dfBase$age,
      donorId = dfBase$donorId,
      sampleType = sampleType,
      celltype = ct,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(finalResultsList) == 0L) {
    return(data.frame())
  }
  
  res <- do.call(rbind, finalResultsList)
  rownames(res) <- NULL
  
  res
}





