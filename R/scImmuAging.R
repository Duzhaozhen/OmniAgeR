#' @title
#' Run the scImmuAging prediction pipeline for multiple cell types
#'
#' @description
#' This function takes as input a log-normalized matrix of scRNA-seq data and
#' will return cell type-specific predicted age.
#'
#' @param x A \code{SingleCellExperiment}, \code{SummarizedExperiment}, Seurat
#' object, matrix, or data.frame containing log-normalized expression values.
#' Rows should be genes and columns should be cells. For Bioconductor workflows,
#' \code{SingleCellExperiment} input is recommended.
#' @param cellTypes Character vector specifying the cell types for which to
#' predict age. Valid model names include \code{"CD4T"}, \code{"CD8T"},
#' \code{"MONO"}, \code{"NK"}, and \code{"B"} when present in the model object.
#' @param metadata Optional metadata data.frame. Required when \code{x} is a
#' matrix or data.frame.
#' @param assayName Assay name to use for \code{SingleCellExperiment} or
#' \code{SummarizedExperiment} input. Default is \code{"logcounts"}.
#' @param donorCol Column in metadata containing donor IDs.
#' @param ageCol Column in metadata containing donor ages.
#' @param cellTypeCol Column in metadata containing cell type labels.
#' @param pseudocellSize Number of cells sampled to generate each pseudocell.
#' @param pseudocellN Number of pseudocells generated per donor.
#' @param replace Sampling strategy passed to \code{.pseudocellScImmuAging()}.
#' @param minCoverage Numeric between 0 and 1. Minimum required feature coverage.
#' @param verbose Logical. Whether to print status messages.
#' @param seuratAssay Assay name to use for Seurat input.
#' @param seuratLayer Layer or slot name to use for Seurat input.
#'
#' @details
#' This function is designed to work directly with Bioconductor data containers,
#' especially \code{SingleCellExperiment}. Expression values are extracted from
#' an assay, and cell-level metadata are extracted from \code{colData()}.
#'
#' Seurat input is also supported when the Seurat package is installed, but
#' Seurat is not required for Bioconductor workflows.
#' @return A list where each element is named by a cell type. Each element is
#' itself a list with two data.frames:
#' \describe{
#'   \item{bootstrapCell}{Pseudocell-level age predictions.}
#'   \item{donor}{Donor-level aggregated age predictions.}
#' }
#'
#'
#' @references
#' Li W, Zhang Z, Kumar S, et al.
#' Single-cell immune aging clocks reveal inter-individual heterogeneity during
#' infection and vaccination.
#' \emph{Nat Aging} 2025
#'
#' @export
#'
#'
#' @examples
#' # ====================================================================
#' # Example 1: SingleCellExperiment input, recommended for Bioconductor
#' # workflows
#' # ====================================================================
#' 
#' data("ScPbmcExample")
#' set.seed(42)
#' sce_res <- scImmuAging(
#'     x = ScPbmcExample,
#'     cellTypes = c("CD4T"),
#'     assayName = "logcounts",
#'     donorCol = "donor_id",
#'     ageCol = "age",
#'     cellTypeCol = "celltype",
#'     verbose = FALSE
#' )
#'
#' names(sce_res)
#' sce_res$CD4T$donor
#'
#' # ====================================================================
#' # Example 2: Matrix input with separate metadata
#' # ====================================================================
#'
#' \dontrun{
#'  expr_mat <- as.matrix(SummarizedExperiment::assay(ScPbmcExample, "logcounts"))
#'  cell_meta <- as.data.frame(SummarizedExperiment::colData(ScPbmcExample))
#'  set.seed(42)
#'  matrix_res <- scImmuAging(
#'       x = expr_mat,
#'       metadata = cell_meta,
#'       cellTypes = "CD4T",
#'       donorCol = "donor_id",
#'       ageCol = "age",
#'       cellTypeCol = "celltype",
#'       verbose = FALSE
#' )
#'
#'  matrix_res$CD4T$donor
#' }
#'
#' # ====================================================================
#' # Example 3: Optional Seurat input
#' # ====================================================================
#'
#' \dontrun{
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     cell_meta <- as.data.frame(SummarizedExperiment::colData(ScPbmcExample))
#'
#'     seurat_obj <- Seurat::CreateSeuratObject(
#'         counts = SummarizedExperiment::assay(ScPbmcExample, "counts"),
#'         meta.data = cell_meta
#'     )
#'
#'     seurat_obj <- Seurat::SetAssayData(
#'         object = seurat_obj,
#'         assay = "RNA",
#'         layer = "data",
#'         new.data = SummarizedExperiment::assay(ScPbmcExample, "logcounts")
#'     )
#'     set.seed(42)
#'     seurat_res <- scImmuAging(
#'         x = seurat_obj,
#'         cellTypes = "CD4T",
#'         seuratAssay = "RNA",
#'         seuratLayer = "data",
#'         donorCol = "donor_id",
#'         ageCol = "age",
#'         cellTypeCol = "celltype",
#'         verbose = FALSE
#'     )
#'
#'     seurat_res$CD4T$donor
#' }
#' }
#' 
# -------------------------------------------------------------------------
# CODE ATTRIBUTION NOTE:
# The core logic of this function was adapted from the original script 
# provided by https://github.com/CiiM-Bioinformatics-group/scImmuAging
# under the Apache License 2.0.
# Modifications: Added generic object support (SingleCellExperiment/Seurat), 
# refactored the extraction pipeline, and standardized variable names.
# -------------------------------------------------------------------------
scImmuAging <- function(x,
                        cellTypes,
                        metadata = NULL,
                        assayName = "logcounts",
                        donorCol = "donor_id",
                        ageCol = "age",
                        cellTypeCol = "celltype",
                        pseudocellSize = 15,
                        pseudocellN = 100,
                        replace = "dynamic",
                        minCoverage = 0,
                        verbose = TRUE,
                        seuratAssay = "RNA",
                        seuratLayer = "data") {
  if (verbose) {
    message("[scImmuAging] Extracting expression matrix and metadata...")
  }
  
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
  
  if (verbose) {
    message("[scImmuAging] Loading scImmuAging model resources...")
  }
  
  scimmuagingModel <- OmniAgeRData::getOmniAgeRData(
    "omniager_scimmuaging_model",
    verbose = verbose
  )
  
  if (!all(c("model_set", "feature_set") %in% names(scimmuagingModel))) {
    stop(
      "The loaded scImmuAging model resource must contain ",
      "'model_set' and 'feature_set'."
    )
  }
  
  allResults <- list()
  
  for (ct in cellTypes) {
    if (verbose) {
      message("\n--- Processing cell type: ", ct, " ---")
    }
    
    if (!ct %in% names(scimmuagingModel$model_set)) {
      warning("Cell type '", ct, "' model not found. Skipping.")
      next
    }
    
    if (!ct %in% names(scimmuagingModel$feature_set)) {
      warning("Cell type '", ct, "' feature set not found. Skipping.")
      next
    }
    
    currentModel <- scimmuagingModel$model_set[[ct]]
    currentMarkerGenes <- unique(scimmuagingModel$feature_set[[ct]])
    
    keepCells <- metadata[[cellTypeCol]] == ct
    keepCells[is.na(keepCells)] <- FALSE
    
    if (!any(keepCells)) {
      warning("No cells found for cell type '", ct, "'. Skipping.")
      next
    }
    
    preprocessedData <- .scImmuAgingMakePseudocells(
      expr = expr,
      metadata = metadata,
      cellType = ct,
      markerGenes = currentMarkerGenes,
      donorCol = donorCol,
      ageCol = ageCol,
      cellTypeCol = cellTypeCol,
      pseudocellSize = pseudocellSize,
      pseudocellN = pseudocellN,
      replace = replace,
      verbose = verbose
    )
    
    predResults <- .scImmuAgingCalculator(
      preprocessedData = preprocessedData,
      model = currentModel,
      markerGenes = currentMarkerGenes,
      minCoverage = minCoverage,
      verbose = verbose
    )
    
    predResults$celltype <- ct
    
    agePerDonor <- .ageDonor(predResults)
    agePerDonor$celltype <- ct
    
    allResults[[ct]] <- list(
      bootstrapCell = predResults,
      donor = agePerDonor
    )
  }
  
  return(allResults)
}