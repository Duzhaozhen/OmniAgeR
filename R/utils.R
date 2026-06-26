#' Load Pre-trained Models and Example Data for OmniAgeR
#'
#' @description
#' This function seamlessly loads specific aging omic clock models, weights, or
#' example datasets required by the \pkg{OmniAgeR} package.
#'
#' @param title A character string specifying the exact name of the model
#' or resource to load (e.g., \code{"omniager_horvath2013_coef"}).
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' To comply with Bioconductor guidelines and minimize the software
#' package size, heavy data files are managed externally via ExperimentHub.
#'
#' @return An R object (typically a \code{data.frame}, \code{matrix},
#' or \code{list}) containing the requested model parameters or reference data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the Horvath2013 model weights
#' horvath2013Model <- loadOmniAgeRdata("omniager_horvath2013_coef")
#' head(horvath2013Model)
#' }
loadOmniAgeRdata <- function(title, verbose = TRUE) {
  # 1. Make sure that the underlying Hub dependency packages have been installed
  if (!requireNamespace("ExperimentHub", quietly = TRUE) ||
      !requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("[OmniAgeR] 'ExperimentHub' and 'AnnotationHub' are required to load data.")
  }
  
  # 2. Instantiate ExperimentHub and query
  eh <- ExperimentHub::ExperimentHub()
  res <- AnnotationHub::query(eh, c("OmniAgeRData", title))
  
  if (length(res) == 0) {
    stop(sprintf("[OmniAgeR] Resource '%s' not found in ExperimentHub.", title))
  }
  
  # 3. Instantiate ExperimentHub and query
  exact_idx <- which(res$title == title)
  if (length(exact_idx) == 0) {
    exact_idx <- 1
  } else {
    exact_idx <- exact_idx[1]
  }
  
  hubTitle <- res$title[exact_idx]
  
  if (verbose) {
    message("[OmniAgeR] Retrieving resource: ", hubTitle)
  }
  
  # 4. Instantiate ExperimentHub and query
  dataObjOrPath <- res[[exact_idx]]
  
  # 5. Special parsing for the qs2 format
  if (is.character(dataObjOrPath) && grepl("\\.qs2?$", hubTitle)) {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop("[OmniAgeR] Package 'qs2' is required to read this resource.")
    }
    dataObjOrPath <- qs2::qs_read(dataObjOrPath)
  }
  
  if (verbose) {
    message("[OmniAgeR] Successfully loaded '", title, "'.")
  }
  
  return(dataObjOrPath)
}

#' Developmental Age Transformation
#'
#' @param x A vector of sample ages
#' @param adultAge The age considered to be the cutoff for adulthood
#'
#' @return transformed age prediction
#' @noRd
.antiTrafo <- function(x, adultAge = 20) {
    ifelse(x < 0, (1 + adultAge) * exp(x) - 1, (1 + adultAge) * x + adultAge)
}




#' Extract Assay Matrix from Various Object Types
#'
#' @description An internal helper function to gracefully extract a numeric
#' matrix from matrix, data.frame.
#'
#' @param x The input object.
#' @param assayName Optional character specifying which assay to extract.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
.extractAssayMatrix <- function(x, assayName = NULL) {
  # 1. Bioconductor Object (SummarizedExperiment)
  if (inherits(x, "SummarizedExperiment")) {
    # Ensure the package is available
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required to process this object.")
    }
    if (is.null(assayName)) {
      return(SummarizedExperiment::assay(x)) # defaults to the first assay
    } else {
      return(SummarizedExperiment::assay(x, assayName))
    }
  } 
  # 2. Base R Matrix or Data Frame
  else if (inherits(x, c("matrix", "data.frame"))) {
    return(as.matrix(x))
  } 
  
  # 3. Unknown Object
  else {
    stop("Input must be a matrix, data.frame, or SummarizedExperiment.")
  }
}


#' Extract and standardize single-cell expression data
#'
#' @description
#' This internal utility function processes various single-cell input formats,
#' including \code{SingleCellExperiment}, \code{Seurat} objects, or raw matrices,
#' into a standardized format. It performs necessary data coercion, column
#' validation, and alignment between the expression matrix and cell metadata.
#'
#' @param x A \code{SingleCellExperiment}, \code{SummarizedExperiment},
#'   \code{Seurat} object, \code{matrix}, or \code{data.frame} containing
#'   the single-cell expression data.
#' @param metadata A \code{data.frame} containing cell-level metadata.
#'   Required only if \code{x} is a \code{matrix} or \code{data.frame}.
#'   Defaults to \code{NULL}.
#' @param assayName A string specifying the assay name to extract if \code{x}
#'   is a \code{SingleCellExperiment} or \code{SummarizedExperiment}.
#'   Defaults to \code{"logcounts"}.
#' @param donorCol A string specifying the column name in metadata representing
#'   the donor identity. Defaults to \code{"donor_id"}.
#' @param ageCol A string specifying the column name in metadata representing
#'   the subject age. Defaults to \code{"age"}.
#' @param cellTypeCol A string specifying the column name in metadata representing
#'   the cell type. Defaults to \code{"celltype"}.
#' @param seuratAssay A string specifying the assay to extract if \code{x}
#'   is a \code{Seurat} object. Defaults to \code{"RNA"}.
#' @param seuratLayer A string specifying the layer or slot to extract if \code{x}
#'   is a \code{Seurat} object. Defaults to \code{"data"}.
#'
#' @return A \code{list} containing two elements:
#' \itemize{
#'   \item \code{expr}: A numeric matrix of gene expression values.
#'   \item \code{metadata}: A \code{data.frame} of the aligned cell metadata.
#' }
#'
#' @details
#' This function performs strict input validation:
#' \enumerate{
#'   \item Ensures the expression values are numeric.
#'   \item Verifies that row names (genes) and column names (cells) are present.
#'   \item Matches the expression matrix columns with the metadata rows.
#'   \item Checks for the presence of mandatory columns (\code{donorCol}, \code{ageCol}).
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames colData
#' @keywords internal

.extractScInput <- function(x,
                            metadata = NULL,
                            assayName = "logcounts",
                            donorCol = "donor_id",
                            ageCol = "age",
                            cellTypeCol = "celltype",
                            seuratAssay = "RNA",
                            seuratLayer = "data") {
  if (inherits(x, "SingleCellExperiment") ||
      inherits(x, "SummarizedExperiment")) {
    
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop(
        "The SummarizedExperiment package is required for ",
        "SingleCellExperiment/SummarizedExperiment input."
      )
    }
    
    assayNames <- SummarizedExperiment::assayNames(x)
    
    if (!assayName %in% assayNames) {
      stop(
        "Assay '", assayName, "' was not found in 'x'. ",
        "Available assays: ", paste(assayNames, collapse = ", ")
      )
    }
    
    expr <- SummarizedExperiment::assay(x, assayName)
    metadata <- as.data.frame(SummarizedExperiment::colData(x))
    
  } else if (inherits(x, "Seurat")) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop(
        "The Seurat package is required for Seurat input. ",
        "For Bioconductor workflows, please use a ",
        "SingleCellExperiment object."
      )
    }
    
    expr <- tryCatch(
      Seurat::GetAssayData(
        x,
        assay = seuratAssay,
        layer = seuratLayer
      ),
      error = function(e) {
        Seurat::GetAssayData(
          x,
          assay = seuratAssay,
          slot = seuratLayer
        )
      }
    )
    
    metadata <- x[[]]
    
  } else if (is.matrix(x) || is.data.frame(x)) {
    
    if (is.null(metadata)) {
      stop(
        "When 'x' is a matrix or data.frame, 'metadata' must also ",
        "be provided."
      )
    }
    
    expr <- as.matrix(x)
    metadata <- as.data.frame(metadata)
    
  } else {
    stop(
      "'x' must be a SingleCellExperiment, SummarizedExperiment, ",
      "Seurat object, matrix, or data.frame."
    )
  }
  
  expr <- as.matrix(expr)
  
  if (!is.numeric(expr)) {
    stop("The selected assay must contain numeric expression values.")
  }
  
  if (is.null(rownames(expr))) {
    stop("The expression matrix must have gene names as row names.")
  }
  
  if (is.null(colnames(expr))) {
    stop("The expression matrix must have cell names as column names.")
  }
  
  if (nrow(metadata) != ncol(expr)) {
    stop(
      "The number of rows in metadata must match the number of ",
      "columns in the expression matrix."
    )
  }
  
  if (!is.null(rownames(metadata)) &&
      all(colnames(expr) %in% rownames(metadata))) {
    metadata <- metadata[colnames(expr), , drop = FALSE]
  } else {
    rownames(metadata) <- colnames(expr)
  }
  
  requiredCols <- c(donorCol, ageCol)
  if (!is.null(cellTypeCol)) {
    requiredCols <- c(requiredCols, cellTypeCol)
  }
  
  missingCols <- setdiff(requiredCols, colnames(metadata))
  
  if (length(missingCols) > 0L) {
    stop(
      "The metadata must contain the following columns: ",
      paste(missingCols, collapse = ", ")
    )
  }
  
  list(
    expr = expr,
    metadata = metadata
  )
}







#' Extract Phenotype Data from Various Object Types
#'
#' @description An internal helper function to gracefully extract metadata
#' (phenotype data) from SummarizedExperiment, Seurat, or list objects.
#'
#' @param x The input object.
#' @return A data.frame containing sample metadata, or NULL if not applicable.
#' @keywords internal
#' @noRd
.extractPhenoData <- function(x) {
  # 1. Bioconductor Object (SummarizedExperiment / SingleCellExperiment)
  if (inherits(x, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required to process this object.")
    }
    return(as.data.frame(SummarizedExperiment::colData(x)))
  } 
  
  # 2. Seurat Object
  else if (inherits(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required to process Seurat objects. Please install it.")
    }
    return(x[[]]) # Extacts the meta.data data.frame
  } 
  
  # 3. Base R list (like your dnamExample)
  else if (inherits(x, "list") && length(x) >= 2 && is.data.frame(x[[2]])) {
    return(x[[2]])
  } 
  
  # 4. Matrix or Data Frame (No attached phenotype data)
  else if (inherits(x, c("matrix", "data.frame"))) {
    return(NULL) # Matrices don't store separate phenotype data
  } 
  
  else {
    return(NULL)
  }
}








