
#' Preprocess and Impute DNA Methylation Data for PC Clocks
#'
#' @description
#' This internal function ensures the input methylation matrix contains all 
#' required CpGs for PC-based clock projection. It performs a hierarchical 
#' quality control and imputation:
#' \enumerate{
#'   \item Checks if the proportion of found CpGs meets the \code{minCoverage} 
#'   threshold.
#'   \item Performs column-mean imputation for any \code{NA} values within 
#'   present CpGs.
#'   \item Appends entirely missing CpGs using reference values (medians/means) 
#'         provided in the \code{requiredCpGs} vector.
#' }
#'
#' @param betaM A numeric matrix of DNA methylation values.
#' @param requiredCpGs A named numeric vector of required CpG sites and their 
#'   reference values.
#' @param minCoverage A numeric value (0-1) representing the minimum required 
#'   proportion of CpGs.
#' @param verbose A logical flag to print progress messages.
#'
#' @return A numeric matrix with dimensions (N_samples x N_required_CpGs), 
#'   where all columns match the \code{requiredCpGs} order. Returns \code{NULL} 
#'   if the coverage check fails.
#'
#' @keywords internal
#' @noRd


.preprocessPcData <- function(betaM, requiredCpGs, minCoverage, verbose) {
  
  targetNames <- names(requiredCpGs)
  presentCpGs <- intersect(colnames(betaM), targetNames)
  
  # 1. Coverage Check
  coverageRatio <- length(presentCpGs) / length(targetNames)
  if (verbose) {
    message(sprintf("Probe Check: Found %d / %d required CpGs (%.1f%%).", 
                    length(presentCpGs), length(targetNames), coverageRatio * 100))
  }
  
  if (coverageRatio < minCoverage) {
    if (verbose) warning("Coverage below threshold. Calculation aborted.")
    return(NULL)
  }
  
  # 2. Extract and Mean Impute existing CpGs
  betaSub <- betaM[, presentCpGs, drop = FALSE]
  if (any(is.na(betaSub))) {
    if (verbose) message("[PCClocks] Imputing missing values within present CpGs using column means...")
    # Vectorized column mean imputation
    colMeansV <- colMeans(betaSub, na.rm = TRUE)
    naIdx <- which(is.na(betaSub), arr.ind = TRUE)
    betaSub[naIdx] <- colMeansV[naIdx[, 2]]
  }
  
  # 3. Append missing CpGs using Reference Medians/Means
  missingCpGs <- setdiff(targetNames, presentCpGs)
  if (length(missingCpGs) > 0) {
    if (verbose) message(sprintf("[PCClocks] Appending %d missing CpGs with reference values...", length(missingCpGs)))
    imputeMat <- matrix(rep(requiredCpGs[missingCpGs], each = nrow(betaM)), 
                        nrow = nrow(betaM), ncol = length(missingCpGs))
    colnames(imputeMat) <- missingCpGs
    betaSub <- cbind(betaSub, imputeMat)
  }
  
  return(betaSub[, targetNames, drop = FALSE])
}



#' Robustly Load Serialized Package Data (.qs2)
#'
#' An internal helper function to safely load .qs2 data files from a local
#' path, validate their existence, and provide user-friendly error messages.
#'
#' @param object_name The base name of the data object (e.g., "PCClocks_data")
#'   to be loaded, excluding the ".qs2" extension.
#' @param path An optional character string specifying either the full path
#'   to the .qs2 file or the path to the directory containing it. If NULL,
#'   defaults to the path from \code{getOmniAgeRPath()}.
#'
#' @return The deserialized R object (e.g., a list or data.frame)
#'   from the .qs2 file.
#' @export
#'
#'

#load_PCClocks_data
load_OmniAgeR_data <- function (object_name,
                                      path = NULL){
  
  checkmate::assert_string(object_name)
  
  full_object_name <- paste0(object_name, ".qs2")
  checkmate::assert_character(path, len = 1, null.ok = TRUE)
  
  file_path <- if (!is.null(path)) {
    stopifnot(`Provided file/dir doesn't exists` = dir.exists(path) ||
                file.exists(path))
    info <- file.info(path)
    if (info$isdir) {
      file.path(path, full_object_name)
    }
    else {
      path
    }
  }
  else {
    file.path(getOmniAgeRPath(), full_object_name)
  }
  
  
  tryCatch({
    obj <- qs2::qs_read(file_path, validate_checksum = TRUE)
  }, error = function(e) {
    
    stop(sprintf(paste("Failed to read '%s'.",
                       "Did you download '%s' and passed the path or loaded object to this function?",
                       "See function `download_biomarker_data()` to download the required data.",
                       "Function failed: %s", sep = "\n"),
                 full_object_name, full_object_name, e))
  })
  
  return(obj)
}




#' Get Default Data Storage Directory
#'
#' Retrieves the OS-specific, cross-platform data directory for the
#' 'OmniAgeR' package using \code{tools::R_user_dir}.
#' Creates the directory if it does not already exist.
#'
#' @return A character string representing the file path to the
#'   package's data directory.
#'
#' @export
#'
getOmniAgeRPath <- function() {
  path <- tools::R_user_dir("OmniAgeR", which = "data")

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  return(path)
}


#' Download OmniAgeR Data
#'
#' Downloads large data files (.qs2) required for clock calculations.
#' Data is downloaded from the Zenodo repository.
#'
#' @param clocks A character vector of clock datasets to download.
#'   Valid options (based on the function's default) are "SystemsAge" and "PCClocks".
#' @param path A character string specifying the directory to save the data.
#'   If NULL (default), uses the path from \code{getOmniAgeRPath()}.
#' @param force A logical value. If TRUE, data will be re-downloaded
#'   and will overwrite any existing files. Default is FALSE.
#' @param ... Additional arguments passed to \code{zen4R::ZenodoRecord$downloadFiles()}.
#'
#' @return Invisibly, a logical vector where TRUE indicates
#'   a successful download for each file.
#' @export
download_OmniAgeR_data <- function(clocks = c("SystemsAge", "PCClocks"),
                                   path = NULL,
                                   force = FALSE,
                                   ...) {
  BIOMARKER_MANIFEST <- data.frame(
    clock_name = c("SystemsAge", "PCClocks"),
    file_name = c("SystemsAge_data.qs2", "PCClocks_data.qs2")
  )

  ZENODO_DOI <- "10.5281/zenodo.17162604" #


  clocks <- match.arg(clocks, several.ok = TRUE)

  if (is.null(path)) {
    path <- getOmniAgeRPath()
    if (!dir.exists(path)) {
      message("Creating a new folder at ", path)
      if (!dir.create(path, showWarnings = TRUE, recursive = TRUE)) {
        stop("Failed to create directory at ", path)
      }
    }
  }
  checkmate::assert_directory_exists(path, access = "rw", .var.name = "path")


  if ("all" %in% clocks) {
    clocks <- BIOMARKER_MANIFEST$clock_name
  }
  clocks <- unique(clocks)


  manifest_subset <- BIOMARKER_MANIFEST[BIOMARKER_MANIFEST$clock_name %in% clocks, ]

  if (nrow(manifest_subset) == 0 && length(clocks) > 0) {
    warning("No valid clock names provided. Check spelling.")
    return(invisible(FALSE))
  }

  download_name <- manifest_subset$file_name
  download_to <- file.path(path, download_name)

  exists <- if (force) {
    rep(FALSE, times = length(download_to))
  } else {
    file.exists(download_to)
  }

  if (any(exists)) {
    message(paste("Skipping", sum(exists), "file(s) that already exist. Use force = TRUE to re-download."))

    clocks <- clocks[!exists]
    download_name <- download_name[!exists]
    download_to <- download_to[!exists]
  }

  if (length(clocks) == 0) {
    message("No new files to download.")
    return(invisible(TRUE))
  }


  if (!requireNamespace("zen4R", quietly = TRUE)) {
    stop("Please install the 'zen4R' package to download these files (`install.packages('zen4R')`).")
  }

  message("Connecting to Zenodo...")
  zenodo <- zen4R::ZenodoManager$new(logger = "INFO")

  rec <- zenodo$getRecordByDOI(ZENODO_DOI)

  if (is.null(rec)) {
    stop(paste("Could not find Zenodo record for DOI:", ZENODO_DOI))
  }


  message(paste("Attempting to download", length(clocks), "file(s) to:", path))
  download_success <- logical(length(clocks))
  names(download_success) <- download_name

  for (i in seq_along(clocks)) {
    tryCatch(
      {
        current_file_name <- download_name[i]
        message("Downloading ", current_file_name, "...")


        rec$downloadFiles(
          path = path,
          files = list(current_file_name),
          ...
        )

        if (!file.exists(download_to[i])) {
          stop(paste("File does not exist after download:", download_to[i]))
        }

        message("Successfully downloaded ", current_file_name)
        download_success[i] <- TRUE
      },
      error = function(e) {
        message("Error: ", e$message)
        warning(paste("Failed to download:", download_name[i]), call. = FALSE)
        download_success[i] <- FALSE
      }
    )
  }

  return(invisible(download_success))
}
