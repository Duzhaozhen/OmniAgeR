#' Load Pre-trained Models and Example Data for OmniAgeR
#'
#' @description
#' This function seamlessly loads specific aging omic clock models, weights, or 
#' example datasets required by the \pkg{OmniAgeR} package. It automatically 
#' retrieves the data from the companion data package \pkg{OmniAgeRData} via 
#' ExperimentHub.
#'
#' @param title A character string specifying the exact name of the model 
#' or resource to load (e.g., \code{"omniager_horvath2013_coef"}).
#'
#' @details 
#' To comply with Bioconductor guidelines and minimize the software 
#' package size, heavy data files are managed externally.
#'
#' @return An R object (typically a \code{data.frame}, \code{matrix}, 
#' or \code{list}) containing the requested model parameters or reference data.
#'
#' @export
#'
#' @examples
#' # Ensure OmniAgeRData is installed before running
#' # Load the Horvath2013 model weights
#' horvath2013Model <- loadOmniAgeRdata("omniager_horvath2013_coef")
#' 
#' # View the first few rows of the loaded model parameters
#' head(horvath2013Model)

loadOmniAgeRdata <- function(title) {
  
  # 1. Make sure that the development version of OmniAgeRData has been 
  # installed locally.
  if (!requireNamespace("OmniAgeRData", quietly = TRUE)) {
    stop("[OmniAgeR] OmniAgeRData package is required but not installed. ",
         "Please install it first.")
  }
  
  # 2. Ultra-fast development mode: Read directly from the local folder
  # Define the local test path and check if the directory exists
  devPath <- "~/AgingBiomarker_work/GitHub/OmniAgeR_backup_data/data28Feb26_rds"
  
  if (dir.exists(devPath)) {
    # Attempt to load locally; return if successful, otherwise proceed to the next step
    res <- tryCatch({
      OmniAgeRData::getOmniAgeData(title, localTest = TRUE, devPath)
    }, error = function(e) NULL)
    
    if (!is.null(res)) return(res)
  }
  
  # 3. Attempt to use the standard ExperimentHub interface of the data package
  res <- tryCatch({
    OmniAgeRData::getOmniAgeData(title)
  }, error = function(e) NULL)
  
  if (!is.null(res)) return(res)
  
  # 4. Simulate ExperimentHub mode (core fallback mechanism)
  # Read metadata.csv from the locally installed OmniAgeRData package
  metaFile <- system.file("extdata", "metadata.csv", package = "OmniAgeRData")
  
  if (file.exists(metaFile)) {
    metaData <- read.csv(metaFile, stringsAsFactors = FALSE)
    
    # Find the corresponding resource row based on the title
    targetRow <- metaData[metaData$Title == title, ]
    
    if (nrow(targetRow) > 0) {
      # Construct the complete download URL
      targetUrl <- paste0(targetRow$Location_Prefix, targetRow$RDataPath) 
      
      # Ensure BiocFileCache is installed (must be in Imports of DESCRIPTION)
      if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
        stop("[OmniAgeR] BiocFileCache is required for downloading data.")
      }
      
      # Use BiocFileCache to simulate cached download (downloads once, 
      # loads instantly thereafter)
      bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
      cachedPath <- BiocFileCache::bfcrpath(bfc, targetUrl)
      
      return(readRDS(cachedPath))
    }
  }
  
  # 5. If all the above methods fail, throw a clear error
  stop("[OmniAgeR] Model ", title, 
       " could not be loaded via local, Hub, or Cache.")
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



#' Download OmniAgeR Example Data
#'
#' Downloads large example datasets (.rda) required for vignettes and examples.
#' Data is stored in the package's user data directory.
#'
#' @param datasets A character vector of datasets to download (e.g., "seu_gabitto_2024_filtered").
#' @param path Path to save data. Defaults to \code{getOmniAgeRPath()}.
#' @param force Logical. If TRUE, overwrites existing files.
#' @param ... Additional arguments passed to \code{zen4R::ZenodoRecord$downloadFiles}.
#' @return Invisible TRUE upon successful completion or if files already exist
#' @export
#' @examples
#' downloadOmniAgeRExample("Hannum_example")
downloadOmniAgeRExample <- function(datasets = c(
                                        "Hannum_example.rda", "GA_example.rda", "LungInv.rda", "Tursiops_example.rda", "TZH_example_CTF.rda", "brain_frohlich_control_example_15donors.rda",
                                        "seu_gabitto_2024_filtered", "Yazar_CD4T_CD8T_example.rda", "CTS_ExampleData_Liver.rda", "CTS_MurphyGSE88890.rda",
                                        "CTS_PaiGSE112179.rda"
                                      ),
                                      path = NULL,
                                      force = FALSE, ...) {
  # List defining the sample data
  exampleManifest <- data.frame(
    name = c(
      "Hannum_example", "GA_example", "LungInv", "Tursiops_example",
      "TZH_example_CTF", "brain_frohlich_control_example_15donors",
      "seu_gabitto_2024_filtered", "Yazar_CD4T_CD8T_example", "CTS_ExampleData_Liver", "CTS_MurphyGSE88890",
      "CTS_PaiGSE112179"
    ),
    file = c(
      "Hannum_example.rda", "GA_example.rda", "LungInv.rda", "Tursiops_example.rda",
      "TZH_example_CTF.rda", "brain_frohlich_control_example_15donors.rda",
      "seu_gabitto_2024_filtered.rda", "Yazar_CD4T_CD8T_example.rda", "CTS_ExampleData_Liver.rda", "CTS_MurphyGSE88890.rda",
      "CTS_PaiGSE112179.rda"
    )
  )
  
  zenodoDoi <- "10.5281/zenodo.18287372"
  
  if (is.null(path)) path <- getOmniAgeRPath()
  
  filesToDownload <- exampleManifest$file[exampleManifest$name %in% datasets]
  
  if (length(filesToDownload) == 0) {
    filesToDownload <- exampleManifest$file[exampleManifest$file %in% datasets]
  }
  
  if (length(filesToDownload) == 0) {
    stop(paste("The dataset", datasets, "is not in the manifest. Please check the spelling."))
  }
  
  if (!force) {
    filesToDownload <- filesToDownload[!file.exists(file.path(path, filesToDownload))]
  }
  
  if (length(filesToDownload) == 0) {
    message("All requested example datasets already exist.")
    return(invisible(TRUE))
  }
  
  message("Connecting to Zenodo to download example datasets...")
  
  if (!requireNamespace("zen4R", quietly = TRUE)) {
    stop("Package 'zen4R' is required to download data. Please install it.")
  }
  
  zenodo <- zen4R::ZenodoManager$new(logger = "INFO")
  rec <- zenodo$getRecordByDOI(zenodoDoi)
  
  for (f in filesToDownload) {
    message("Downloading example: ", f, "...")
    rec$downloadFiles(path = path, files = list(f), ...)
  }
  
  return(invisible(TRUE))
}
#' Load OmniAgeR Example Data
#'
#' Safely load .rda example datasets from the local storage path.
#'
#' @param dataset_name Base name of the dataset (e.g., "seu_gabitto_2024_filtered").
#' @param envir The environment where the data should be loaded. Defaults to parent.frame().
#' @return Invisible NULL. The dataset is loaded into the specified environment as a side effect.
#' @export
#' @examples
#' loadOmniAgeRExample("Hannum_example")
loadOmniAgeRExample <- function(datasetName, envir = parent.frame()) {
  fullName <- paste0(datasetName, ".rda")
  # 假设 getOmniAgeRPath 也改为了 getOmniAgeRPath
  filePath <- file.path(getOmniAgeRPath(), fullName)

  if (!file.exists(filePath)) {
    stop(sprintf(
      "Example data '%s' not found. Please run `downloadOmniAgeRExample('%s')` first.",
      fullName, datasetName
    ))
  }

  # load .rda file
  load(filePath, envir = envir)
  message(sprintf("Dataset '%s' loaded into environment.", datasetName))
  return(invisible(NULL))
}
