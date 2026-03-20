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
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
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
loadOmniAgeRdata <- function(title, verbose = TRUE) {
    # 1. Make sure that the development version of OmniAgeRData has been
    # installed locally.
    if (!requireNamespace("OmniAgeRData", quietly = TRUE)) {
        stop(
            "[OmniAgeR] OmniAgeRData package is required but not installed. ",
            "Please install it first."
        )
    }

    # 2. Ultra-fast development mode: Read directly from the local folder
    # Define the local test path and check if the directory exists
    dockerPath <- "/data/OmniAgeData"
    localDevPath <- "~/AgingBiomarker_work/GitHub/OmniAgeR_backup_data/data28Feb26_rds"
    if (dir.exists(dockerPath)) {
      devPath <- dockerPath
    } else {
      devPath <- localDevPath
    }
    
    # devtools::install("~/AgingBiomarker_work/GitHub/OmniAgeRData")
    if (dir.exists(devPath)) {
        # Attempt to load locally; return if successful, otherwise proceed to the next step
        res <- tryCatch(
            {
                OmniAgeRData::getOmniAgeData(title, localTest = TRUE, devPath, verbose)
            },
            error = function(e) NULL
        )

        if (!is.null(res)) {
            return(res)
        }
    }

    # 3. Attempt to use the standard ExperimentHub interface of the data package
    res <- tryCatch(
        {
            OmniAgeRData::getOmniAgeData(title)
        },
        error = function(e) NULL
    )

    if (!is.null(res)) {
        return(res)
    }

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
    stop(
        "[OmniAgeR] Model ", title,
        " could not be loaded via local, Hub, or Cache."
    )
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
