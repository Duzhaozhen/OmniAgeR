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
