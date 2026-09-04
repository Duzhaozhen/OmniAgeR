#' Retrieve Feature Weights and Intercepts for Epigenetic Clocks
#'
#' @description
#' A high-level API to access the underlying parameters 
#' for all aging omic markers implemented in the \code{OmniAgeR} package. This 
#' function supports batch retrieval by individual clock names or predefined 
#' biological categories.
#'
#' @param clockNames A character vector of marker identifiers. 
#' Options include:
#' \itemize{
#'   \item \strong{Individual Names:} e.g., \code{"Horvath2013"}, \code{"PhenoAge"}.
#'   \item \strong{Categories:} e.g., \code{"biological"}, \code{"chronological"}, \code{"mitotic"}.
#'   \item \strong{Keyword:} \code{"all"} (default) to retrieve all available clocks.
#' }
#' @param verbose A logical flag. If `FALSE` (default), prints status messages.
#' 
#' 
#' @return A named \code{list}. Each element corresponds to a clock and contains a 
#' \code{data.frame} with the following columns:
#' \describe{
#'   \item{feature/probe}{Character. The Illumina probe identifier or feature name.}
#'   \item{coef}{Numeric. The regression coefficient or relative weight.}
#' }
#'
#'
#' @examples
#' # Get weights for specific clocks
#' weights <- getMarkerWeights(c("Horvath2013", "Hannum"))
#' 
#' # Access specific weights
#' head(weights$Horvath2013)
#' @export
getMarkerWeights <- function(clockNames = "all", verbose = FALSE) {
  
  # 1. Resolve clock categories and identifiers
  clockMap <- .getMarkerCategories()
  clocksToGet <- .resolveClockLogic(clockNames, clockMap)
  
  # 2. Check for PC-based clocks to manage memory-intensive data loading
  # These clocks require projection from PC-space back to CpG-space
  pcClockNames <- c("PCHorvath2013", "PCHorvath2018", "PCHannum", 
                    "PCPhenoAge", "PCDNAmTL", "PCGrimAge1")
  pcRequested <- intersect(clocksToGet, pcClockNames)
  
  pcClockData <- NULL
  if (length(pcRequested) > 0) {
    pcClockData <- OmniAgeRData::getOmniAgeRData("PCClocks_data", verbose = verbose)
  }
  
  # 3. Dispatch and aggregate weights across requested clocks
  weightsList <- lapply(clocksToGet, function(clockLabel) {
    .dispatchClockWeights(clockLabel, pcClockData, verbose)
  })
  names(weightsList) <- clocksToGet
  
  # 4. Cleanup: Remove NULLs (unimplemented or missing clocks)
  weightsList <- weightsList[!vapply(weightsList, is.null, logical(1))]
  
  return(weightsList)
}
#' Internal Central Dispatcher for Clock Coefficients and Models
#' 
#' @description
#' A centralized internal routing function that retrieves pre-trained model 
#' coefficients, weights, and parameters for various epigenetic clocks. It 
#' acts as an abstraction layer to handle diverse clock architectures—including 
#' standard linear models, Principal Component (PC) clocks, cell-type-specific 
#' models, and ensemble methods—standardizing them for downstream calculations.
#' 
#' @param clockLabel A character string specifying the unique identifier of the clock.
#' @param pcClockData A list containing the pre-calculated PCA rotation and 
#' model data, specifically required for PC-based clocks.
#' @param verbose A logical flag. If `FALSE` (default), prints status messages
#' @return Depending on the \code{clockLabel}, returns either:
#' \itemize{
#'   \item A \code{data.frame} with features, weights, and intercepts.
#'   \item A \code{list} for multi-component or cell-type-specific clocks (e.g., \code{scImmuAging}, \code{PASTA}).
#'   \item \code{NULL} if the clock label is unrecognized.
#' }
#' 
#' @keywords internal
.dispatchClockWeights <- function(clockLabel, pcClockData, verbose) {

  ctsClockNames <- c("Neu-In", "Glia-In", "Brain", "Neu-Sin", 
                     "Glia-Sin", "Hep", "Liver")
  if(clockLabel %in% ctsClockNames){
    ctsCoefs <- OmniAgeRData::getOmniAgeRData("omniager_cts_clocks_coef", verbose = verbose)
  }
  
  
  res <- switch(clockLabel,
                "epiTOC1" = OmniAgeRData::getOmniAgeRData("omniager_epitoc1_model", verbose),
                "epiTOC2" = OmniAgeRData::getOmniAgeRData("omniager_epitoc2_model", verbose),
                "epiTOC3" = OmniAgeRData::getOmniAgeRData("omniager_epitoc3_model", verbose),
                "stemTOCvitro" = OmniAgeRData::getOmniAgeRData("omniager_stemtocvitro_cpg", verbose),
                "stemTOC" = OmniAgeRData::getOmniAgeRData("omniager_stemtoc_cpg", verbose),
                "HypoClock" = OmniAgeRData::getOmniAgeRData("omniager_hypoclock_cpg", verbose),
                "RepliTali" = OmniAgeRData::getOmniAgeRData("omniager_replitali_coef", verbose),
                "EpiCMIT_Hyper" = {
                  epiCMITdf <- OmniAgeRData::getOmniAgeRData("omniager_epicmit_model", verbose)
                  hyperProbes <- epiCMITdf[grep("hyper", epiCMITdf$epiCMIT.class), 1]
                  
                  n_probes <- length(hyperProbes)
                  data.frame(
                    CpG = hyperProbes,
                    Weight = 1 / n_probes,  
                    Intercept = 0,
                    stringsAsFactors = FALSE
                  )
                },
                
                "EpiCMIT_Hypo" = {
                  epiCMITdf <- OmniAgeRData::getOmniAgeRData("omniager_epicmit_model", verbose)
                  hypoProbes <- epiCMITdf[grep("hypo", epiCMITdf$epiCMIT.class), 1]
                  
                  n_probes <- length(hypoProbes)
                  data.frame(
                    CpG = hypoProbes,
                    Weight = -1 / n_probes, 
                    Intercept = 1,       
                    stringsAsFactors = FALSE
                  )
                },

                "DNAmTL" = OmniAgeRData::getOmniAgeRData("omniager_dnamtl_coef", verbose),
                "Horvath2013" = OmniAgeRData::getOmniAgeRData("omniager_horvath2013_coef", verbose),
                "Hannum" = OmniAgeRData::getOmniAgeRData("omniager_hannum", verbose),
                "Lin" = OmniAgeRData::getOmniAgeRData("omniager_lin_coef", verbose),
                "VidalBralo" = OmniAgeRData::getOmniAgeRData("omniager_vidalbralo_coef", verbose),
                "ZhangClock" = OmniAgeRData::getOmniAgeRData("omniager_zhang_clock_coef", verbose),
                "Horvath2018" = OmniAgeRData::getOmniAgeRData("omniager_horvath2013_coef", verbose),
                "Bernabeu_cAge" = OmniAgeRData::getOmniAgeRData("omniager_bernabeu_cage_coef", verbose),
                "PedBE"  = OmniAgeRData::getOmniAgeRData("omniager_pedbe_coef", verbose),
                "CorticalClock" = OmniAgeRData::getOmniAgeRData("omniager_cortical_clock_coef", verbose),
                
                "CentenarianClock" = OmniAgeRData::getOmniAgeRData("omniager_centenarian_coef", verbose),
                "Retro_age" = OmniAgeRData::getOmniAgeRData("omniager_retroage_coef", verbose),
                "ABEC" = OmniAgeRData::getOmniAgeRData("omniager_abec_coef", verbose),
                "eABEC" = OmniAgeRData::getOmniAgeRData("omniager_eabec_coef", verbose),
                "cABEC" = OmniAgeRData::getOmniAgeRData("omniager_cabec_coef", verbose),
                "PipekElasticNet" = OmniAgeRData::getOmniAgeRData("omniager_pipek_elasticnet_coef", verbose),
                "PipekFilteredh"  = OmniAgeRData::getOmniAgeRData("omniager_pipek_filteredh_coef", verbose),
                "PipekRetrainedh" = OmniAgeRData::getOmniAgeRData("omniager_pipek_retrainedh_coef", verbose),
                "WuClock" = OmniAgeRData::getOmniAgeRData("omniager_wu_clock_coef", verbose),
                "Weidner" = OmniAgeRData::getOmniAgeRData("omniager_weidner_coef", verbose),
                "IntrinClock" = OmniAgeRData::getOmniAgeRData("omniager_intrin_clock_coef", verbose),
                "Garagnani" = OmniAgeRData::getOmniAgeRData("omniager_garagnani_coef", verbose),
                
                "Zhang10" = OmniAgeRData::getOmniAgeRData("omniager_zhang10_coef", verbose),
                "PhenoAge" = OmniAgeRData::getOmniAgeRData("omniager_phenoage_coef", verbose),
                "DunedinPACE" = OmniAgeRData::getOmniAgeRData("omniager_dunedinpace_model", verbose),
                "GrimAge1" = OmniAgeRData::getOmniAgeRData("omniager_grimage1_model", verbose),
                "GrimAge2" = OmniAgeRData::getOmniAgeRData("omniager_grimage2_model", verbose),
                "IC_Clock" = OmniAgeRData::getOmniAgeRData("omniager_ic_clock_coef", verbose),
                
                "DNAmFitAge" = OmniAgeRData::getOmniAgeRData("omniager_dnamfitage_coef", verbose),
                # --- PC Clocks ---
                "PCHorvath2013" = pcClockData$CalcPCHorvath1,
                "PCHorvath2018" = pcClockData$CalcPCHorvath2,
                "PCHannum" = pcClockData$CalcPCHannum,
                "PCPhenoAge" = pcClockData$CalcPCPhenoAge,
                "PCDNAmTL" = pcClockData$CalcPCDNAmTL,
                "PCGrimAge1" = pcClockData$CalcPCGrimAge,
                "SystemsAge" = OmniAgeRData::getOmniAgeRData("SystemsAge_data", verbose),

                "CausalAge" = .formatLinearWeights(OmniAgeRData::getOmniAgeRData("omniager_causal_clocks_coef", FALSE)[[1]]),
                "DamAge" = .formatLinearWeights(OmniAgeRData::getOmniAgeRData("omniager_causal_clocks_coef", FALSE)[[2]]),
                "AdaptAge" = .formatLinearWeights(OmniAgeRData::getOmniAgeRData("omniager_causal_clocks_coef", FALSE)[[3]]),
                
                "StocH" = {
                  stocAll <- OmniAgeRData::getOmniAgeRData("omniager_stoch_clocks", verbose)
                  .extractGlmnetWeights(stocAll[["H"]])
                },
                
                "StocP" = {
                  stocAll <- OmniAgeRData::getOmniAgeRData("omniager_stoch_clocks", verbose)
                  .extractGlmnetWeights(stocAll[["P"]])
                },
                
                "StocZ" = {
                  stocAll <- OmniAgeRData::getOmniAgeRData("omniager_stoch_clocks", verbose)
                  .extractGlmnetWeights(stocAll[["Z"]])
                },
    
                "Neu-In" = ctsCoefs[["Neu-In"]],
                "Neu-Sin" = ctsCoefs[["Neu-Sin"]],
                "Glia-In" = ctsCoefs[["Glia-In"]],
                "Glia-Sin" = ctsCoefs[["Glia-Sin"]],
                "Brain" = ctsCoefs[["Brain"]],
                "Hep" = .extractGlmnetWeights(ctsCoefs[["Hep"]]),
                "Liver" = ctsCoefs[["Liver"]],
                
                "StocP" = {
                  stocAll <- OmniAgeRData::getOmniAgeRData("omniager_stoch_clocks", verbose)
                  .extractGlmnetWeights(stocAll[["P"]])
                },
                
                "StocZ" = {
                  stocAll <- OmniAgeRData::getOmniAgeRData("omniager_stoch_clocks", verbose)
                  .extractGlmnetWeights(stocAll[["Z"]])
                },
                
                "scImmuAging" = {
                  scModelData <- OmniAgeRData::getOmniAgeRData("omniager_scimmuaging_model", verbose)
                  scWeightsList <- list()

                  for (ct in names(scModelData$model_set)) {
                    glmObj <- scModelData$model_set[[ct]]
                    if (is.null(glmObj)) next
                    scWeightsList[[ct]] <- .extractGlmnetWeights(glmObj, s = "lambda.min")
                  }
                  
                  return(scWeightsList)
                },
                
                "Brain_CT_Clock" = OmniAgeRData::getOmniAgeRData("omniager_brain_celltype_specific_clocks_coef", verbose),

                "PASTA" = {
                      list(
                      REG = .extractGlmnetWeights(OmniAgeRData::getOmniAgeRData("omniager_cvfit_reg", verbose)),
                      PASTA = .extractGlmnetWeights(OmniAgeRData::getOmniAgeRData("omniager_cvfit_pasta", verbose)),
                      CT46 = .extractGlmnetWeights(OmniAgeRData::getOmniAgeRData("omniager_cvfit_c46", verbose))
                    )
                },
                
                "BohlinGA" = OmniAgeRData::getOmniAgeRData("omniager_bohlin_ga_coef", verbose),
                "EPICGA" = OmniAgeRData::getOmniAgeRData("omniager_epic_ga_coef", verbose),
                "KnightGA" = OmniAgeRData::getOmniAgeRData("omniager_knight_ga_coef", verbose),
                "LeeGA" = OmniAgeRData::getOmniAgeRData("omniager_lee_ga_coef", verbose),
                "MayneGA" = OmniAgeRData::getOmniAgeRData("omniager_mayne_ga_coef", verbose),
                
                "CRP" = OmniAgeRData::getOmniAgeRData("omniager_crp_cpg", verbose),
                "CHIP" = OmniAgeRData::getOmniAgeRData("omniager_chip_cpg", verbose),
                "IL6" = OmniAgeRData::getOmniAgeRData("omniager_il6_coef", verbose),
                "EpiScores" = OmniAgeRData::getOmniAgeRData("omniager_episcores_coef", verbose),
                "McCartney_Trait" = OmniAgeRData::getOmniAgeRData("omniager_mccartney_trait_coef", verbose),
                "SmokeIndex" = OmniAgeRData::getOmniAgeRData("omniager_coeff_smk_idx", verbose),
                "HepatoXuRisk" = OmniAgeRData::getOmniAgeRData("omniager_hepato_xu_coef", verbose),
                
                "EnsembleAge_HumanMouse" = {
                  ensembleAgeAll <- OmniAgeRData::getOmniAgeRData("omniager_ensembleage_coef", verbose)
                  ensembleAgeAll$HumanMouse
                },
                
                "EnsembleAge_Static" = {
                  ensembleAgeAll <- OmniAgeRData::getOmniAgeRData("omniager_ensembleage_coef", verbose)
                  ensembleAgeAll$Static
                },
                
                "EnsembleAge_Dynamic" = {
                  ensembleAgeAll <- OmniAgeRData::getOmniAgeRData("omniager_ensembleage_coef", verbose)
                  ensembleAgeAll$Dynamic
                },
                
                "UniversalPanMammalianClocks" = OmniAgeRData::getOmniAgeRData("omniager_pan_mammalian_clock_coef", verbose),
                "PanMammalianBlood" = OmniAgeRData::getOmniAgeRData("omniager_pan_mammalian_blood_coef", verbose),
                "PanMammalianSkin" = OmniAgeRData::getOmniAgeRData("omniager_pan_mammalian_skin_coef", verbose),
                {
                  NULL
                }
  )
  return(res)
}




#' Extract Non-Zero Weights from glmnet or cv.glmnet Objects
#' 
#' @description 
#' Internal helper to extract the active features, their weights, and the intercept 
#' from a fitted glmnet or cv.glmnet model. It automatically filters out features 
#' with a weight of exactly zero.
#' 
#' @param glmObj A fitted glmnet or cv.glmnet model object.
#' @param s Value of the penalty parameter 'lambda'. For `cv.glmnet`, "lambda.min" 
#'   (default) is used. For standard `glmnet`, it defaults to the last column 
#'   (smallest lambda) in the path.
#' @return A \code{data.frame} with two columns:
#' \describe{
#'   \item{feature}{Character vector containing "(Intercept)" and active features.}
#'   \item{coef}{Numeric vector of the corresponding regression coefficients.}
#' }
#' @keywords internal
.extractGlmnetWeights <- function(glmObj, s = "lambda.min") {
  
  # 1. Handle different glmnet object classes elegantly
  if (inherits(glmObj, "cv.glmnet")) {
    # For CV objects, use the specified string (e.g., "lambda.min")
    fullCoef <- stats::coef(glmObj, s = s)
  } else if (inherits(glmObj, "glmnet")) {
    # For standard glmnet, extract the entire path and take the last column
    # This perfectly mirrors your original stochClocks logic
    allCoefs <- stats::coef(glmObj)
    fullCoef <- allCoefs[, ncol(allCoefs), drop = FALSE]
  } else {
    # Fallback just in case
    fullCoef <- stats::coef(glmObj)
  }
  
  # 2. Extract intercept (first row) and feature weights (remaining rows)
  # Using drop = FALSE above ensures it remains a matrix, so [, 1] is safe
  intercept <- as.numeric(fullCoef[1, 1])
  weightsSparse <- fullCoef[-1, 1] 
  featureNames <- rownames(fullCoef)[-1]
  
  # 3. Filter for active features (non-zero coefficients)
  weightsDense <- as.numeric(weightsSparse)
  keepIdx <- weightsDense != 0
  
  tmpDf <- data.frame(
    feature = "(Intercept)", 
    coef = intercept, 
    stringsAsFactors = FALSE
  )
  
  # 4. Construct the standardized output data frame
  res <- data.frame(
    feature = featureNames[keepIdx],
    coef = weightsDense[keepIdx],
    stringsAsFactors = FALSE
  )
  
  return(rbind(tmpDf, res))
}


#' Standardize Linear Clock Weights Format
#' 
#' @description 
#' Internal helper to ensure all linear clock weights are returned as a 
#' standard data.frame with 'feature' and 'coef' columns.
#' 
#' @param rawCoef A data.frame containing raw clock coefficients.
#' @return A standardized data.frame.
#' @keywords internal
.formatLinearWeights <- function(rawCoef) {
  colnames(rawCoef)[seq_len(2)] <- c("feature", "coef")
  return(rawCoef)
}





#' Marker Category Mapping
#' 
#' @description
#' Defines the hierarchical structure and groups of all supported clocks. 
#' This centralizes the clock registry for easier maintenance.
#'
#' @return A named list where each element is a character vector of clock names.
#' @keywords internal
.getMarkerCategories <- function() {
  mitotic <- c("epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
               "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo")
  dnamtl <- c("DNAmTL", "PCDNAmTL")
  stochastic <- c("StocH", "StocZ", "StocP")
  ensembleAge <- c("EnsembleAge_HumanMouse", 
                   "EnsembleAge_Static", 
                   "EnsembleAge_Dynamic")
  panMammalian <- c("UniversalPanMammalianClocks", "PanMammalianBlood", "PanMammalianSkin")
  cts <- c("Neu-In", "Glia-In", "Neu-Sin", "Glia-Sin", "Hep")
  causal <- c("CausalAge", "DamAge", "AdaptAge")
  transcriptomic <- c("scImmuAging", "Brain_CT_Clock","PASTA")
  
  
  list(
    mitotic = mitotic,
    dnamtl = dnamtl,
    cellularAging = c(mitotic, dnamtl),
    chronological = c(
      "Horvath2013", "Hannum", "Lin", "VidalBralo", "ZhangClock",
      "Horvath2018", "Bernabeu_cAge", "CorticalClock", "PedBE",
      "ABEC", "eABEC", "cABEC", "PipekElasticNet",
      "PipekFilteredh", "PipekRetrainedh", "WuClock",
      "Weidner", "IntrinClock", "Garagnani",
      "CentenarianClock", "Retro_age", "PCHorvath2013", "PCHorvath2018", 
      "PCHannum"),
    biological = c(
      "Zhang10", "PhenoAge", "DunedinPACE", "GrimAge1", "GrimAge2", "PCPhenoAge", 
      "PCGrimAge1", "DNAmFitAge", "IC_Clock", "SystemsAge"),
    causal= causal,  
    stochastic = stochastic,
    cellTypeSpecific = cts, 
    gestationalAge = c("BohlinGA", "EPICGA", "KnightGA", "LeeGA", "MayneGA"),
    surrogateBiomarkers = c("CRP", "CHIP", "IL6", "EpiScores"),
    traitPred = c("McCartney_Trait"),
    diseaseRisk = c("SmokeIndex","HepatoXuRisk"),
    crossSpecies =  c(ensembleAge, panMammalian),
    transcriptomic = transcriptomic,
    
    pcClocks = c("PCHorvath2013", "PCHorvath2018", "PCHannum", "PCPhenoAge", "PCGrimAge1", "PCDNAmTL"),
    stochClocks = stochastic,
    panMammalian = panMammalian,
    ensembleAge = ensembleAge,
    epicmitGroup = c("EpiCMIT_Hyper", "EpiCMIT_Hypo")
  )
}
