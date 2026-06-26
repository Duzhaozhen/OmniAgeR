#' Example DNA methylation data
#'
#' @description
#' A lightweight example dataset containing whole blood DNA methylation profiles
#' and matched phenotypic metadata from a random subset of the Hannum cohort.
#' It is provided to test DNAm aging clocks in the \code{OmniAgeR} package.
#'
#' @format
#' This resource is a \code{list} of length 2, containing:
#' \describe{
#'   \item{betaM}{A numeric matrix of BMIQ-normalized DNA methylation beta values.}
#'   \item{phenotype}{A \code{data.frame} of 6 phenotypic variables associated with the samples: \code{Age}, \code{Source}, \code{Plate}, \code{Sex}, \code{Ethn} (Ethnicity), and \code{SmkStat} (Smoking Status).}
#' }
#'
#' @details
#' \strong{How it was generated:}
#' The original Hannum cohort dataset (Illumina 450k platform) was downloaded and
#' preprocessed as previously described by Luo et al. (\emph{Genome Med}, 2023).
#' Following BMIQ normalization to correct for probe-type bias, I randomly extracted
#' a subset of 3 samples to create a lightweight testing object. Subsequently, I 
#' retained only the CpGs commonly utilized by the majority of epigenetic clocks. 
#' The beta matrix and clinical variables were then bundled into a standard R list 
#' to seamlessly interface with the \code{OmniAgeR} pipeline.
#'
#' @source
#' The original data is derived from the Hannum et al. study (GEO Accession: GSE40279).
#' Preprocessing pipelines reference: Luo, Q. et al. \emph{Genome Med} 15, 59 (2023).
#'
#' \strong{License:}
#' The original data was deposited in GEO as a public resource. de
"dnamExample"


#' @title Example scRNA-seq Dataset for Transcriptomic Clocks (Yazar Cohort)
#'
#' @description
#' A lightweight single-cell RNA sequencing (scRNA-seq) example dataset
#' containing expression profiles of CD4+ T cells from a subset
#' of the Yazar (OneK1K) cohort. This \code{SingleCellExperiment} object is specifically
#' provided to demonstrate and test cell-type-specific transcriptomic aging
#' clocks within the \code{OmniAgeR} package. 
#'
#' @format
#' A \code{SingleCellExperiment} object containing 3116 features (genes) across 8035 cells
#' (1 assay: "RNA" with "counts" and "data" layers).
#'
#' Key metadata columns include:
#' \describe{
#'   \item{\code{celltype}}{Cell type annotations, specifically subsetted to "CD4T".}
#'   \item{\code{donor_id}}{Anonymized individual donor identifiers.}
#' }
#'
#' @details
#' \strong{How it was generated:}
#' The original, full-scale scRNA-seq dataset was accessed via the CZ CELLxGENE portal.
#' To create a lightweight and functional testing object suitable for the \code{OmniAgeR}
#' pipeline, I computationally subsetted the original data to retain only the CD4T
#' cell populations. Subsequently, I randomly downsampled the cohort to
#' include exactly 5 unique donors. I then retained only the specific genes required for 
#' scImmuAging and PASTA predictions, and finally converted the object into the standard 
#' \code{SingleCellExperiment} format to preserve both raw counts and normalized data layers.
#'
#' @source
#' The original comprehensive dataset is hosted on the CZ CELLxGENE Discover portal.
#' \strong{Collection URL:} \url{https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1}
#'
#' \strong{License:}
#' Data hosted on CELLxGENE are generally distributed under permissive open-access
#' licenses. This highly subsetted and downsampled \code{SingleCellExperiment}
#' object is distributed here strictly for academic testing and reproducibility.
#'
#'
"ScPbmcExample"


#' @title Example scRNA-seq Brain Dataset (Fröhlich Cohort - Oligodendrocytes)
#'
#' @description
#' A lightweight single-cell RNA sequencing (scRNA-seq) example dataset
#' containing expression profiles of oligodendrocytes from 5 healthy
#' control donors. This \code{SingleCellExperiment} object is provided to demonstrate
#' and test cell-type-specific transcriptomic aging clocks
#' within the \code{OmniAgeR} package.
#'
#' @format
#' A \code{SingleCellExperiment} object containing scRNA-seq expression data.
#'
#' Key metadata columns include:
#' \describe{
#'   \item{\code{celltype}}{Cell type annotations, containing exclusively "Oligodendrocytes".}
#'   \item{\code{donor_id}}{Anonymized individual donor identifiers (randomly subsetted to 5 healthy control donors).}
#' }
#'
#' @details
#' \strong{How it was generated:}
#' I obtained the original single-cell dataset from the Gene Expression Omnibus
#' (GEO accession: GSE254569). To create an efficient testing object, I subsetted
#' the data to include only the "Oligodendrocytes" cell type from healthy control
#' individuals, subsequently downsampling it to exactly 5 unique donors. Finally,
#' I converted the original Python-based AnnData (\code{.h5ad}) format into a
#' standard \code{SingleCellExperiment} object to ensure seamless integration with the R
#' environment and the \code{OmniAgeR} analytical pipeline.
#'
#' @source
#' The original dataset is hosted on the Gene Expression Omnibus (GEO).
#' \strong{GEO Accession:} \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254569}
#'
#' \strong{License:}
#' The original data was deposited in GEO as a public resource. 
#'
"ScBrainExample"