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
#' Hannum, G. et al. (2013).
#' \emph{Genome-wide methylation profiles reveal quantitative views of human
#' aging rates}. \emph{Molecular Cell}, 49, 359--367.
#' \doi{10.1016/j.molcel.2012.10.016}
#'
#' Original data: NCBI Gene Expression Omnibus, accession GSE40279:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279}
#'
#' The preprocessing procedure followed Luo, Q. et al. (2023),
#' \emph{Genome Medicine}, 15, 59.
#' \doi{10.1186/s13073-023-01211-5}
#'
#' \strong{Data-use terms:}
#' NCBI places no restrictions on the use or distribution of data deposited in
#' GEO, while noting that submitters may retain applicable patent, copyright,
#' or other intellectual-property rights. See the GEO disclaimer:
#' \url{https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html}.
#'
#'
#' \strong{License:}
#' The original data was deposited in GEO as a public resource. 
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
#' The original dataset is derived from:
#' \itemize{
#'   \item \strong{Title:} Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease
#'   \item \strong{Authors:} Seyhan Yazar et al. (\emph{Science}, 2022)
#'   \item \strong{Source Portal:} CZ CELLxGENE Discover (\url{https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1})
#'   \item \strong{License:} Creative Commons Attribution 4.0 International (CC BY 4.0) (\url{https://creativecommons.org/licenses/by/4.0/})
#' }
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
#' Fröhlich, A. S. et al. (2024).
#' \emph{Single-nucleus transcriptomic profiling of human orbitofrontal cortex
#' reveals convergent effects of aging and psychiatric disease}.
#' \emph{Nature Neuroscience}, 27, 2021--2032.
#' \doi{10.1038/s41593-024-01742-z}
#'
#' Original snRNA-seq data, including the processed AnnData object, are
#' available from NCBI GEO under accession GSE254569:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254569}
#'
#' \strong{Data-use terms:}
#' NCBI places no restrictions on the use or distribution of data deposited in
#' GEO, while noting that submitters may retain applicable patent, copyright,
#' or other intellectual-property rights. See the GEO disclaimer:
#' \url{https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html}.
#'
#' \strong{License:}
#' The original data was deposited in GEO as a public resource. 
#'
"ScBrainExample"