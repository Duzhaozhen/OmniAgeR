# OmniAgeR

> [!IMPORTANT]
> **This repository contains the version of OmniAgeR being developed for
> Bioconductor submission.**
>
> In this Bioconductor version, the software functions and model data are
> distributed separately as `OmniAgeR` and `OmniAgeRData`.
>
> Users following the publication
> *OmniAge: a compendium of aging omic biomarkers links mitotic clocks to
> clonal hematopoiesis and causality*
> should use the complete version available from:
>
> https://github.com/Duzhaozhen/OmniAge

OmniAgeR provides a comprehensive suite of tools for calculating and evaluating
aging-related biomarkers from multi-omics data.

## Installation

### Complete version associated with the OmniAge publication

Users following the publication
*OmniAge: a compendium of aging omic biomarkers links mitotic clocks to
clonal hematopoiesis and causality*
should install the complete version from the `Duzhaozhen/OmniAge` repository:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github(
    "Duzhaozhen/OmniAge",
    subdir = "OmniAgeR"
)
```

This version contains most model coefficients directly within the package and
corresponds to the implementation associated with the OmniAge publication.

### Bioconductor version

The Bioconductor version of OmniAgeR is currently under review. Once it becomes
available through Bioconductor, it can be installed using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("OmniAgeR")
```

The accompanying model and example data are provided through the
`OmniAgeRData` package and Bioconductor ExperimentHub.

### Bioconductor development version

The latest Bioconductor development version can be installed from this GitHub
repository:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("Duzhaozhen-BioC/OmniAgeR")
```

The accompanying data package is maintained separately at:

https://github.com/Duzhaozhen-BioC/OmniAgeRData

## Quick Start

```r
library(OmniAgeR)

# Load the Hannum example dataset.
# This requires access to the accompanying OmniAgeRData resources.
data_list <- loadOmniAgeRdata(
    "omniager_hannum_example"
)

beta_matrix <- data_list[[1]]

# Calculate selected aging biomarkers.
results <- epiMarker(
    betaM = beta_matrix,
    clockNames = c(
        "Horvath2013",
        "Hannum",
        "PhenoAge"
    )
)

head(results)
```

## Tutorials

For detailed information on function usage, parameter specifications and
benchmarking examples, please refer to the package vignette:

- [OmniAgeR: User Guide and Tutorials](vignettes/OmniAgeR.Rmd)

## Acknowledgements and AI Usage Statement

During the development of the *OmniAgeR* package, AI-assisted technologies,
specifically Google Gemini, were used to assist with code refactoring,
structural optimization and the drafting of documentation. All AI-assisted
contributions were supervised, reviewed, edited and validated by the human
authors to ensure scientific accuracy, code security and compliance with
Bioconductor standards.