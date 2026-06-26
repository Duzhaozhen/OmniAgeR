# OmniAgeR

Provides a comprehensive suite of tools for calculating and evaluating various aging biomarkers from multi-omics data

## Installation

**Bioconductor Release Version** (Recommended)
You can install the official release version of OmniAgeR from Bioconductor once it is available:

```r, eval=FALSE
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("OmniAgeR")
```

**Development Version** (GitHub)
You can install the latest development version directly from GitHub:

```r, eval=FALSE
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("Duzhaozhen/OmniAgeR")
```

## 📖 Quick Start

```r, eval=FALSE
library(OmniAgeR)

# 1. Load example data (requires OmniAgeRData package installed)
# Using the Hannum example dataset as a demonstration
data_list <- loadOmniAgeRdata("omniager_hannum_example")
beta_matrix <- data_list[[1]]

# 2. Call the core function to calculate aging scores
# You can specify multiple clocks of interest in the clockNames argument
results <- epiMarker(
    betaM = beta_matrix, 
    clockNames = c("Horvath2013", "Hannum", "PhenoAge")
)

# 3. View the calculation results
head(results)
```

## 📖 Tutorials
For comprehensive details on function usage, parameter specifications, and benchmarking case studies, please refer to the package vignette:
* [OmniAgeR: User Guide and Tutorials](vignettes/OmniAgeR.Rmd) - Comprehensive guide for the R-based workflow.


## Acknowledgements / AI Usage Statement

During the development of the *OmniAgeR* package, AI-assisted technologies (specifically, Google Gemini) were utilized to aid in code refactoring, structural optimization, and the drafting of documentation. All AI-assisted contributions were strictly supervised, thoroughly reviewed, edited, and validated by the human authors to ensure scientific accuracy, code security, and compliance with Bioconductor standards.