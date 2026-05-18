## ----mototic-clock-example----------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
lungInv <- loadOmniAgeRdata(
    "omniager_lung_inv",
    verbose = FALSE
)
lungInvM <- lungInv$bmiq_m
phenoDf <- lungInv$PhenoTypes

my_comparisons <- list(c("N\nN=21", "LCIS\nN=13"), c("LCIS\nN=13", "LCIS->LC\nN=22"))
table(phenoDf$Group)
## Check available epigenetic clocks
listEpiMarker()

