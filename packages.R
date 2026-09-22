#!/usr/bin/env Rscript

# CRAN packages,
# add any missing/required
# CRAN packages to the list
# directly below. This script
# will install them if they
# are not already installed.
# Install via: 
#   install.packages('packageName', repos='http://cran.r-project.org')
cran_packages <- c(
    "BiocManager",
    "corrplot",
    "ggpmisc",
    "ggpubr",
    "ggh4x",
    "here",
    "lmerTest",
    "mixtools",
    "readxl",
    "SomaDataIO",
    "tidyverse"
)

# Install missing CRAN packages 
install.packages(
    setdiff(
        cran_packages, rownames(installed.packages())
    ),
    repos = 'http://cran.r-project.org'
)

# Bioconductor packages,
# add any missing/required
# Bioconductor packages to the list
# directly below. This script will
# install them if they are missing.
# Install via: 
#   BiocManager::install('packageName')
bioc_packages <- c(
    "Biobase",
    "limma",
    "variancePartition"
)

# Install missing Bioconductor packages 
for (p in bioc_packages){
    if (!p %in% rownames(installed.packages())){
        # Package missing
        print(paste("Installing missing Bioconductor package...", p))
        BiocManager::install(p, ask = FALSE, update = FALSE)
    }
}
