MQ-DATAMIND ECR workshop: Proteomics and Mental Health (environment set
up)
================
X Shen
29 September, 2025

Welcome! This tutorial guides through the main steps to set up your
local environment for the MQ-DATAMIND ECR workshop, ‘Proteomics and
Mental health’.

## Introduction

All scripts are written in R. You can choose to run the practical on a
jupyter notebook provided to you and run the analysis on Google Cloud,
or you could download scripts to your local computer or cluster and
perform analysis on your own device. This Setup tutorial shows you steps
to set up or local environment.

## Download R and RStudio

R and RStudio are available at the URL:
<https://posit.co/download/rstudio-desktop/>

In this page, scroll down and you will find two links for R and RStudio
respectively. Download and install the latest versions of both following
the prompts.

![alt text](image.png)

## Install R packages

Open RStudio, run the following commands in the console:

``` r
check_and_install_package <- function(package_name) {
  # Check if the package is already installed.
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    # If not installed, try to install it.
    message(paste("Package '", package_name, "' not found. Attempting to install...", sep = ""))
    
    # Use tryCatch to handle potential installation errors.
    tryCatch({
      if (package_name!='TwoSampleMR'){
        install.packages(package_name, dependencies = TRUE)
      }else if(package_name=='TwoSampleMR'){
        library(remotes)
        remotes::install_github("MRCIEU/")
      }else if(package_name=='coloc'){
        library(remotes)
        remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
      }
      
      # Try to load the package after installation.
      if (require(package_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Package '", package_name, "' installed and loaded successfully.", sep = ""))
      } else {
        # This case is rare, but good to have a fallback.
        stop("Package installed, but could not be loaded. Please check your R environment.")
      }
    }, error = function(e) {
      message(paste("Error installing package '", package_name, "': ", e$message, sep = ""))
    })
  } else {
    # If the package is already installed, just load it.
    message(paste("Package '", package_name, "' is already installed and loaded.", sep = ""))
  }
}

check_and_install_package("dplyr")
check_and_install_package("data.table")
check_and_install_package("readr")
check_and_install_package("here")
check_and_install_package("remotes")
check_and_install_package("TwoSampleMR")
check_and_install_package("coloc")
```

The code should be able to check if you have installed the packages
correctly. Contact Xueyi Shen (<xueyi.shen@ed.ac.uk>) if you run into
any issue.

## Readings

Optional readings on the datasets we will use:

Readings:

-   Schizophrenia GWAS: Eldjarn, G.H., Ferkingstad, E., Lund, S.H. et
    al. Large-scale plasma proteomics comparisons through genetics and
    disease associations. Nature 622, 348–358 (2023).
    <https://doi.org/10.1038/s41586-023-06563-x>

-   Trubetskoy, V., Pardiñas, A.F., Qi, T. et al. Mapping genomic loci
    implicates genes and synaptic biology in schizophrenia. Nature 604,
    502–508 (2022). <https://doi.org/10.1038/s41586-022-04434-5>

## What next

Additional scripts and mock data will become available prior to the
workshop.
