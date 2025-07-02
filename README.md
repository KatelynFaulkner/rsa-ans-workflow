# An automated workflow to standardise taxon names for South African alien species lists
This automated workflow was developed to standardise the names of alien taxa for South African alien species lists. It can be used to standardise digitised lists of taxon names that are in one of two formats: (1) full scientific name with authorship and date information (hereafter ‘scientific name’); and (2) scientific name without authorship and date information (hereafter ‘canonical name’). While the workflow can be used to taxonomically standardise canonical names, it is recommended that information on authorship and date is retained as this decreases errors and makes it easier to resolve issues. The workflow implements taxonomic backbones that are relevant to the South African context.

Three taxonomic backbones are implemented: that of the Global Biodiversity Information Facility (GBIF); that of the Botanical Database of Southern Africa (BODATSA); and that of The World Checklist of Vascular Plants (WCVP), accessed via Plants of the World Online (POWO) (https://powo.science.kew.org/).

The taxonomic backbone used depends on the taxonomic group:

1. The GBIF backbone is used for alien taxa that are not vascular plants (e.g., animals, fungi, and chromista)

2. The BODATSA backbone is used for alien vascular plants recorded outside of captivity or cultivation in South Africa (i.e., those taxa included in BODATSA)

3. The WCVP backbone is used for alien vascular plants not found outside of captivity or cultivation in South Africa (i.e., those taxa not included in BODATSA)

In addition to standardising a list of taxon names based on the three taxonomic backbones, the workflow obtains canonical names and higher taxonomic information for the taxa from GBIF, and detects, corrects, and flags potential issues and errors. 

The workflow has been created in R software, version 4.4.0 (R Core Team 2024). 

## Requirements 

The following are required to execute the workflow: 

1. Installed R software (version 4.4.0 or higher) and Rstudio

2. Installed R packages: "tidyr", "dplyr","rgbif", "stringr", "stringdist", "rWCVP", "purrr", "lubridate", "remotes", "rWCVPdata"

3. A stable internet connection

### The R environment

R (version 4.0.0 or higher) and Rstudio need to be installed, which are freely available at: https://cran.r-project.org and https://www.rstudio.com/

### Installing the required packages

Ten R packages and their dependencies must be installed. Nine of these packages can be obtained through the R CRAN, and the remaining package "rWCVPdata" can be obtained from GitHub, using the package "remotes". 

If executed the following code will load and, if required, install the packages.

Specify the packages required from R CRAN:

```{r}
packages = c("tidyr", "dplyr",
             "rgbif", "stringr", "stringdist", "rWCVP", "purrr", "lubridate")
```

Install R CRAN packages (if required) and load:

```{r}
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
```

Install (if required) the "remotes" package from R CRAN and the "rWCVPdata" package from GitHub, and load these packages as required:

```{r}
rWCVPcheck<-"rWCVPdata" %in% rownames(installed.packages())
if(rWCVPcheck == FALSE){
remotes.check<-"remotes" %in% rownames(installed.packages()) 
if(remotes.check == FALSE){
  install.packages("remotes") # install remotes if required
}
remotes::install_github('matildabrown/rWCVPdata')
}
require('rWCVPdata')
```

### Installing and running the scripts

This repository can be downloaded onto your computer as a zip file. The downloaded zip file will contain all the required folders. The zip file will need to be extracted. The main folder contains the Rstudio project file (rsa-ans-workflow.Rproj) of the workflow. This file should be opened to use the workflow. The subfolder `R/` contains the script an annotated .Rmd file. The folder `manuals/` contains a manual describing the workflow (.Rmd), and a README.Rmd that decribes the outputs. The manual must be read before attempting to run the workflow script.
