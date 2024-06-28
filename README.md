# Drug Combination v1.0

This repository contains the scripts linked to the study "Network-based estimation of therapeutic efficacy and adverse reaction potential for prioritisation of anti-cancer drug combinations". 
The scripts could be used calculate the efficacy and safety estimates of dug combinations followed by developing predictive systems for distinguishing effective from adverse drug combinations. 
It can currently internally compile the data required to develop predictive systems for six cancer types namely breast cancer, kidney cancer, lung cancer, ovary cancer, prostate cancer and skin cancer. 
The repository also contains scripts to validate the predictive system and identify novel drug combinations.

## Running environment

The scripts have been developed and tested on R 4.2.2.
The list of R packages required for running the scripts are listed in the file [Required_R_packages.csv](Environment/Required_R_packages.csv).
The packages including all of its dependencies must be installed prior to execution of the scripts. 
The instructions for installing the CRAN packages can be found [here](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/install.packages) 
and those for BioConductor packages can be found [here](https://www.bioconductor.org/install/). 
