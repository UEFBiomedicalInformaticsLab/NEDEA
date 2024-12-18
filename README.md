# NEDEA: Network-based Evaluation of Drug Efficacy and Adverse effect

This repository contains the scripts linked to the development of network medicine approach for identification of cancer-specific drug combinations. The scripts could be used calculate the efficacy and safety estimates of dug combinations followed by developing predictive systems for distinguishing effective from adverse drug combinations. It can currently internally compile the data required to develop predictive systems for six cancer types namely breast cancer, kidney cancer, lung cancer, ovary cancer, prostate cancer and skin cancer. The repository also contains scripts to validate the predictive system and identify novel drug combinations.

## Publications

The following publications refer to this repository:

(1) [Network-based estimation of therapeutic efficacy and adverse reaction potential for prioritisation of anti-cancer drug combinations](https://doi.org/10.1101/2024.09.17.613439)

## Running environment

The scripts have been developed and tested on R 4.2.2. The list of R packages required for running the scripts are listed in the file [Required_R_packages.csv](Environment/Required_R_packages.csv). The packages including all of its dependencies must be installed prior to execution of the scripts. The instructions for installing the CRAN packages can be found [here](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/install.packages) and those for BioConductor packages can be found [here](https://www.bioconductor.org/install/).

## Usage

The scripts could be divided into five parts as listed below:

-   [Part 1: Preparing the data](Manual/Part_1_manual.md)
-   [Part 2: Preparing the enrichment libraries](Manual/Part_2_manual.md)
-   [Part 3: Feature generation and exploration](Manual/Part_3_manual.md)
-   [Part 4: Training and validation of predictive system](Manual/Part_4_manual.md)
-   [Part 5: De novo predictions](Manual/Part_5_manual.md)

The scripts are detailed in the hyper links in the above list. The scripts are arranged in the way they are supposed to be executed.
The scripts are expected to be executed in the base directory which contains the 'Script' directory. 
Three main directories would be generated while executing the scripts. The 'Databases' directory to store the data downloaded from various online databases, 
the 'InputFiles' directory to store the primary input files and the 'OutputFiles' directory to store all results.

## Contacts

-   Arindam Ghosh (arindam.ghosh@uef.fi)
-   Vittorio Fortino (vittorio.fortino@uef.fi)
