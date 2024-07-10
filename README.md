# Drug Combination v1.0

This repository contains the scripts linked to the development of network medicine approach for identification of cancer-specific drug combinations. The scripts could be used calculate the efficacy and safety estimates of dug combinations followed by developing predictive systems for distinguishing effective from adverse drug combinations. It can currently internally compile the data required to develop predictive systems for six cancer types namely breast cancer, kidney cancer, lung cancer, ovary cancer, prostate cancer and skin cancer. The repository also contains scripts to validate the predictive system and identify novel drug combinations.

## Publications

The following publications refer to this repository:

(1) Network-based estimation of therapeutic efficacy and adverse reaction potential for prioritisation of anti-cancer drug combinations

## Running environment

The scripts have been developed and tested on R 4.2.2. The list of R packages required for running the scripts are listed in the file [Required_R_packages.csv](Environment/Required_R_packages.csv). The packages including all of its dependencies must be installed prior to execution of the scripts. The instructions for installing the CRAN packages can be found [here](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/install.packages) and those for BioConductor packages can be found [here](https://www.bioconductor.org/install/).

## Usage

The scripts could be divided into five parts as listed below:

-   [Part 1: Preparing the data](Manual%20/Part_1_manual.md)
-   [Part 2: Preparing the enrichment libraries](Manual%20/Part_2_manual.md)
-   [Part 3: Feature generation and exploration](Manual%20/Part_3_manual.md)
-   [Part 4: Training and validation of predictive system](Manual%20/Part_4_manual.md)
-   [Part 5: De novo predictions](Manual%20/Part_5_manual.md)

The scripts are detailed in the hyper links in the above list. The scripts are arranged in the way they are supposed to be executed.

## Contacts

-   Arindam Ghosh ([arindam.ghosh\@uef.fi](mailto:arindam.ghosh@uef.fi){.email})
-   Vittorio Fortino ([vittorio.fortino\@uef.fi](mailto:vittorio.fortino@uef.fi){.email})
