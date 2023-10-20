# Install packages

cran_pkgs <- c("caret", "dbparser", "devtools", "dnet", "doParallel", "foreach", "ggpubr", "httr", 
               "igraph", "jsonlite", "msigdbr", "openxlsx", "optparse", "readxl", 
               "sparklyr", "sparklyr.nested", "tidyverse", "unixtools")

bioc_pkgs <- c("BiocParallel", "biomaRt", "fgsea", "org.Hs.eg.db", "RCy3", "Rgraphviz", "supraHex")

cran_pkgs_install <- cran_pkgs[!(cran_pkgs %in% row.names(installed.packages()))]
bioc_pkgs_install <- bioc_pkgs[!(bioc_pkgs %in% row.names(installed.packages()))]



# CRAN packages
if(length(cran_pkgs_install) > 0 ){
  cat(paste0("\n\nInstalling packages: ", paste0(cran_pkgs_install, collapse = ", "), " \n\n"))
  install.packages(cran_pkgs_install, dependencies = TRUE, repos = c("https://cloud.r-project.org/", "http://rforge.net/"))
}



# BioConductor packages
if(length(bioc_pkgs_install) > 0 ){
  cat(paste0("\n\nInstalling packagens: ", paste0(bioc_pkgs_install, collapse = ", "), " \n\n"))
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",  dependencies = TRUE, repos = "https://cloud.r-project.org/")
  BiocManager::install(bioc_pkgs_install)
}


cran_pkgs_install <- cran_pkgs[!(cran_pkgs %in% row.names(installed.packages()))]
bioc_pkgs_install <- bioc_pkgs[!(bioc_pkgs %in% row.names(installed.packages()))]


if(length(c(cran_pkgs_install)) > 0){
  cat(paste0("\n\n\nInstallation failed for CRAN packages: ", paste0(cran_pkgs_install, collapse = ", "), "\n\n"))
}
if(length(c(bioc_pkgs_install)) > 0){
  cat(paste0("\n\n\nInstallation failed for BioConductor packages: ", paste0(bioc_pkgs_install, collapse = ", "), "\n\n"))
}



# Install OmniPath from source
require(devtools)
install_github("saezlab/OmnipathR")
