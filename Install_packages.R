# Install packages

cran_pkgs <- c("caret", "dbparser", "dnet", "doParallel", "foreach", "glmnet", "httr", "igraph", "jsonlite", "MLmetrics", "msigdbr", "naivebayes", "openxlsx", "optparse", "RSQLite", "sparklyr", "sparklyr.nested", "pheatmap", "tidyverse", "TopKLists", "xml2")
bioc_pkgs <- c("XVector", "AnnotationDbi", "BiocParallel", "biomaRt", "Biostrings", "fgsea", "org.Hs.eg.db", "Rgraphviz", "supraHex")

cran_pkgs_install <- cran_pkgs[!(cran_pkgs %in% row.names(installed.packages()))]
bioc_pkgs_install <- bioc_pkgs[!(bioc_pkgs %in% row.names(installed.packages()))]



# CRAN packages
if(length(cran_pkgs_install) > 0 ){
  cat(paste0("\n\nInstalling packages: ", paste0(cran_pkgs_install, collapse = ", "), " \n\n"))
  install.packages(cran_pkgs_install, dependencies = TRUE, repos = "https://cloud.r-project.org/")
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



# 
# 
# install.packages("doParallel", repos="http://R-Forge.R-project.org")
# 
# install.packages("unixtools","http://rforge.net/", type = "source")
