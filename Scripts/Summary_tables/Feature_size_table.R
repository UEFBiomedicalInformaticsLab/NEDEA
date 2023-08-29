set.seed(5081)



# Script to compile the library terms and sizes used for enrichment



# Load libraries
library(tidyverse)
library(openxlsx)

library_size <- list()


# Compile efficacy libraries
diseases <- c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer", "LiverCancer")

for(disease in diseases){
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "Size"
  size <- rownames_to_column(size, "Description")
  library_size[[paste0("Efficacy_", disease)]] <- size
}



# Compile safety library
enrichment_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["Safety"]] <- size



# Compile KEGG library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/CHG_keggPath2Gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["KEGG"]] <- size



# Compile SMPDB (Drug Metabolsim) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Metabolism`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["SMPDB_DrugMetabolism"]] <- size



# Compile SMPDB (Drug Action) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Action`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["SMPDB_DrugAction"]] <- size



# Compile miscelleneous library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/miscellaneous_gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["miscellaneous"]] <- size



# Export to file
if(!dir.exists("OutputFiles/Tables/")){
  dir.create("OutputFiles/Tables/", recursive = TRUE)
}
write.xlsx(library_size, "OutputFiles/Tables/Enrichment_lib_size.xlsx")



print(warnings())