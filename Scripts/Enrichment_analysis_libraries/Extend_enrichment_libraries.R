set.seed(5081)



# Script to extend the enrichment library


# Load libraries
library(igraph)
source("Scripts/Functions/Functions_extend_enrichment_library.R")



# Read the network which will be used to extend the network
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))


# Extend the disease libraries
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  print(paste0(" -- Extending ", disease, " library"))
  # Read the enrichment library
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  
  # Extend the library
  # extended_enrichment_lib <- func_extendEnrichmentLib_bySteinerTree(enrichment_library = enrichment_lib, 
  #                                                                   input_network = input_network, 
  #                                                                   ST_type = "KB", 
  #                                                                   nproc = 10)
  
  
  extended_enrichment_lib <- func_extendEnrichmentLib_byRWR(enrichment_library = enrichment_lib, 
                                                            input_network = input_network, 
                                                            nproc = 10)
  
  extended_enrichment_lib <- extended_enrichment_lib[names(which(!unlist(lapply(extended_enrichment_lib, anyNA))))]
  
  if(!dir.exists("InputFiles/Enrichment_analysis_libraries_extended/")){dir.create("InputFiles/Enrichment_analysis_libraries_extended/", recursive = TRUE)}
  saveRDS(extended_enrichment_lib, paste0("InputFiles/Enrichment_analysis_libraries_extended/Disease2Gene_", disease, "_extendedLib.rds"))
}


# Extend the safety library
print(paste0(" -- Extending safety library"))
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
# extended_enrichment_lib <- func_extendEnrichmentLib_bySteinerTree(enrichment_library = enrichment_lib, 
#                                                                   input_network = input_network, 
#                                                                   ST_type = "SP", 
#                                                                   nproc = 10)

extended_enrichment_lib <- func_extendEnrichmentLib_byRWR(enrichment_library = enrichment_lib,
                                                          input_network = input_network,
                                                          nproc = 10)
if(!dir.exists("InputFiles/Enrichment_analysis_libraries_extended/")){dir.create("InputFiles/Enrichment_analysis_libraries_extended/", recursive = TRUE)}
saveRDS(extended_enrichment_lib, "InputFiles/Enrichment_analysis_libraries_extended/curatedAdr2Gene_extendedLib.rds")



print(warnings())