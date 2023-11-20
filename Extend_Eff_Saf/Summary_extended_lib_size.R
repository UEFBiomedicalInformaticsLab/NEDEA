

library(igraph)

# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))



summary_df <- data.frame()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the libraries and extract gene list
  enrichment_lib_efficacy <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  enrichment_lib_efficacy <- unique(unlist(enrichment_lib_efficacy, use.names = FALSE))
  
  enrichment_lib_safety <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
  enrichment_lib_safety <- unique(unlist(enrichment_lib_safety, use.names = FALSE))
  
  
  
  extended_enrichment_lib <- readRDS(paste0("Extend_Eff_Saf/Enrichment_analysis_libraries/", disease, "_extended_EfficacySafety_lib.rds"))
  
  
  # Generate summary
  summary_df <- rbind(summary_df,
                      data.frame("Disease"  = disease, 
                                 "Efficacy_seed" = length(enrichment_lib_efficacy),
                                 "Safety_seed" = length(enrichment_lib_safety),
                                 "Efficacy_seed_inNet" = length(enrichment_lib_efficacy[enrichment_lib_efficacy %in% V(rwr_input_network)$name]),
                                 "Safety_seed_inNet" = length(enrichment_lib_safety[enrichment_lib_safety %in% V(rwr_input_network)$name]),
                                 "Efficacy_extended" = length(extended_enrichment_lib$Efficacy), 
                                 "Safety_extended" = length(extended_enrichment_lib$Safety),
                                 "Efficacy_common" = length(intersect(enrichment_lib_efficacy, extended_enrichment_lib$Efficacy)),
                                 "Safety_common" = length(intersect(enrichment_lib_safety, extended_enrichment_lib$Safety)),
                                 "Efficacy_Safety_seed_common" = length(intersect(enrichment_lib_efficacy, enrichment_lib_safety)),
                                 "Efficacy_Safety_extended_common" = length(intersect(extended_enrichment_lib$Efficacy, extended_enrichment_lib$Safety)))
  )
}



