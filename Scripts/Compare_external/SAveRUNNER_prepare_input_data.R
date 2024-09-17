set.seed(5081)


# Build the input files for SAveRUNNER


# Load libraries
library(org.Hs.eg.db)
library(tidyverse)

#####


# Prepare the interactome

input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
input_network <- igraph::as_data_frame(input_network, "edges")
colnames(input_network) <- c("Gene_A", "Gene_B")


if(!dir.exists("InputFiles/Compare_external/SAveRUNNER")){
  dir.create("InputFiles/Compare_external/SAveRUNNER", recursive = TRUE)
} 
write.table(input_network, file = "InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_Interactome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#####


# Prepare the disease/ADR gene associations


# Compile efficacy and safety libraries
drug_target_type <- "known"

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
  
  enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
  names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
  
  enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)
  
  # Keep only features that were used in the predictive system
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  predictive_system_features <- model$feature_threshold$feature
  enrichment_lib <- enrichment_lib[names(enrichment_lib) %in% predictive_system_features] 
  
  phenotype_gene_link <- bind_rows(lapply(enrichment_lib, as.data.frame), 
                                   .id = "disease") 
  colnames(phenotype_gene_link)[2] <- "GeneID"
  
  write.table(phenotype_gene_link, file = paste0("InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_PhenotypeGenes_", disease, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}


#####


# Prepare the drug combination-target asoications

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "comb_name", "class_EffAdv")]
  
  drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
  drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
  drugCombs_targets <- drugCombs_targets[drugCombs_targets$comb_name %in% drugCombs_cat$comb_name, ]

  drugCombs_targets <- drugCombs_targets %>% 
    select(c(comb_name, drugTarget_ensembl_id)) %>% 
    separate_rows(drugTarget_ensembl_id, sep = ",") %>%
    rename("Drug" = "comb_name", "GeneID" = "drugTarget_ensembl_id")
  
  write.table(drugCombs_targets, file = paste0("InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_DrugCombTargets_", disease, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
}


#####


print(warnings())