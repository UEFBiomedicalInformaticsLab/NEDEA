set.seed(5081)


# Prepare data for deepDDI

# Load libraries
library(tidyverse)


#####

if(!dir.exists("InputFiles/Compare_external/deepDDI")){
  dir.create("InputFiles/Compare_external/deepDDI", recursive = TRUE)
}

# Get the drug structures as SMILES

# Read the drug smiles
DrugBank_drug_smile <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_smile <- DrugBank_drug_smile$drugs$calculated_properties
DrugBank_drug_smile <- DrugBank_drug_smile[DrugBank_drug_smile$kind %in% "SMILES", ]


for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "comb_name", "class_EffAdv")]
  
  drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
  drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
  drugCombs_targets <- drugCombs_targets[drugCombs_targets$comb_name %in% drugCombs_cat$comb_name, ]
  
  drugCombs_targets <- drugCombs_targets %>% 
    select(c(Drug1_DrugBank_id, Drug2_DrugBank_id)) %>% 
    left_join(DrugBank_drug_smile[, c("parent_key", "value")], by = c("Drug1_DrugBank_id" = "parent_key")) %>% 
    rename("Drug1_SMILES" = "value") %>% 
    left_join(DrugBank_drug_smile[, c("parent_key", "value")], by = c("Drug2_DrugBank_id" = "parent_key")) %>% 
    rename("Drug2_SMILES" = "value") %>%
    select(c(Drug1_DrugBank_id, Drug1_SMILES, Drug2_DrugBank_id, Drug2_SMILES)) %>%
    head(10)
  
  write.table(drugCombs_targets, file = paste0("InputFiles/Compare_external/deepDDI/deepDDI_DrugCombSMILES_", disease, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}



#####


print(warnings())