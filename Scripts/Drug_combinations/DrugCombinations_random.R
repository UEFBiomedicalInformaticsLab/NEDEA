# Script to generate random drug combinations based on drugs in DrugBank





# Load libraries
library(tidyverse)





# DrugBank drug groups (small molecule only)
DrugBank_drug_groups <- read.csv("Databases/DrugBank/drug_groups.csv")
DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", ] 
DrugBank_drug_groups <- DrugBank_drug_groups[DrugBank_drug_groups$`drugbank.id` %in% DrugBank_Drugs$primary_key, ]




# Read drug combinations used in training set
drugCombs_training <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", disease, ".rds"))
  drugCombs <- do.call(rbind, drugCombs)
  row.names(drugCombs) <- NULL
  drugCombs_training[[disease]] <- drugCombs
}
drugCombs_training <- do.call(rbind, drugCombs_training)
row.names(drugCombs_training) <- NULL
drugCombs_training <- unique(drugCombs_training)
drugs_training <- unique(c(drugCombs_training$Drug1_DrugBank_drug_id, drugCombs_training$Drug2_DrugBank_drug_id))












# Read drugs approved for cancer
# Terms used:
# Cancer = MONDO_0004992
# Carcinoma = EFO_0000313

OpenTargets_Drug_Disease_Net <- readRDS("InputFiles/Associations/OpenTargets_Drug_Disease_Net.rds")
OpenTargets_approvedDrugs_cancer <- OpenTargets_Drug_Disease_Net[OpenTargets_Drug_Disease_Net$Node2_disease_id %in% c("MONDO_0004992", "EFO_0000313"), ]
OpenTargets_approvedDrugs_cancer <- OpenTargets_approvedDrugs_cancer[OpenTargets_approvedDrugs_cancer$Drug_Disease_clinical_status == "Approved",]
OpenTargets_approvedDrugs_cancer <- unique(OpenTargets_approvedDrugs_cancer$Node1_drugbank_drug_id)






# Random combinations 1
tmp1 <- DrugBank_drug_groups[DrugBank_drug_groups$group == "approved",] 
tmp1 <- tmp1[tmp1$`drugbank.id` %in% OpenTargets_approvedDrugs_cancer, ]
tmp1 <- tmp1[!tmp1$`drugbank.id` %in% drugs_training, ]