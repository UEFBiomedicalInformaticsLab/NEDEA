set.seed(5081)



# Extract all small molecular drug-drug interactions from DrugBank



# Load libraries
library(tidyverse)



if(!dir.exists("InputFiles/ReferenceList/")){
  dir.create("InputFiles/ReferenceList/", recursive = TRUE)
} 

# Adverse drug drug interactions from DrugBank (only small molecules)
DrugBank_drugInteractions <- read.csv("Databases/DrugBank/drug_drug_interactions.csv")
DrugBank_drugInteractions <- DrugBank_drugInteractions[, c("drugbank.id", "parent_key")] 
colnames(DrugBank_drugInteractions) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")

DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", "primary_key"] 

DrugBank_drugInteractions <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in% DrugBank_Drugs & DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in% DrugBank_Drugs), ]
rownames(DrugBank_drugInteractions) <- NULL

saveRDS(DrugBank_drugInteractions, "InputFiles/ReferenceList/DrugBank_drugInteractions.rds")




# Adverse drug drug interactions from DrugBank (only small molecules and involving risk/severity for any side effect)
# Keyword chosen based on literature review;
# Refer: A multimodal deep learning framework for predicting drug–drug interaction events (supplementary material)
DrugBank_drugInteractions <- read.csv("Databases/DrugBank/drug_drug_interactions.csv")
DrugBank_drugInteractions <- DrugBank_drugInteractions[str_detect(DrugBank_drugInteractions$description, "risk|severity"), ]
DrugBank_drugInteractions <- DrugBank_drugInteractions[, c("drugbank.id", "parent_key")] 
colnames(DrugBank_drugInteractions) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")

DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", "primary_key"] 

DrugBank_drugInteractions <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in% DrugBank_Drugs & DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in% DrugBank_Drugs), ]
rownames(DrugBank_drugInteractions) <- NULL

saveRDS(DrugBank_drugInteractions, "InputFiles/ReferenceList/DrugBank_drugInteractions_withRiskSeverity.rds")



print(warnings())