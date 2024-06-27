set.seed(5081)


# Script to compile the list of priority de novo drug combinations


# Load libraries
library(tidyverse)



# Set the parameters
drug_target_type <- "known"


#####


# Compile the complete list
priority_drugCombs <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  priority_drugCombs[[disease]]  <- read.csv(file = paste0("OutputFiles/DeNovo_data_1/Priority_drug_combinations/priorityDrugCombs_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))
  
}
priority_drugCombs <- bind_rows(priority_drugCombs, .id = "Disease")
priority_drugCombs <- priority_drugCombs[, ! colnames(priority_drugCombs) %in% c("DDI_description", "efficacy_score", "safety_score", "final_score", "isSelected_byDiffOfMax")]

if(!dir.exists("OutputFiles/Tables_publication/")){ dir.create("OutputFiles/Tables_publication/", recursive = TRUE) }
write.csv(priority_drugCombs, "OutputFiles/Tables_publication/Compiled_priorityDrugCombs_NES_combinedEfficacySafety.csv", row.names = FALSE)


#####


# Compile a shorter list
priority_drugCombs <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  priority_drugCombs[[disease]]  <- read.csv(file = paste0("OutputFiles/DeNovo_data_1/Priority_drug_combinations/priorityDrugCombs_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))
  
}

priority_drugCombs <- lapply(priority_drugCombs, function(x){ 
  x <- x[order(x$diff_of_maxs, decreasing = TRUE), ] 
  row.names(x) <- NULL
  x[1:3,]
  
  })

priority_drugCombs <- bind_rows(priority_drugCombs, .id = "Disease")

priority_drugCombs <- priority_drugCombs[, c("Disease", "Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                             "Drug1_Year", "Drug1_Product", "Drug1_Indications", 
                                             "Drug2_Year", "Drug2_Product", "Drug2_Indications", 
                                             "which_max_efficacy_NES", "which_max_safety_NES")]

if(!dir.exists("OutputFiles/Tables_publication/")){ dir.create("OutputFiles/Tables_publication/", recursive = TRUE) }
write.csv(priority_drugCombs, "OutputFiles/Tables_publication/shortList_priorityDrugCombs_NES_combinedEfficacySafety.csv", row.names = FALSE)



print(warnings())