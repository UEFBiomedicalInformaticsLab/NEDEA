set.seed(5081)



# Script to summarise the number of drug combinations in validation and de novo data set


# Load libraries
library(tidyverse)


# Validation data
validation_drugCombs <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  if(file.exists(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))){
    validation_drugCombs[["Validation1"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))){
    validation_drugCombs[["Validation2"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))){
    validation_drugCombs[["Validation3"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_4/drugCombs_validation4_", disease, ".rds"))){
    validation_drugCombs[["Validation4"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_4/drugCombs_validation4_", disease, ".rds"))
  }
  
}

validation_drugCombs <- unlist(validation_drugCombs, recursive = FALSE)
validation_drugCombs <- lapply(validation_drugCombs, function(x){ as.data.frame(table(x$class_EffAdv, useNA = "ifany")) })
validation_drugCombs <- bind_rows(validation_drugCombs, .id = "Disease")
validation_drugCombs <- separate(data = validation_drugCombs, col = "Disease", into = c("DataSet", "Disease"))
validation_drugCombs <- pivot_wider(validation_drugCombs, names_from = Var1, values_from = Freq)



# De Novo data

denovo_drugCombs <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  if(file.exists(paste0("InputFiles/DeNovo_data_1/drugCombs_denovo1_", disease, ".rds"))){
    denovo_drugCombs[["DeNovo1"]][[disease]] <- readRDS(paste0("InputFiles/DeNovo_data_1/drugCombs_denovo1_", disease, ".rds"))
  }
  
}

denovo_drugCombs <- unlist(denovo_drugCombs, recursive = FALSE)
denovo_drugCombs <- lapply(denovo_drugCombs, function(x){ data.frame("Unk" = nrow(x)) })
denovo_drugCombs <- bind_rows(denovo_drugCombs, .id = "Disease")
denovo_drugCombs <- separate(data = denovo_drugCombs, col = "Disease", into = c("DataSet", "Disease"))

drugCombs <- bind_rows(validation_drugCombs, denovo_drugCombs)


write.csv(drugCombs, "OutputFiles/Tables/Summary_validNdenovo_drugComb_count.csv", row.names = FALSE)


print(warnings())