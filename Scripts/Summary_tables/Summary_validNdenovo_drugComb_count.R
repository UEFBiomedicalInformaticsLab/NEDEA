set.seed(5081)



# Script to summarise the number of drug combinations in validation and de novo data set


# Load libraries
library(tidyverse)


#####


# Read the Validation data

validation_drugCombs <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  if(file.exists(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))){
    validation_drugCombs[["Validation1"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_1a/drugCombs_validation1a_", disease, ".rds"))){
    validation_drugCombs[["Validation1a"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_1a/drugCombs_validation1a_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_1b/drugCombs_validation1b_", disease, ".rds"))){
    validation_drugCombs[["Validation1b"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_1b/drugCombs_validation1b_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))){
    validation_drugCombs[["Validation2"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_2a/drugCombs_validation2a_", disease, ".rds"))){
    validation_drugCombs[["Validation2a"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_2a/drugCombs_validation2a_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_2b/drugCombs_validation2b_", disease, ".rds"))){
    validation_drugCombs[["Validation2b"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_2b/drugCombs_validation2b_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))){
    validation_drugCombs[["Validation3"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_3a/drugCombs_validation3a_", disease, ".rds"))){
    validation_drugCombs[["Validation3a"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_3a/drugCombs_validation3a_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_3b/drugCombs_validation3b_", disease, ".rds"))){
    validation_drugCombs[["Validation3b"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_3b/drugCombs_validation3b_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_4/drugCombs_validation4_", disease, ".rds"))){
    validation_drugCombs[["Validation4"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_4/drugCombs_validation4_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_6/drugCombs_validation6_", disease, ".rds"))){
    validation_drugCombs[["Validation6"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_6/drugCombs_validation6_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_6a/drugCombs_validation6a_", disease, ".rds"))){
    validation_drugCombs[["Validation6a"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_6a/drugCombs_validation6a_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_6b/drugCombs_validation6b_", disease, ".rds"))){
    validation_drugCombs[["Validation6b"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_6b/drugCombs_validation6b_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_7/drugCombs_validation7_", disease, ".rds"))){
    validation_drugCombs[["Validation7"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_7/drugCombs_validation7_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_7a/drugCombs_validation7a_", disease, ".rds"))){
    validation_drugCombs[["Validation7a"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_7a/drugCombs_validation7a_", disease, ".rds"))
  }
  
  if(file.exists(paste0("InputFiles/Validation_data_7b/drugCombs_validation7b_", disease, ".rds"))){
    validation_drugCombs[["Validation7b"]][[disease]] <- readRDS(paste0("InputFiles/Validation_data_7b/drugCombs_validation7b_", disease, ".rds"))
  }
  
}

validation_drugCombs <- unlist(validation_drugCombs, recursive = FALSE)
validation_drugCombs <- lapply(validation_drugCombs, function(x){ as.data.frame(table(x$class_EffAdv, useNA = "ifany")) })
validation_drugCombs <- bind_rows(validation_drugCombs, .id = "Disease")
validation_drugCombs <- separate(data = validation_drugCombs, col = "Disease", into = c("DataSet", "Disease"))
validation_drugCombs <- pivot_wider(validation_drugCombs, names_from = Var1, values_from = Freq)


#####


# Read the de Novo data

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


#####


# Merge data and save counts

drugCombs <- bind_rows(validation_drugCombs, denovo_drugCombs)

drugCombs <- pivot_wider(data = drugCombs,
                         names_from = c("Disease"), 
                         values_from = c("Eff", "Adv", "Unk"), 
                         names_glue = "{Disease}_{.value}")

drugCombs <- column_to_rownames(drugCombs, "DataSet")
drugCombs <- drugCombs[,sort(colnames(drugCombs))]

write.csv(drugCombs, "OutputFiles/Tables/Summary_validNdenovo_drugComb_count.csv", row.names = TRUE)


#####


print(warnings())