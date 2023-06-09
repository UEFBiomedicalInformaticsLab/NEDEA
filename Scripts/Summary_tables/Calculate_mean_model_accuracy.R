# Script for calculate the mean model accuracy



# Load libraries
library(openxlsx)
library(tidyverse)
library(data.table)





model_stats <- data.frame()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  cat(paste0("\nReading data for: ", disease, "\n"))
  
  # Get the statistics for unbalanced models
  files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                      pattern = "^models_none_[a-zA-Z_]+.xlsx", 
                      ignore.case = TRUE, full.names = TRUE)
  
  none_model_stats <- list()
  
  for(file in files){
    sheet_names_none <- getSheetNames(file)
    
    ## Read files
    for(name in sheet_names_none){
      tmp <- strsplit(x = file, split = "\\/")[[1]][4]
      tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
      none_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
    }
    if(grepl(pattern = "BarabasiProx", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
      names(none_model_stats[[tmp]]) <- paste(prox_comp, names(none_model_stats[[tmp]]), sep = "_")
    }
    if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
      prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
      names(none_model_stats[[tmp]]) <- gsub(pattern = "\\.", replacement = paste0("_", prox_comp, "."), x = names(none_model_stats[[tmp]]), )
    }
  }
  rm(tmp)
  
  none_model_stats <- unlist(none_model_stats, recursive = FALSE)
  none_model_stats <- bind_rows(none_model_stats, .id = "model")
  none_model_stats <- separate(none_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")
  
  none_model_stats <- none_model_stats[, c("featureType", "model", "Fold", 
                                           "PRAUC_train", "PRAUC_test",
                                           "F1_train", "F1_test")]
  
  none_model_stats <- none_model_stats[none_model_stats$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]
  none_model_stats$imbalance <- "none"
  
  
  
  
  
  # Get the statistics for SMOTE models
  files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                      pattern = "^models_SMOTE_[a-zA-Z_]+.xlsx", 
                      ignore.case = TRUE, full.names = TRUE)
  
  smote_model_stats <- list()
  
  for(file in files){
    sheet_names_smote <- getSheetNames(file)
    
    ## Read files
    for(name in sheet_names_smote){
      tmp <- strsplit(x = file, split = "\\/")[[1]][4]
      tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
      smote_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
    }
    if(grepl(pattern = "BarabasiProx", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
      names(smote_model_stats[[tmp]]) <- paste(prox_comp, names(smote_model_stats[[tmp]]), sep = "_")
    }
    if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
      prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
      names(smote_model_stats[[tmp]]) <- gsub(pattern = "\\.", 
                                              replacement = paste0("_", prox_comp, "."), 
                                              x = names(smote_model_stats[[tmp]]), )
    }
  }
  rm(tmp)
  
  smote_model_stats <- unlist(smote_model_stats, recursive = FALSE)
  smote_model_stats <- bind_rows(smote_model_stats, .id = "model")
  smote_model_stats <- separate(smote_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")
  
  smote_model_stats <- smote_model_stats[, c("featureType", "model", "Fold", 
                                             "PRAUC_train", "PRAUC_test",
                                             "F1_train", "F1_test")]
  smote_model_stats <- smote_model_stats[smote_model_stats$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]  
  smote_model_stats$imbalance <- "SMOTE"
  
  
  
  
  
  # Get the statistics for upSample models
  files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                      pattern = "^models_upSample_[a-zA-Z_]+.xlsx", 
                      ignore.case = TRUE, full.names = TRUE)
  
  upSample_model_stats <- list()
  
  for(file in files){
    sheet_names_upSample <- getSheetNames(file)
    
    ## Read files
    for(name in sheet_names_upSample){
      tmp <- strsplit(x = file, split = "\\/")[[1]][4]
      tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
      upSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
    }
    if(grepl(pattern = "BarabasiProx", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
      names(upSample_model_stats[[tmp]]) <- paste(prox_comp, names(upSample_model_stats[[tmp]]), sep = "_")
    }
    if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
      prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
      names(upSample_model_stats[[tmp]]) <- gsub(pattern = "\\.", 
                                                 replacement = paste0("_", prox_comp, "."), 
                                                 x = names(upSample_model_stats[[tmp]]), )
    }
  }
  rm(tmp)
  
  upSample_model_stats <- unlist(upSample_model_stats, recursive = FALSE)
  upSample_model_stats <- bind_rows(upSample_model_stats, .id = "model")
  upSample_model_stats <- separate(upSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")
  
  upSample_model_stats <- upSample_model_stats[, c("featureType", "model", "Fold", 
                                                   "PRAUC_train", "PRAUC_test",
                                                   "F1_train", "F1_test")]
  upSample_model_stats <- upSample_model_stats[upSample_model_stats$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]  
  upSample_model_stats$imbalance <- "upSample"
  
  
  
  
  
  # Get the statistics for downSample models
  files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                      pattern = "^models_downSample_[a-zA-Z_]+.xlsx", 
                      ignore.case = TRUE, full.names = TRUE)
  
  downSample_model_stats <- list()
  
  for(file in files){
    sheet_names_downSample <- getSheetNames(file)
    
    ## Read files
    for(name in sheet_names_downSample){
      tmp <- strsplit(x = file, split = "\\/")[[1]][4]
      tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
      downSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
    }
    if(grepl(pattern = "BarabasiProx", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
      names(downSample_model_stats[[tmp]]) <- paste(prox_comp, names(downSample_model_stats[[tmp]]), sep = "_")
    }
    if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
      prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
      names(downSample_model_stats[[tmp]]) <- gsub(pattern = "\\.", 
                                                   replacement = paste0("_", prox_comp, "."), 
                                                   x = names(downSample_model_stats[[tmp]]), )
    }
  }
  rm(tmp)
  
  downSample_model_stats <- unlist(downSample_model_stats, recursive = FALSE)
  downSample_model_stats <- bind_rows(downSample_model_stats, .id = "model")
  downSample_model_stats <- separate(downSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")
  
  downSample_model_stats <- downSample_model_stats[, c("featureType", "model", "Fold", 
                                                       "PRAUC_train", "PRAUC_test",
                                                       "F1_train", "F1_test")]
  downSample_model_stats <- downSample_model_stats[downSample_model_stats$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]  
  downSample_model_stats$imbalance <- "downSample"
  
  
  
  
  
  # Merge all model stats and rearrange for plotting
  tmp <- rbind(none_model_stats, smote_model_stats, upSample_model_stats, downSample_model_stats) 
  tmp$disease <- disease
  model_stats <- rbind(model_stats, tmp)
  
}






# Perform the group-wise calculations using data.table
model_stats_dt <- as.data.table(model_stats)

mean_accuracy_scores <- data.table()

mean_accuracy_scores <- model_stats_dt[, .(mean_PRAUC_train = round(mean(PRAUC_train, na.rm = TRUE),3),
                                           median_PRAUC_train = round(median(PRAUC_train, na.rm = TRUE),3),
                                           sd_PRAUC_train = round(sd(PRAUC_train, na.rm = TRUE),3),
                                           mean_PRAUC_test = round(mean(PRAUC_test, na.rm = TRUE),3),
                                           median_PRAUC_test = round(median(PRAUC_test, na.rm = TRUE),3),
                                           sd_PRAUC_test = round(sd(PRAUC_test, na.rm = TRUE),3),
                                           mean_F1_train = round(mean(F1_train, na.rm = TRUE),3),
                                           median_F1_train = round(median(F1_train, na.rm = TRUE),3),
                                           sd_F1_train = round(sd(F1_train, na.rm = TRUE),3),
                                           mean_F1_test = round(mean(F1_test, na.rm = TRUE),3),
                                           median_F1_test = round(median(F1_test, na.rm = TRUE),3),
                                           sd_F1_test = round(sd(F1_test, na.rm = TRUE),3)),
                                       by = .(featureType, model, imbalance, disease)]


write.csv(mean_accuracy_scores, "OutputFiles/Tables/panCancer_meanModelAccuracy.csv", row.names = FALSE)



print(warnings())