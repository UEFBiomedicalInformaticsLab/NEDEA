set.seed(5081)



# Compile the consensus feature importance across diseases for selected feature type,
# model and data balance methods


# Load libraries
library(openxlsx)
library(tidyverse)

source("Scripts/Functions/Functions_data_manipulation.R")



# Read all files with consensus ranking of features across folds
files <- list.files(path = "OutputFiles/Tables",
                    pattern = "^consensusVarImp_models_",
                    recursive = TRUE,
                    ignore.case = TRUE, full.names = TRUE)

consensusVarImp <- list()
 
for(file in files){
  
  tmp1 <- strsplit(x = file, split = "\\/")[[1]][5]
  tmp2 <- strsplit(x = tmp1, split = "_")[[1]][3]
  tmp3 <- strsplit(x = tmp1, split = "_")[[1]][4]
  tmp3 <- strsplit(x = tmp3, split = "\\.")[[1]][1]
  
  sheet_names <- getSheetNames(file)
  
  for(sheet in sheet_names){
    
    tmp4 <- paste(tmp2, tmp3, sheet, sep = ".")
    consensusVarImp[[tmp4]] <- read.xlsx(file, sheet = sheet)
  }
}
rm(list = c("tmp1", "tmp2", "tmp3", "tmp4"))

feature_list <- unique(unlist(lapply(strsplit(x = names(consensusVarImp), split = "\\."), function(x)x[3])))



if(!dir.exists(paste0("OutputFiles/Tables/featureImportance_xDis/"))){
  dir.create(paste0("OutputFiles/Tables/featureImportance_xDis/"), recursive = TRUE)
}


# featureType <- "DrgDisAdr"
# data_balance_method <- "none"
# model <- "rf"


# Compile the feature importace per feature type, model, data balance method
# for(featureType in c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene", "keggPath")){
for(featureType in feature_list){
  for(data_balance_method in c("SMOTE", "downSample", "upSample", "none")){
    for(model in c("rf", "glmnet", "svmRd", "knn", "nb")){
      
      print(paste(featureType, data_balance_method, model, sep = "__"))
      
      selected_consensusVarImp <- names(consensusVarImp)[grep(featureType, names(consensusVarImp))]
      selected_consensusVarImp <- selected_consensusVarImp[grep(data_balance_method, selected_consensusVarImp)]
      selected_consensusVarImp <- selected_consensusVarImp[grep(model, selected_consensusVarImp)]
      if(length(selected_consensusVarImp) != 0){
              selected_consensusVarImp <- consensusVarImp[selected_consensusVarImp]
              selected_consensusVarImp <- lapply(selected_consensusVarImp, function(x){x[,grep("_median$", colnames(x))]})
              selected_consensusVarImp <- lapply(selected_consensusVarImp, function(x){x[order(x$Scores_median),]})
              names(selected_consensusVarImp) <- sapply(names(selected_consensusVarImp), function(x){strsplit(x, split = "\\.")[[1]][2]})

              for(disease in names(selected_consensusVarImp)){
                colnames(selected_consensusVarImp[[disease]]) <- paste(disease, colnames(selected_consensusVarImp[[disease]]), sep = ".")
              }

              selected_consensusVarImp <- do.call(cbind.fill, selected_consensusVarImp)

              write.table(selected_consensusVarImp, 
                          paste0("OutputFiles/Tables/featureImportance_xDis/", featureType, "_", model, "_", data_balance_method, ".tsv"), 
                          sep = "\t", row.names = FALSE)
      }

    }
  }
}


print(warnings())