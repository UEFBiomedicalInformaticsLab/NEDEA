set.seed(5081)



# Script for plotting the variable importance of the model



# Load libraries
library(optparse)
library(openxlsx)
library(caret)
library(openxlsx)
library(pheatmap)
library(tidyverse)
library(TopKLists)



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: none. Default: none.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$data_balance_method %in% c("none")){
  stop("--data_balance_method: No data balancing methods currently included. Default: none")
}


# Define global options for this script 
disease <- opt$disease
data_balance_method <- opt$data_balance_method



cat(paste0("\n\nExtracting feature importance for: ", disease, "-", data_balance_method, "\n\n"))



# Read models files
files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                    pattern = paste0("models_", data_balance_method, "_[a-zA-Z_]+.rds"), 
                    ignore.case = TRUE, full.names = TRUE)



# Extract feature importance
variable_importance <- list()
for(file in files){
  tmp1 <- strsplit(x = file, split = "\\/")[[1]][4]
  tmp1 <- strsplit(x = tmp1, split = "\\.")[[1]][1]
  
  model <- readRDS(file)
  for(featureType in names(model)){
    for(fold in names(model[[featureType]])) {
      for(select_model in names(model[[featureType]][[fold]])){
        
        tmp2 <- model[[featureType]][[fold]][[select_model]][["model"]] 
        if(class(tmp2) == "train"){
          tmp3 <- varImp(tmp2)
          if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
            tmp4 <- strsplit(x = file, split = "\\/")[[1]][4]
            tmp5 <- strsplit(x = tmp4, split = "\\_")[[1]][4]
            tmp4 <- strsplit(x = tmp4, split = "\\_")[[1]][6]
            tmp4 <- strsplit(x = tmp4, split = "\\.")[[1]][1]
            tmp4 <- paste(tmp5, featureType, tmp4, sep = "_")
          }else if(grepl("BarabasiProx", file)){
            tmp4 <- strsplit(x = file, split = "\\/")[[1]][4]
            tmp4 <- strsplit(x = tmp4, split = "\\_")[[1]][4]
            tmp4 <- paste(tmp4, featureType, sep = "_")
          }else{
            tmp4 <- featureType
          }
          if(ncol(tmp3$importance) == 1){
            variable_importance[[tmp1]][[tmp4]][[select_model]][[fold]] <- tmp3$importance
          }else{
            variable_importance[[tmp1]][[tmp4]][[select_model]][[fold]] <- tmp3$importance[, "Eff", drop = FALSE]
          }
        }
      }
    }
  }
}
rm(list = c("tmp1", "tmp2", "tmp3", "tmp4", "tmp5"))

variable_importance <- unlist(variable_importance, recursive = FALSE)



# Turn variable importance into data frame
variable_importance_df <- list()
for(featureType in names(variable_importance)){
  for (select_model in names(variable_importance[[featureType]])){
    tmp1 <- do.call(rbind, variable_importance[[featureType]][[select_model]])
    tmp1$features <- rownames(tmp1)
    rownames(tmp1) <- NULL
    tmp1 <- separate(data = tmp1, col = "features", into = c("Fold", "Features"), sep = "\\.", extra = "merge")
    tmp1 <- reshape(tmp1, direction = "wide", idvar = "Features", timevar = "Fold")
    colnames(tmp1) <- gsub("^[a-zA-Z]+\\.", "", colnames(tmp1))
    tmp2 <- strsplit(x = featureType, split = "\\.")[[1]][2]
    variable_importance_df[[tmp2]][[select_model]] <- tmp1
  }
}
rm(list = c("tmp1", "tmp2"))

# Export to xlsx file
tmp <- unlist(variable_importance_df, recursive = FALSE)
names(tmp) <- gsub("BbsiProx_", "BBSI_", names(tmp))
names(tmp) <- gsub("DrugDrug", "DrgDrg", names(tmp))
names(tmp) <- gsub("DrugDisease", "DrgDis", names(tmp))
names(tmp) <- gsub("DrugAdr", "DrgAdr", names(tmp))
names(tmp) <- gsub("svmRadial", "svmRd", names(tmp))
names(tmp) <- gsub("closest", "clo", names(tmp))
names(tmp) <- gsub("centre", "cen", names(tmp))
names(tmp) <- gsub("kernel", "ker", names(tmp))
names(tmp) <- gsub("separation", "sep", names(tmp))
names(tmp) <- gsub("shortest", "sho", names(tmp))


if(!dir.exists(paste0("OutputFiles/Tables/", disease, "/featureImportance/"))){
  dir.create(paste0("OutputFiles/Tables/", disease, "/featureImportance/"), recursive = TRUE)
}

write.xlsx(tmp, paste0("OutputFiles/Tables/", disease, "/featureImportance/varImp_models_", data_balance_method, "_", disease, ".xlsx"), 
           overwrite = TRUE)
rm(tmp)



# Create consensus ranking for the features
consensus_rank <- list()
for(featureType in names(variable_importance)){
  for (select_model in names(variable_importance[[featureType]])){
    tmp1 <- variable_importance[[featureType]][[select_model]]
    tmp1 <- lapply(tmp1, function(x){x[order(x[,1], decreasing = TRUE),, drop = FALSE]})
    tmp1 <- lapply(tmp1, row.names)
    tmp2 <- Borda(tmp1)
    for(i in names(tmp2)){
      colnames(tmp2[[i]]) <- paste(i, colnames(tmp2[[i]]), sep = "_")
    }
    tmp2 <- merge(tmp2$TopK, tmp2$Scores, by = 0, all = TRUE)
    tmp2 <- tmp2[, -1]
    tmp2 <- tmp2[order(tmp2$Scores_mean),]
    tmp3 <- strsplit(x = featureType, split = "\\.")[[1]][2]
    consensus_rank[[tmp3]][[select_model]] <- tmp2
  }
}
rm(list = c("tmp1", "tmp2", "tmp3"))

# Export to xlsx file
tmp <- unlist(consensus_rank, recursive = FALSE)
names(tmp) <- gsub("BbsiProx_", "BBSI_", names(tmp))
names(tmp) <- gsub("DrugDrug", "DrgDrg", names(tmp))
names(tmp) <- gsub("DrugDisease", "DrgDis", names(tmp))
names(tmp) <- gsub("DrugAdr", "DrgAdr", names(tmp))
names(tmp) <- gsub("svmRadial", "svmRd", names(tmp))
names(tmp) <- gsub("closest", "clo", names(tmp))
names(tmp) <- gsub("centre", "cen", names(tmp))
names(tmp) <- gsub("kernel", "ker", names(tmp))
names(tmp) <- gsub("separation", "sep", names(tmp))
names(tmp) <- gsub("shortest", "sho", names(tmp))

if(!dir.exists(paste0("OutputFiles/Tables/", disease, "/featureImportance/"))){
  dir.create(paste0("OutputFiles/Tables/", disease, "/featureImportance/"), recursive = TRUE)
}

write.xlsx(tmp, paste0("OutputFiles/Tables/", disease, "/featureImportance/consensusVarImp_models_", data_balance_method, "_", disease, ".xlsx"), 
           overwrite = TRUE)
rm(tmp)



# Create heatmaps
if(!dir.exists(paste0("OutputFiles/Plots/", disease, "/featureImportance/"))){
  dir.create(paste0("OutputFiles/Plots/", disease, "/featureImportance/"), recursive = TRUE)
}
for(featureType in names(variable_importance_df)){
  for (select_model in names(variable_importance_df[[featureType]])){
    # print(paste(disease, featureType, select_model, data_balance_method, sep = " - "))
    tmp1 <- variable_importance_df[[featureType]][[select_model]]
    tmp1 <- tmp1[!rowSums(is.na(tmp1)) >= length(colnames(tmp1)[grep("^Repeat", colnames(tmp1))]), ]
    if(nrow(tmp1) != 0){
      row.names(tmp1) <- NULL
      tmp1 <- column_to_rownames(tmp1, "Features")
      pheatmap(mat = as.matrix(tmp1), 
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               main = paste(disease, featureType, select_model, data_balance_method, sep = " - "),
               filename = paste0("OutputFiles/Plots/", disease, 
                                 "/featureImportance/varImp_models_",  
                                 disease, "_", featureType, "_", select_model, "_", data_balance_method, ".pdf"),
               width = 20,
               height = nrow(tmp1) * 0.5)
    }
    
  }
}



print(warnings())