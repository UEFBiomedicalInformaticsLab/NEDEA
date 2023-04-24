set.seed(5081)
rm(list = ls())



# Script for plotting the variable importance of the model



# Load libraries
library(optparse)
library(openxlsx)
library(svglite)
library(caret)
library(openxlsx)
library(pheatmap)
library(tidyverse)
library(svglite)






# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--model"), type = "character", default = NULL, 
              help = "The model for which to plot feature importance. Possible values: rf, glmnet, svmRadial, knn, nb", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: SMOTE, downSample, upSample, or none. Default: none.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(is.null(opt$model)){
  print_help(opt_parser)
  stop("--model argument needed", call.=FALSE)
}
# if(is.null(opt$feature_type)){
#   print_help(opt_parser)
#   stop("--feature_type argument needed", call.=FALSE)
# }

if(!opt$model %in% c("rf", "glmnet", "svmRadial", "knn", "nb")){
  print_help(opt_parser)
  stop("Possible values for --model: rf, glmnet, svmRadial, knn, nb", call.=FALSE)
}
# if(!opt$feature_type %in% c("BarabasiProx_DrugDrug", "rwrFgsea", "SteinerTopol_AdrDisDrg")){
#   print_help(opt_parser)
#   stop("Possible values for --feature_type: BarabasiProx_DrugDrug, rwrFgsea, SteinerTopol_AdrDisDrg", call.=FALSE)
# }

if(!opt$data_balance_method %in% c("SMOTE", "downSample", "upSample", "none")){
  stop("--data_balance_method must be SMOTE, downSample, upSample or none")
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
model <- opt$model
# feature_type <- opt$feature_type
data_balance_method <- opt$data_balance_method
nproc <- opt$nproc



cat(paste0("\n\nExtracting feature importance for: ", disease, "-", model, "-", data_balance_method, "\n\n"))



# Read models files
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = paste0("models_", data_balance_method, "_[a-zA-Z_]+.rds"), 
                    ignore.case = TRUE, full.names = TRUE)









print("A")
# Extract feature importance
variable_importance <- list()
for(file in files){
  print(file)
  tmp1 <- strsplit(x = file, split = "\\/")[[1]][5]
  tmp1 <- strsplit(x = tmp1, split = "\\.")[[1]][1]
  
  model <- readRDS(file)
  for(featureType in names(model)){
    for(fold in names(model[[featureType]])) {
      for(select_model in names(model[[featureType]][[fold]])){
        
        print(paste(featureType, fold, select_model, sep = " - "))
        tmp2 <- model[[featureType]][[fold]][[select_model]][["model"]] 
        if(!is.na(tmp2)){
          tmp3 <- varImp(tmp2)
          if(grepl("BarabasiProx", file)){
            tmp4 <- strsplit(x = file, split = "\\/")[[1]][5]
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
rm(list = c("tmp1", "tmp2", "tmp3", "tmp4"))

variable_importance <- unlist(variable_importance, recursive = FALSE)


print("B")

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




print("C")

# Export to xlsx file
tmp <- unlist(variable_importance_df, recursive = FALSE)
names(tmp) <- gsub("BbsiProx_", "BBSI_", names(tmp))
names(tmp) <- gsub("DrugDrug", "DrgDrg", names(tmp))
names(tmp) <- gsub("DrugDisease", "DrgDis", names(tmp))
names(tmp) <- gsub("DrugAdr", "DrgAdr", names(tmp))
names(tmp) <- gsub("svmRadial", "svmRd", names(tmp))

if(!dir.exists(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/featureImportance/"))){
  dir.create(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/featureImportance/"), recursive = TRUE)
}

write.xlsx(tmp, paste0("Analysis/STRING/DrugCombs_v5/", disease, "/featureImportance/varImp_models_", data_balance_method, "_", disease, ".xlsx"), 
           overwrite = TRUE)
rm(tmp)


print("D")

# Create heatmaps
for(featureType in names(variable_importance_df)){
  for (select_model in names(variable_importance_df[[featureType]])){
    print(paste(disease, featureType, select_model, data_balance_method, sep = " - "))
    tmp1 <- variable_importance_df[[featureType]][[select_model]]
    tmp1 <- tmp1[!rowSums(is.na(tmp1)) >= length(colnames(tmp1)[grep("^Fold", colnames(tmp1))]), ]
    if(nrow(tmp1) != 0){
      tmp1 <- column_to_rownames(tmp1, "Features")
         pheatmap(mat = as.matrix(tmp1), 
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               main = paste(disease, featureType, select_model, data_balance_method, sep = " - "),
               filename = paste0("Analysis/STRING/DrugCombs_v5/", disease, 
                                 "/featureImportance/varImp_models_",  
                                 disease, "_", featureType, "_", select_model, "_", data_balance_method, ".pdf"),
              width = 15,
              height = nrow(tmp1) * 0.5)
    }
    
  }
}