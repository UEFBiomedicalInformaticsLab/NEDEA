set.seed(5081)



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(caret)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known", 
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type
nproc <- opt$nproc

cat(paste0("\n\nTraining model for: ", disease, "\n"))


# Read the FGSEA result
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
fgsea_result <- fgsea_result[["combinedEfficacySafety"]]


# Read the drug combination category
drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
drugCombs_cat <- drugCombs_cat[, c("comb_name", "class_EffAdv")]
drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), ]


# Read the list of selected features
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_known.csv"))


# Sub-set the FGSEA results for only the selected drug combinations and features
fgsea_result_select <- fgsea_result[row.names(fgsea_result) %in% selected_features$feature, 
                                    colnames(fgsea_result) %in% drugCombs_cat$comb_name]
fgsea_result_select <- as.data.frame(t(fgsea_result_select))
fgsea_result_select$category <-  drugCombs_cat$class_EffAdv[match(row.names(fgsea_result_select), drugCombs_cat$comb_name)]


# Prepare the data for training
trainData <- fgsea_result_select[, !colnames(fgsea_result_select) %in% c("category")]
trainClass <- fgsea_result_select$category


# # Tune the hyper-parameters 
# # Doing separately as caret does not allow tuning the number of trees in random forest
# tune_grid <- expand.grid(mtry = c(1:10), 
#                          ntree = seq(50, 1000, 50))
# 
# tuning_folds <- createFolds(y = trainClass, 
#                             k = 5, 
#                             list = TRUE, 
#                             returnTrain = TRUE)
# 
# cl <- makeCluster(nproc)
# registerDoParallel(cl) 
# 
# tune_result <- data.frame()
# for(i in 1:nrow(tune_grid)){
#   tmp1 <- train(x = trainData,
#                 y = trainClass,
#                 method = "rf",
#                 metric = "F",
#                 ntree = tune_grid[i, "ntree"],
#                 allowParallel = TRUE,
#                 tuneGrid = data.frame(mtry = tune_grid[i, "mtry"]),
#                 trControl = trainControl(index = tuning_folds, 
#                                          method = "cv",
#                                          number = 5,
#                                          summaryFunction = prSummary,
#                                          classProbs = TRUE,
#                                          savePredictions = TRUE))
#   tmp1 <- tmp1$results
#   tmp1$ntree <- tune_grid[i, "ntree"]
#   tune_result <- rbind(tune_result, tmp1)
# }
# rm(tmp1)
# 
# stopCluster(cl)
# unregister_dopar()
# 
# 
# 
# # Select the best tune 
# best_tune <- tune_result[tune_result$F == max(tune_result$F, na.rm = TRUE), ]
# if(nrow(best_tune) > 1){
#   best_tune <- best_tune[best_tune$FSD == min(best_tune$F, na.rm = TRUE), ]
# }
# 
# best_tune <- best_tune[, c("mtry", "ntree")]
# 
# cat(paste0("\n\nSelected mtry: ", best_tune$mtry, "\n"))
# cat(paste0("\n\nSelected ntree: ", best_tune$ntree, "\n"))


# Final model training
# model <- train(x = trainData,
#                y = trainClass,
#                method = "rf",
#                metric = "F",
#                ntree = best_tune$ntree,
#                allowParallel = TRUE,
#                tuneGrid = data.frame(mtry = best_tune$mtry),
#                trControl = trainControl(method = "none"))

# Using mtry = 2 since only one efficacy and one safety feature should be sufficient 

model <- train(x = trainData,
               y = trainClass,
               method = "rf",
               metric = "F",
               ntree = 500,
               allowParallel = TRUE,
               tuneGrid = data.frame(mtry = 2),
               trControl = trainControl(method = "none"))


# Save the models
if(!dir.exists("OutputFiles/Predictive_model/")){ dir.create("OutputFiles/Predictive_model/", recursive = TRUE) }
saveRDS(model, file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))



print(warnings())