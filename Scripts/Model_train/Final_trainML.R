# Script to train final model
# Note:
# (1) Uses the complete dataset to train a final model using the selected feature type





# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(openxlsx)
source("Scripts/Functions/Functions_modelTraning_drugComb.R")
source("Scripts/Functions/Functions_parallelprocesses.R")





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: SMOTE, downSample, upSample, or none. Default: none.", metavar = "character"),
  make_option(c("--model"), type = "character", default = NULL, 
            help = "The modelling technique to use. Possible values: glmnet, knn, nb, rf, svmRadial. Default: NULL", metavar = "character"),
  make_option(c("--feature_type"), type = "character", default = NULL, 
            help = "The feature type to use for modelling. Possible values: Disease2Gene, WithdrawalAdr2Gene, CombinedDisAdr2Gene, keggPath, SMPDbPath_DrugMet, SMPDbPath_DrugAction, miscGeneSet. Default: NULL", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
  
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$data_balance_method %in% c("SMOTE", "downSample", "upSample", "none")){
  stop("--data_balance_method must be SMOTE, downSample, upSample or none")
}

if(is.null(opt$model)){
  print_help(opt_parser)
  stop("--model argument needed", call.=FALSE)
}

if(!opt$model %in% c("glmnet", "knn", "nb", "rf", "svmRadial")){
  stop("--model must be glmnet, knn, nb, rf, svmRadial")
}

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}


# Define global options for this script 
disease <- opt$disease
data_balance_method <- opt$data_balance_method
model <- opt$model
feature_type <- opt$feature_type
nproc <- opt$nproc


# disease <- "LungCancer"
# data_balance_method <- "none"
# model <- "rf"
# feature_type <- "Disease2Gene"
# nproc <- 10




cat(paste0("\n\nTraining final model\n"))
cat(paste0("Disease: ", disease, "\n"))
cat(paste0("Data balance method: ", data_balance_method, "\n"))
cat(paste0("Model: ", model, "\n"))
cat(paste0("Feature type: ", feature_type, "\n"))





# Read the features
switch(feature_type,
       "Disease2Gene" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_Disease2Gene
       }, 
       "WithdrawalAdr2Gene" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_WithdrawalAdr2Gene
       }, 
       "CombinedDisAdr2Gene" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
       }, 
       "keggPath" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_keggPath
       }, 
       "SMPDbPath_DrugMet" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_SMPDbPath_DrugMet
       }, 
       "SMPDbPath_DrugAction" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_SMPDbPath_DrugAction
       }, 
       "miscGeneSet" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_miscGeneSet
       }
      )





# Prepare the input data
data <- as.data.frame(t(feature_matrix[,-1]))
colnames(data) <- feature_matrix[,1]
data$Class <- factor(substr(row.names(data), 1, 3)) 
  




# Balance the data (if used)
switch(data_balance_method,
       "SMOTE" = {
         trainData_smote <- SMOTE(Class ~ ., data, perc.over = 100, perc.under = 200)
         trainData <- trainData_smote[ , !(colnames(trainData_smote) %in% c("Class")), drop = FALSE]
         trainClass <- as.factor(trainData_smote[, c("Class")])
       },

       "downSample" = {
         trainData_downSample <- downSample(x = data[ , !(colnames(data) %in% c("Class")), drop = FALSE],
                                            y = data$Class)

         trainData <- trainData_downSample[ , !(colnames(trainData_downSample) %in% c("Class")), drop = FALSE]
         trainClass <- as.factor(trainData_downSample[, c("Class")])
       },

       "upSample" = {
         trainData_upSample <- upSample(x = data[ , !(colnames(data) %in% c("Class")), drop = FALSE],
                                        y = data$Class)

         trainData <- trainData_upSample[ , !(colnames(trainData_upSample) %in% c("Class")), drop = FALSE]
         trainClass <- as.factor(trainData_upSample[, c("Class")])
       },

       "none" = {
         trainData <- data[, !(colnames(data) %in% c("Class")), drop = FALSE]
         trainClass <- as.factor(data[, c("Class")])
       })





# Preprocess the data
preProcess <- preProcess(trainData, method = c("zv", "center", "scale"))
trainData <- predict(object = preProcess, newdata = trainData)
 



# Train model

if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 


switch(model,
        "glmnet" = { trained_model <- train_glmnet_model(x = trainData, y = trainClass) }, 
        "knn" = { trained_model <- train_knn_model(x = trainData, y = trainClass) }, 
        "nb" = { trained_model <- train_nb_model(x = trainData, y = trainClass) }, 
        "rf" = { trained_model <- train_rf_model(x = trainData, y = trainClass) },  
        "svmRadial" = { trained_model <- train_svmRadial_model(x = trainData, y = trainClass) })


stopCluster(cl)
unregister_dopar()
Sys.sleep(60)


# Calculate the accuracy metrices
confusionMatrix <- confusionMatrix(table(trained_model$pred$pred, trained_model$pred$obs), positive = "Eff")
rocauc <- AUC(y_pred = na.exclude(trained_model$pred)$Eff, y_true = ifelse(na.exclude(trained_model$pred)$obs == "Eff", 1, 0))
prauc <- PRAUC(y_pred = na.exclude(trained_model$pred)$Eff, y_true = ifelse(na.exclude(trained_model$pred)$obs == "Eff", 1, 0))
sensitivity <- unname(confusionMatrix$byClass["Sensitivity"])
specificity <- unname(confusionMatrix$byClass["Specificity"])
precision <- unname(confusionMatrix$byClass["Precision"])
recall <- unname(confusionMatrix$byClass["Recall"])
f1 <- unname(confusionMatrix$byClass["F1"])
balancedAccuracy <- unname (confusionMatrix$byClass["Balanced Accuracy"])



                            
                            
# Save as file
if(!dir.exists("OutputFiles/Final_model/")){
  dir.create("OutputFiles/Final_model/", recursive = TRUE)
}                          
saveRDS(trained_model, paste0("OutputFiles/Final_model/FinalModel_", disease, "_", model, "_", data_balance_method, "_", feature_type, ".rds"))          

                            
                            
print(warnings())