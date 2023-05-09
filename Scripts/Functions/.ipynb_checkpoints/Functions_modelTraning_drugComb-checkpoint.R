# Function to train model for drug combinations
# based on different types of feature matrix



# Load libraries
library(caret)
library(MLmetrics)
library(DMwR)
source("Scripts/Functions/Functions_parallelprocesses.R")


# MLmetrics::AUC (modified function)
# AUC calcution returns NA due to large value of n_pos * n_neg which cannot be handled by R
AUC <- function(y_pred, y_true){
            rank <- rank(y_pred)
            n_pos <- as.numeric(sum(y_true == 1))
            n_neg <- as.numeric(sum(y_true == 0))
            auc <- (sum(rank[y_true == 1]) - n_pos * (n_pos + 1)/2)/(n_pos * n_neg)
            return(auc)
}


# Parameters to train Random Forest model
train_rf_model <- function(x, y){
  caret::train(x = x,
               y = y,
               method = "rf",
               metric = "F",
               ntree = 1000,
               allowParallel = TRUE,
               tuneGrid = expand.grid(mtry = c(1:10)),
               trControl = trainControl(method = "repeatedcv",
                                        number = 5,
                                        repeats = 2,
                                        summaryFunction = prSummary,
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}


# Parameters to train glmnet model
train_glmnet_model <- function(x, y){
  caret::train(x = x,
               y = y,
               method = "glmnet",
               metric = "F",
               allowParallel = TRUE,
               tuneGrid = expand.grid(alpha = seq(0.1, 1, 0.1),
                                      lambda = 10^seq(-6, 3, by = 0.5)),
               trControl = trainControl(method = "repeatedcv",
                                        number = 5,
                                        repeats = 2,
                                        summaryFunction = prSummary,
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}


# Parameters to train SVM model
train_svmRadial_model <- function(x, y){
  caret::train(x = x, 
               y = y,
               method = "svmRadial",
               metric = "F",
               allowParallel = TRUE,
               tuneGrid = expand.grid(sigma = 10^(-4:2),
                                      C = 10^(-3:3)),
               trControl = trainControl(method = "repeatedcv",
                                        number = 5,
                                        repeats = 2,
                                        summaryFunction = prSummary, 
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}


# Parameters to train Naive Bayes model
train_nb_model <- function(x, y){
  caret::train(x = x,
               y = y,
               method = "naive_bayes",
               metric = "F",
               tuneGrid = expand.grid(laplace = seq(0, 10, 0.1), 
                                      adjust = seq(0, 1, 0.1),
                                      usekernel = TRUE),
               trControl = trainControl(method = "repeatedcv",
                                        number = 5,
                                        repeats = 2,
                                        summaryFunction = prSummary,
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}



# Parameters to train k-Nearest Neighbors model
train_knn_model <- function(x, y){
  caret::train(x = x,
               y = y,
               method = "kknn",
               metric = "F",
               tuneGrid = expand.grid(kmax = seq(0, 25, 5),
                                     distance = seq(1, 10, 1),
                                     kernel = c("optimal")), # "triangular", "rectangular", "epanechnikov", 
               trControl = trainControl(method = "repeatedcv",
                                        number = 5,
                                        repeats = 2,
                                        summaryFunction = prSummary,
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}


# Function to train all models
func_repeated_train <- function(feature_matrix, 
                                train_test_split, 
                                data_balance_method = "none", 
                                allow_parallel = TRUE, 
                                nproc = NULL){
  
  require(foreach)
  require(doParallel)
  
  
  if(!data_balance_method %in% c("SMOTE", "downSample", "upSample", "none")){
    stop("data_balance_method must be SMOTE, downSample, upSample or none")
  }
  
  if(allow_parallel){
    if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
    cl <- makeCluster(nproc)
    registerDoParallel(cl) 
  } 
  
  
  rf_result_table <- glmnet_result_table <- svmRadial_result_table <- nb_result_table <- knn_result_table <- data.frame()
  modelling_results <- list()
  
  
  # Prepare the input data
  data <- as.data.frame(t(feature_matrix[,-1]))
  colnames(data) <- feature_matrix[,1]
  data$Class <- factor(substr(row.names(data), 1, 3)) 
  
  
  for(i in names(train_test_split)){
    
    cat(paste0("-- ", i, "\n"))
    
    
    # Balance the training data (if used)
    switch(data_balance_method,
           "SMOTE" = {
             trainData <- data[row.names(data) %in% train_test_split[[i]]$train$Name,]
             trainData_smote <- SMOTE(Class ~ ., trainData, perc.over = 100, perc.under = 200)
             trainData <- trainData_smote[ , !(colnames(trainData_smote) %in% c("Class")), drop = FALSE]
             trainClass <- as.factor(trainData_smote[, c("Class")])
           },
           
           "downSample" = {
             trainData <- data[row.names(data) %in% train_test_split[[i]]$train$Name,]
             trainData_downSample <- downSample(x = trainData[ , !(colnames(trainData) %in% c("Class")), drop = FALSE],
                                                y = trainData$Class)
             
             trainData <- trainData_downSample[ , !(colnames(trainData_downSample) %in% c("Class")), drop = FALSE]
             trainClass <- as.factor(trainData_downSample[, c("Class")])
           },
           
           "upSample" = {
             trainData <- data[row.names(data) %in% train_test_split[[i]]$train$Name,]
             trainData_upSample <- upSample(x = trainData[ , !(colnames(trainData) %in% c("Class")), drop = FALSE],
                                            y = trainData$Class)
             
             trainData <- trainData_upSample[ , !(colnames(trainData_upSample) %in% c("Class")), drop = FALSE]
             trainClass <- as.factor(trainData_upSample[, c("Class")])
           },
           
           "none" = {
             trainData <- data[row.names(data) %in% train_test_split[[i]]$train$Name , !(colnames(data) %in% c("Class")), drop = FALSE]
             trainClass <- as.factor(data[row.names(data) %in% train_test_split[[i]]$train$Name, c("Class")])
           })
    
    testData <- data[row.names(data) %in% train_test_split[[i]]$test$Name, !(colnames(data) %in% c("Class")), drop = FALSE]
    testClass  <-  as.factor(data[row.names(data) %in% train_test_split[[i]]$test$Name, c("Class")])
    
    trainClass_count <- rbind(table(trainClass))
    colnames(trainClass_count) <- paste0("Train_", colnames(trainClass_count))
    testClass_count <- rbind(table(testClass))
    colnames(testClass_count) <- paste0("Test_", colnames(testClass_count))
      
    # Preprocess the data
    preProcess <- preProcess(trainData, method = c("zv", "center", "scale"))
    trainData <- predict(object = preProcess, newdata = trainData)
    testData <- predict(object = preProcess, newdata = testData)
    
    
      
      
    # Train Random Forest model
    cat("--- Random Forest\n")
    rf_model <- train_rf_model(x = trainData, y = trainClass)
    rf_predictions <- predict(object = rf_model, newdata = testData, type = "raw")
    rf_predictions_prob <- predict(object = rf_model, newdata = testData, type = "prob") 
    rf_confusionMatrix_train <- confusionMatrix(table(rf_model$pred$pred, rf_model$pred$obs), positive = "Eff")
    rf_confusionMatrix_test <- confusionMatrix(table(rf_predictions, testClass), positive = "Eff")
    rf_rocauc_train = AUC(y_pred = na.exclude(rf_model$pred)$Eff, y_true = ifelse(na.exclude(rf_model$pred)$obs == "Eff", 1, 0))
    rf_rocauc_test = AUC(y_pred = rf_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    rf_prauc_train = PRAUC(y_pred = na.exclude(rf_model$pred)$Eff, y_true = ifelse(na.exclude(rf_model$pred)$obs == "Eff", 1, 0))
    rf_prauc_test = PRAUC(y_pred = rf_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    rf_result_table <- rbind(rf_result_table,
                             data.frame(
                               Fold = i,
                               trainClass_count,
                               testClass_count,
                               BestTune_mtry = rf_model$bestTune$mtry,
                               PositiveClass = rf_confusionMatrix_test$positive,
                               ROCAUC_train = as.numeric(rf_rocauc_train),
                               ROCAUC_test = as.numeric(rf_rocauc_test),
                               PRAUC_train = as.numeric(rf_prauc_train),
                               PRAUC_test = as.numeric(rf_prauc_test),
                               Accuracy_train =  unname (rf_confusionMatrix_train$overall["Accuracy"]),
                               Accuracy_test =  unname (rf_confusionMatrix_test$overall["Accuracy"]),
                               Kappa_train =  unname (rf_confusionMatrix_train$overall["Kappa"]),
                               Kappa_test =  unname (rf_confusionMatrix_test$overall["Kappa"]),
                               Sensitivity_train =  unname (rf_confusionMatrix_train$byClass["Sensitivity"]),
                               Sensitivity_test =  unname (rf_confusionMatrix_test$byClass["Sensitivity"]),
                               Specificity_train =  unname (rf_confusionMatrix_train$byClass["Specificity"]),
                               Specificity_test =  unname (rf_confusionMatrix_test$byClass["Specificity"]),
                               Precision_train =  unname (rf_confusionMatrix_train$byClass["Precision"]),
                               Precision_test =  unname (rf_confusionMatrix_test$byClass["Precision"]),
                               Recall_train =  unname (rf_confusionMatrix_train$byClass["Recall"]),
                               Recall_test =  unname (rf_confusionMatrix_test$byClass["Recall"]),
                               F1_train =  unname (rf_confusionMatrix_train$byClass["F1"]),
                               F1_test =  unname (rf_confusionMatrix_test$byClass["F1"]),
                               BalancedAccuracy_train =  unname (rf_confusionMatrix_train$byClass["Balanced Accuracy"]),
                               BalancedAccuracy_test =  unname (unname (rf_confusionMatrix_test$byClass["Balanced Accuracy"]))))
    
    
    
    
    # Train glmnet model
    cat("--- glmnet\n")
    if(ncol(trainData) > 1){
      glmnet_model <- train_glmnet_model(x = trainData, y = trainClass)
      glmnet_predictions <- predict(object = glmnet_model, newdata = testData, type = "raw")
      glmnet_predictions_prob <- predict(object = glmnet_model, newdata = testData, type = "prob") 
      glmnet_confusionMatrix_train <- confusionMatrix(table(glmnet_model$pred$pred, glmnet_model$pred$obs), positive = "Eff")
      glmnet_confusionMatrix_test <- confusionMatrix(table(glmnet_predictions, testClass), positive = "Eff")
      glmnet_rocauc_train = AUC(y_pred = na.exclude(glmnet_model$pred)$Eff, y_true = ifelse(na.exclude(glmnet_model$pred)$obs == "Eff", 1, 0))
      glmnet_rocauc_test = AUC(y_pred = glmnet_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
      glmnet_prauc_train = PRAUC(y_pred = na.exclude(glmnet_model$pred)$Eff, y_true = ifelse(na.exclude(glmnet_model$pred)$obs == "Eff", 1, 0))
      glmnet_prauc_test = PRAUC(y_pred = glmnet_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
      glmnet_result_table <- rbind(glmnet_result_table,
                                   data.frame(
                                     Fold = i,
                                     trainClass_count,
                                     testClass_count,
                                     BestTune_alpha = glmnet_model$bestTune$alpha,
                                     BestTune_lambda = glmnet_model$bestTune$lambda,
                                     PositiveClass = glmnet_confusionMatrix_test$positive,
                                     ROCAUC_train = as.numeric(glmnet_rocauc_train),
                                     ROCAUC_test = as.numeric(glmnet_rocauc_test),
                                     PRAUC_train = as.numeric(glmnet_prauc_train),
                                     PRAUC_test = as.numeric(glmnet_prauc_test),
                                     Accuracy_train =  unname (glmnet_confusionMatrix_train$overall["Accuracy"]),
                                     Accuracy_test =  unname (glmnet_confusionMatrix_test$overall["Accuracy"]),
                                     Kappa_train =  unname (glmnet_confusionMatrix_train$overall["Kappa"]),
                                     Kappa_test =  unname (glmnet_confusionMatrix_test$overall["Kappa"]),
                                     Sensitivity_train =  unname (glmnet_confusionMatrix_train$byClass["Sensitivity"]),
                                     Sensitivity_test =  unname (glmnet_confusionMatrix_test$byClass["Sensitivity"]),
                                     Specificity_train =  unname (glmnet_confusionMatrix_train$byClass["Specificity"]),
                                     Specificity_test =  unname (glmnet_confusionMatrix_test$byClass["Specificity"]),
                                     Precision_train =  unname (glmnet_confusionMatrix_train$byClass["Precision"]),
                                     Precision_test =  unname (glmnet_confusionMatrix_test$byClass["Precision"]),
                                     Recall_train =  unname (glmnet_confusionMatrix_train$byClass["Recall"]),
                                     Recall_test =  unname (glmnet_confusionMatrix_test$byClass["Recall"]),
                                     F1_train =  unname (glmnet_confusionMatrix_train$byClass["F1"]),
                                     F1_test =  unname (glmnet_confusionMatrix_test$byClass["F1"]),
                                     BalancedAccuracy_train =  unname (glmnet_confusionMatrix_train$byClass["Balanced Accuracy"]),
                                     BalancedAccuracy_test =  unname (unname (glmnet_confusionMatrix_test$byClass["Balanced Accuracy"]))))
    }else{
      # Refer: https://stackoverflow.com/a/59414707/9450714
      warning("Training data contains only one feature column. Will return NA for glmnet")
      glmnet_model <- NA
      glmnet_predictions <- NA
      glmnet_predictions_prob <- NA
      glmnet_confusionMatrix_test <- NA
      glmnet_result_table <- NA
    }
    
    
    
    
    # Train Support Vector Machines with Radial Basis Function Kernel model
    cat("--- Support Vector Machine\n")
    svmRadial_model <- train_svmRadial_model(x = trainData, y = trainClass)
    svmRadial_predictions <- predict(object = svmRadial_model, newdata = testData, type = "raw")
    svmRadial_predictions_prob <- predict(object = svmRadial_model, newdata = testData, type = "prob") 
    svmRadial_confusionMatrix_train <- confusionMatrix(table(svmRadial_model$pred$pred, svmRadial_model$pred$obs), positive = "Eff")
    svmRadial_confusionMatrix_test <- confusionMatrix(table(svmRadial_predictions, testClass), positive = "Eff")
    svmRadial_rocauc_train = AUC(y_pred = na.exclude(svmRadial_model$pred)$Eff, y_true = ifelse(na.exclude(svmRadial_model$pred)$obs == "Eff", 1, 0))
    svmRadial_rocauc_test = AUC(y_pred = svmRadial_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    svmRadial_prauc_train = PRAUC(y_pred = na.exclude(svmRadial_model$pred)$Eff, y_true = ifelse(na.exclude(svmRadial_model$pred)$obs == "Eff", 1, 0))
    svmRadial_prauc_test = PRAUC(y_pred = svmRadial_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    svmRadial_result_table <- rbind(svmRadial_result_table,
                                    data.frame(
                                      Fold = i,
                                      trainClass_count,
                                      testClass_count,
                                      BestTune_sigma = svmRadial_model$bestTune$sigma,
                                      BestTune_C = svmRadial_model$bestTune$C,
                                      PositiveClass = svmRadial_confusionMatrix_test$positive,
                                      ROCAUC_train = as.numeric(svmRadial_rocauc_train),
                                      ROCAUC_test = as.numeric(svmRadial_rocauc_test),
                                      PRAUC_train = as.numeric(svmRadial_prauc_train),
                                      PRAUC_test = as.numeric(svmRadial_prauc_test),
                                      Accuracy_train =  unname (svmRadial_confusionMatrix_train$overall["Accuracy"]),
                                      Accuracy_test =  unname (svmRadial_confusionMatrix_test$overall["Accuracy"]),
                                      Kappa_train =  unname (svmRadial_confusionMatrix_train$overall["Kappa"]),
                                      Kappa_test =  unname (svmRadial_confusionMatrix_test$overall["Kappa"]),
                                      Sensitivity_train =  unname (svmRadial_confusionMatrix_train$byClass["Sensitivity"]),
                                      Sensitivity_test =  unname (svmRadial_confusionMatrix_test$byClass["Sensitivity"]),
                                      Specificity_train =  unname (svmRadial_confusionMatrix_train$byClass["Specificity"]),
                                      Specificity_test =  unname (svmRadial_confusionMatrix_test$byClass["Specificity"]),
                                      Precision_train =  unname (svmRadial_confusionMatrix_train$byClass["Precision"]),
                                      Precision_test =  unname (svmRadial_confusionMatrix_test$byClass["Precision"]),
                                      Recall_train =  unname (svmRadial_confusionMatrix_train$byClass["Recall"]),
                                      Recall_test =  unname (svmRadial_confusionMatrix_test$byClass["Recall"]),
                                      F1_train =  unname (svmRadial_confusionMatrix_train$byClass["F1"]),
                                      F1_test =  unname (svmRadial_confusionMatrix_test$byClass["F1"]),
                                      BalancedAccuracy_train =  unname (svmRadial_confusionMatrix_train$byClass["Balanced Accuracy"]),
                                      BalancedAccuracy_test =  unname (unname (svmRadial_confusionMatrix_test$byClass["Balanced Accuracy"]))))
    
    
    
    
    # Train Naive Bayes model
    cat("--- Naive Bayes\n")
    nb_model <- train_nb_model(x = trainData, y = trainClass)
    nb_predictions <- predict(object = nb_model, newdata = testData, type = "raw")
    nb_predictions_prob <- predict(object = nb_model, newdata = testData, type = "prob") 
    nb_confusionMatrix_train <- confusionMatrix(table(nb_model$pred$pred, nb_model$pred$obs), positive = "Eff")
    nb_confusionMatrix_test <- confusionMatrix(table(nb_predictions, testClass), positive = "Eff")
    nb_rocauc_train = AUC(y_pred = na.exclude(nb_model$pred)$Eff, y_true = ifelse(na.exclude(nb_model$pred)$obs == "Eff", 1, 0))
    nb_rocauc_test = AUC(y_pred = nb_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    nb_prauc_train = PRAUC(y_pred = na.exclude(nb_model$pred)$Eff, y_true = ifelse(na.exclude(nb_model$pred)$obs == "Eff", 1, 0))
    nb_prauc_test = PRAUC(y_pred = nb_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    nb_result_table <- rbind(nb_result_table,
                             data.frame(
                               Fold = i,
                               trainClass_count,
                               testClass_count,
                               BestTune_laplace = nb_model$bestTune$laplace,
                               BestTune_adjust = nb_model$bestTune$adjust,
                               BestTune_usekernel = nb_model$bestTune$usekernel,
                               PositiveClass = nb_confusionMatrix_test$positive,
                               ROCAUC_train = as.numeric(nb_rocauc_train),
                               ROCAUC_test = as.numeric(nb_rocauc_test),
                               PRAUC_train = as.numeric(nb_prauc_train),
                               PRAUC_test = as.numeric(nb_prauc_test),
                               Accuracy_train =  unname (nb_confusionMatrix_train$overall["Accuracy"]),
                               Accuracy_test =  unname (nb_confusionMatrix_test$overall["Accuracy"]),
                               Kappa_train =  unname (nb_confusionMatrix_train$overall["Kappa"]),
                               Kappa_test =  unname (nb_confusionMatrix_test$overall["Kappa"]),
                               Sensitivity_train =  unname (nb_confusionMatrix_train$byClass["Sensitivity"]),
                               Sensitivity_test =  unname (nb_confusionMatrix_test$byClass["Sensitivity"]),
                               Specificity_train =  unname (nb_confusionMatrix_train$byClass["Specificity"]),
                               Specificity_test =  unname (nb_confusionMatrix_test$byClass["Specificity"]),
                               Precision_train =  unname (nb_confusionMatrix_train$byClass["Precision"]),
                               Precision_test =  unname (nb_confusionMatrix_test$byClass["Precision"]),
                               Recall_train =  unname (nb_confusionMatrix_train$byClass["Recall"]),
                               Recall_test =  unname (nb_confusionMatrix_test$byClass["Recall"]),
                               F1_train =  unname (nb_confusionMatrix_train$byClass["F1"]),
                               F1_test =  unname (nb_confusionMatrix_test$byClass["F1"]),
                               BalancedAccuracy_train =  unname (nb_confusionMatrix_train$byClass["Balanced Accuracy"]),
                               BalancedAccuracy_test =  unname (unname (nb_confusionMatrix_test$byClass["Balanced Accuracy"]))))
    
    
    
    
    # Train k-Nearest Neighbors model
    cat("--- k-Nearest Neighbors\n")
    knn_model <- train_knn_model(x = trainData, y = trainClass)
    knn_predictions <- predict(object = knn_model, newdata = testData, type = "raw")
    knn_predictions_prob <- predict(object = knn_model, newdata = testData, type = "prob") 
    knn_confusionMatrix_train <- confusionMatrix(table(knn_model$pred$pred, knn_model$pred$obs), positive = "Eff")
    knn_confusionMatrix_test <- confusionMatrix(table(knn_predictions, testClass), positive = "Eff")
    knn_rocauc_train = AUC(y_pred = na.exclude(knn_model$pred)$Eff, y_true = ifelse(na.exclude(knn_model$pred)$obs == "Eff", 1, 0))
    knn_rocauc_test = AUC(y_pred = knn_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    knn_prauc_train = PRAUC(y_pred = na.exclude(knn_model$pred)$Eff, y_true = ifelse(na.exclude(knn_model$pred)$obs == "Eff", 1, 0))
    knn_prauc_test = PRAUC(y_pred = knn_predictions_prob$Eff, y_true = ifelse(testClass == "Eff", 1, 0))
    knn_result_table <- rbind(knn_result_table,
                             data.frame(
                               Fold = i,
                               trainClass_count,
                               testClass_count,
                               BestTune_kmax = knn_model$bestTune$kmax,
                               BestTune_distance = knn_model$bestTune$distance,
                               BestTune_kernel = knn_model$bestTune$kernel,
                               PositiveClass = knn_confusionMatrix_test$positive,
                               ROCAUC_train = as.numeric(knn_rocauc_train),
                               ROCAUC_test = as.numeric(knn_rocauc_test),
                               PRAUC_train = as.numeric(knn_prauc_train),
                               PRAUC_test = as.numeric(knn_prauc_test),
                               Accuracy_train =  unname (knn_confusionMatrix_train$overall["Accuracy"]),
                               Accuracy_test =  unname (knn_confusionMatrix_test$overall["Accuracy"]),
                               Kappa_train =  unname (knn_confusionMatrix_train$overall["Kappa"]),
                               Kappa_test =  unname (knn_confusionMatrix_test$overall["Kappa"]),
                               Sensitivity_train =  unname (knn_confusionMatrix_train$byClass["Sensitivity"]),
                               Sensitivity_test =  unname (knn_confusionMatrix_test$byClass["Sensitivity"]),
                               Specificity_train =  unname (knn_confusionMatrix_train$byClass["Specificity"]),
                               Specificity_test =  unname (knn_confusionMatrix_test$byClass["Specificity"]),
                               Precision_train =  unname (knn_confusionMatrix_train$byClass["Precision"]),
                               Precision_test =  unname (knn_confusionMatrix_test$byClass["Precision"]),
                               Recall_train =  unname (knn_confusionMatrix_train$byClass["Recall"]),
                               Recall_test =  unname (knn_confusionMatrix_test$byClass["Recall"]),
                               F1_train =  unname (knn_confusionMatrix_train$byClass["F1"]),
                               F1_test =  unname (knn_confusionMatrix_test$byClass["F1"]),
                               BalancedAccuracy_train =  unname (knn_confusionMatrix_train$byClass["Balanced Accuracy"]),
                               BalancedAccuracy_test =  unname (unname (knn_confusionMatrix_test$byClass["Balanced Accuracy"]))))
    
    
    # Compile results
    tmp <- list(rf = list(model = rf_model, 
                          test_predictions = rf_predictions, 
                          test_probabilities = rf_predictions_prob, 
                          test_confusionMatrix = rf_confusionMatrix_test),
                glmnet = list(model = glmnet_model, 
                              test_predictions = glmnet_predictions, 
                              test_probabilities = glmnet_predictions_prob, 
                              test_confusionMatrix = glmnet_confusionMatrix_test),
                svmRadial = list(model = svmRadial_model, 
                                 test_predictions = svmRadial_predictions, 
                                 test_probabilities = svmRadial_predictions_prob, 
                                 test_confusionMatrix = svmRadial_confusionMatrix_test),
                 nb = list(model = nb_model, 
                           test_predictions = nb_predictions, 
                           test_probabilities = nb_predictions_prob, 
                           test_confusionMatrix = nb_confusionMatrix_test),
                 knn = list(model = knn_model, 
                            test_predictions = knn_predictions, 
                            test_probabilities = knn_predictions_prob, 
                            test_confusionMatrix = knn_confusionMatrix_test))
    
    modelling_results[[i]] <- tmp
    
  }

  result_summary_tables <- list(rf = rf_result_table,
                                glmnet = glmnet_result_table,
                                svmRadial = svmRadial_result_table,
                                nb = nb_result_table,
                                knn = knn_result_table)
  
  
  return(list(result_summary_tables = result_summary_tables, 
              modelling_results = modelling_results))
  
  
  if(allow_parallel){
    stopCluster(cl)
    unregister_dopar()
    Sys.sleep(300)
  } 
  
}