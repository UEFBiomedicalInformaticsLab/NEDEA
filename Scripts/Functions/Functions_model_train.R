# Function to train model for drug combinations



# Load libraries
library(caret)
library(MLmetrics)






# MLmetrics::AUC (modified function)
# AUC calcution returns NA due to large value of n_pos * n_neg which cannot be handled by R
AUC <- function(y_pred, y_true){
  rank <- rank(y_pred)
  n_pos <- as.numeric(sum(y_true == 1))
  n_neg <- as.numeric(sum(y_true == 0))
  auc <- (sum(rank[y_true == 1]) - n_pos * (n_pos + 1)/2)/(n_pos * n_neg)
  return(auc)
}







# Parameters to train glmnet model
train_glmnet_model <- function(x, y){
  caret::train(x = x,
               y = y,
               method = "glmnet",
               metric = "F",
               allowParallel = TRUE,
               tuneGrid = expand.grid(alpha = seq(0, 1, 0.1),
                                      lambda = 10^seq(-6, 3, by = 0.5)),
               trControl = trainControl(method = "cv",
                                        number = 5,
                                        summaryFunction = prSummary,
                                        classProbs = TRUE,
                                        savePredictions = TRUE))
}


















###########################################################




data <- fgsea_result_select 
train_test_split <- train_test_split
data_balance_method <- data_balance_method



require(foreach)
require(doParallel)


# Prepare the input data



for(fold in names(train_test_split)){}


cat(paste0("-- ", fold, "\n"))


# Balance the training data (if used)
switch(data_balance_method,
       "none" = {
         trainData <- data[row.names(data) %in% train_test_split[[fold]]$train$Name , !(colnames(data) %in% c("Class")), drop = FALSE]
         trainClass <- as.factor(data[row.names(data) %in% train_test_split[[i]]$train$Name, c("Class")])
       })
















