set.seed(5081)



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
# library(pROC)
# library(caret)
# library(foreach)
# library(doParallel)
library(yardstick)
source("Scripts/Functions/Functions_parallelprocesses.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


# # Get arguments
# option_list = list(
#   make_option(c("--disease"), type = "character", default = NULL, 
#               help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
#   make_option(c("--drug_target_type"), type = "character", default = "known", 
#               help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
#   make_option(c("--nproc"), type = "numeric", default = NULL, 
#               help = "Number of processes to use. Default: NULL", metavar = "numeric")
# )
# 
# opt_parser = OptionParser(option_list = option_list)
# opt = parse_args(opt_parser)
# 
# 
# 
# # Define global options for this script 
# disease <- opt$disease
# drug_target_type <- opt$drug_target_type
# nproc <- opt$nproc

disease <- "BreastCancer"
drug_target_type <- "known"
nproc <- 3

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




method1_prediction_mat <- data.frame("drugCombs" = row.names(fgsea_result_select))
method2_prediction_mat <- data.frame("drugCombs" = row.names(fgsea_result_select))
method3_prediction_mat <- data.frame("drugCombs" = row.names(fgsea_result_select))
method4_prediction_mat <- data.frame("drugCombs" = row.names(fgsea_result_select))

selected_threshold <- as.data.frame(matrix(NA, nrow = length(selected_features$feature), ncol = 2, dimnames = list(selected_features$feature, c("M2", "M3"))))

plot_method_2 <- list()
plot_method_3 <- list()

metrics_df <- data.frame()

for(feature_name in selected_features$feature){
  
  
  
  #####
  
  
  
  # Prepare the data for fitting
  fit_data <- fgsea_result_select[, c(feature_name, "category")]
  colnames(fit_data)[1] <- "feature"
  tmp1 <- c("Eff" = 1, "Adv" = 0)
  fit_data$category <- unname(tmp1[fit_data$category])
  
  
  # Fit the model
  model <- glm(category ~ feature , data = fit_data, family = "binomial")
  fit_data$predicted_probs <- predict(model, newdata = fit_data, type = "response")
  
  
  
  #####
  
  
  
  # METHOD 1: Average of predicted probabilities (by group)
  # Save the predictions into the final dataframe
  tmp1 <- fit_data[, c("predicted_probs"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  method1_prediction_mat <- full_join(method1_prediction_mat, tmp1, by = "drugCombs")
  rm(tmp1)
  
  
  
  #####
  
  
  
  # METHOD 2: Use ROC directly on NES to find threshold  
  
  roc_obj <- pROC::roc(fit_data$category, fit_data$feature)
  best_feature_threshold_2 <- pROC::coords(roc_obj, x = "best")
  selected_threshold[feature_name, "M2"] <- best_feature_threshold_2$threshold
  tmp1 <- fit_data
  if(grepl("^\\[DISEASE\\]", feature_name)){
    tmp1$method2_predicted_class <- ifelse(fit_data$feature >= best_feature_threshold_2$threshold, 1, 0)
  }
  if(grepl("^\\[ADR\\]", feature_name)){
    tmp1$method2_predicted_class <- ifelse(fit_data$feature >= best_feature_threshold_2$threshold, 0, 1)
  }
  
  
  
  # plot_data <- tmp1
  # plot_data$feature_name <- feature_name
  # 
  # plot_data$category <- as.factor(plot_data$category)
  # plot_data$method2_predicted_class <- as.factor(plot_data$method2_predicted_class)
  # 
  # plot_method_2[[feature_name]] <- ggplot(plot_data, aes(x = feature_name, y = feature, group = category)) +
  #   geom_boxplot(outliers = FALSE) +
  #   geom_jitter(aes(shape = category, color = method2_predicted_class)) +
  #   geom_hline(yintercept = best_feature_threshold_2$threshold, 
  #              linetype="dashed",
  #              color = "red") +
  #   scale_x_discrete(labels = function(x) scales::label_wrap(30)(x))
  
  
  
  tmp1 <- tmp1[, c("method2_predicted_class"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  method2_prediction_mat <- full_join(method2_prediction_mat, tmp1, by = "drugCombs")
  rm(tmp1)
  
  
  
  #####
  
  
  
  # METHOD 3: Use ROC on predicted probabilities to find threshold
  
  roc_obj <- pROC::roc(fit_data$category, fit_data$predicted_probs)
  best_predictedProb_threshold <- pROC::coords(roc_obj, x = "best")
  best_feature_threshold_3 <- (log(best_predictedProb_threshold$threshold / (1 - best_predictedProb_threshold$threshold)) - as.numeric(coef(model)[1])) / as.numeric(coef(model)[2])
  selected_threshold[feature_name, "M3"] <- best_feature_threshold_3
  
  tmp1 <- fit_data
  if(grepl("^\\[DISEASE\\]", feature_name)){
    tmp1$method3_predicted_class <- ifelse(fit_data$feature >= best_feature_threshold_3, 1, 0)
  }
  if(grepl("^\\[ADR\\]", feature_name)){
    tmp1$method3_predicted_class <- ifelse(fit_data$feature >= best_feature_threshold_3, 0, 1)
  }  
  
  
  
  # plot_data <- tmp1
  # plot_data$feature_name <- feature_name
  # 
  # plot_data$category <- as.factor(plot_data$category)
  # plot_data$method3_predicted_class <- as.factor(plot_data$method3_predicted_class)
  # 
  # plot_method_3[[feature_name]] <- ggplot(plot_data, aes(x = feature_name, y = feature, group = category)) +
  #   geom_boxplot(outliers = FALSE) +
  #   geom_jitter(aes(shape = category, color = method3_predicted_class)) +
  #   geom_hline(yintercept = best_feature_threshold_3, 
  #              linetype="dashed",
  #              color = "red") +
  #   scale_x_discrete(labels = function(x) scales::label_wrap(30)(x))
  
  
  tmp1 <- tmp1[, c("method3_predicted_class"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  method3_prediction_mat <- full_join(method3_prediction_mat, tmp1, by = "drugCombs")
  rm(tmp1)
  
  
  
  #####
  
  
  # METHOD 4: Average of predicted probabilities
  # Save the predictions into the final dataframe
  tmp1 <- fit_data[, c("predicted_probs"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  method4_prediction_mat <- full_join(method4_prediction_mat, tmp1, by = "drugCombs")
  rm(tmp1)
  
  
}






# METHOD 1: Average of predicted probabilities

method1_prediction_mat <- column_to_rownames(method1_prediction_mat, "drugCombs")
method1_prediction_mat <- method1_prediction_mat %>% 
  mutate(efficacy_prob = rowMeans(select(., starts_with("[DISEASE]"))), 
         safety_prob = rowMeans(select(., starts_with("[ADR]")))) %>%
  mutate(final_prob = rowMeans(select(., c("efficacy_prob", "safety_prob"))))
method1_prediction_mat$final_class <- ifelse(method1_prediction_mat$final_prob > 0.5, "Eff", "Adv")
method1_prediction_mat$final_class <- factor(method1_prediction_mat$final_class, levels = c("Eff", "Adv"))


method1_prediction_mat$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(method1_prediction_mat), drugCombs_cat$comb_name)]
method1_prediction_mat$class_EffAdv <- factor(method1_prediction_mat$class_EffAdv, levels = c("Eff", "Adv"))

metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy, ppv, npv)
tmp2 <- metrics(data = method1_prediction_mat, truth = class_EffAdv, estimate = final_class)
tmp2$method <- "method1"
metrics_df <- rbind(metrics_df, tmp2)
rm(tmp2)







# METHOD 2: Use ROC directly on NES to find threshold  

method2_prediction_mat <- column_to_rownames(method2_prediction_mat, "drugCombs")
method2_prediction_mat <- method2_prediction_mat %>% 
  mutate(efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
         safety_score = rowMeans(select(., starts_with("[ADR]")))) 
method2_prediction_mat$final_score <- method2_prediction_mat$efficacy_score - method2_prediction_mat$safety_score


method2_prediction_mat$final_class <- ifelse(method2_prediction_mat$final_score >= 0, "Eff", "Adv")
method2_prediction_mat$final_class <- factor(method2_prediction_mat$final_class, levels = c("Eff", "Adv"))


method2_prediction_mat$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(method2_prediction_mat), drugCombs_cat$comb_name)]
method2_prediction_mat$class_EffAdv <- factor(method2_prediction_mat$class_EffAdv, levels = c("Eff", "Adv"))

metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy, ppv, npv)
tmp2 <- metrics(data = method2_prediction_mat, truth = class_EffAdv, estimate = final_class)
tmp2$method <- "method2"
metrics_df <- rbind(metrics_df, tmp2)
rm(tmp2)






# METHOD 3: Use ROC on predicted probabilities to find threshold
method3_prediction_mat <- column_to_rownames(method3_prediction_mat, "drugCombs")

method3_prediction_mat <- method3_prediction_mat %>% 
  mutate(efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
         safety_score = rowMeans(select(., starts_with("[ADR]")))) 
method3_prediction_mat$final_score <- method3_prediction_mat$efficacy_score - method3_prediction_mat$safety_score


method3_prediction_mat$final_class <- ifelse(method3_prediction_mat$final_score >= 0, "Eff", "Adv")
method3_prediction_mat$final_class <- factor(method3_prediction_mat$final_class, levels = c("Eff", "Adv"))


method3_prediction_mat$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(method3_prediction_mat), drugCombs_cat$comb_name)]
method3_prediction_mat$class_EffAdv <- factor(method3_prediction_mat$class_EffAdv, levels = c("Eff", "Adv"))

metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy, ppv, npv)
tmp2 <- metrics(data = method3_prediction_mat, truth = class_EffAdv, estimate = final_class)
tmp2$method <- "method3"
metrics_df <- rbind(metrics_df, tmp2)
rm(tmp2)






# METHOD 4: Average of predicted probabilities

method4_prediction_mat <- column_to_rownames(method4_prediction_mat, "drugCombs")
method4_prediction_mat$final_prob <- rowMeans(method4_prediction_mat)

method4_prediction_mat$final_class <- ifelse(method4_prediction_mat$final_prob > 0.5, "Eff", "Adv")
method4_prediction_mat$final_class <- factor(method4_prediction_mat$final_class, levels = c("Eff", "Adv"))


method4_prediction_mat$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(method4_prediction_mat), drugCombs_cat$comb_name)]
method4_prediction_mat$class_EffAdv <- factor(method4_prediction_mat$class_EffAdv, levels = c("Eff", "Adv"))

metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy, ppv, npv)
tmp2 <- metrics(data = method4_prediction_mat, truth = class_EffAdv, estimate = final_class)
tmp2$method <- "method4"
metrics_df <- rbind(metrics_df, tmp2)
rm(tmp2)


# ggpubr::ggarrange(plotlist = plot_method_2)
# ggpubr::ggarrange(plotlist = plot_method_3)


