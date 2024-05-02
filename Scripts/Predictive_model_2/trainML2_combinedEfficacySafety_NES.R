set.seed(5081)


# Script to train predictive model (type 2)



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


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known",
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type
# nproc <- opt$nproc


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


#####


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
selected_features <- selected_features[order(selected_features$feature), ]

# Sub-set the FGSEA results for only the selected drug combinations and features
fgsea_result_select <- fgsea_result[row.names(fgsea_result) %in% selected_features$feature, 
                                    colnames(fgsea_result) %in% drugCombs_cat$comb_name]
fgsea_result_select <- as.data.frame(t(fgsea_result_select))
fgsea_result_select$category <-  drugCombs_cat$class_EffAdv[match(row.names(fgsea_result_select), drugCombs_cat$comb_name)]


#####


# Create different folds of data to verify the stability
data_folds <- list()
for(i in seq(1:10)){
  set.seed(5081 + i*2)
  data_folds[[paste0("Repeat", i)]] <- caret::createFolds(y = fgsea_result_select$category,
                                                                    k = 3,
                                                                    list = TRUE,
                                                                    returnTrain = TRUE)
}

data_folds <- unlist(data_folds, recursive = FALSE)
names(data_folds) <- gsub("\\.", "_", names(data_folds))
set.seed(5081)


#####


metrics_list <- list()
metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy, ppv, npv)


for(fold in names(data_folds)){
  
  train_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select[data_folds[[fold]], ]))
  test_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select[-data_folds[[fold]], ]))
  
  
  for(feature_name in selected_features$feature){
    
    # Prepare the data for fitting
    fit_data <- fgsea_result_select[data_folds[[fold]], c(feature_name, "category")]
    colnames(fit_data)[1] <- "feature"
    tmp1 <- c("Eff" = 1, "Adv" = 0)
    fit_data$category <- unname(tmp1[fit_data$category])
    rm(tmp1)
    
    
    # Fit the model
    model <- glm(category ~ feature , data = fit_data, family = "binomial")
    
    
    # Predict on the training data
    fit_data$predicted_probability <- predict(model, newdata = fit_data, type = "response")
    tmp1 <- fit_data[, c("predicted_probability"), drop = FALSE]
    colnames(tmp1) <- feature_name
    tmp1 <- rownames_to_column(tmp1, "drugCombs")
    train_prediction_df <- full_join(train_prediction_df, tmp1, by = "drugCombs")
    rm(tmp1)
    
    
    # Predict on the test data
    predict_data <- fgsea_result_select[-data_folds[[fold]], c(feature_name, "category")]
    colnames(predict_data)[1] <- "feature"
    predict_data$predicted_probability <- predict(model, newdata = predict_data, type = "response")
    tmp1 <- predict_data[, c("predicted_probability"), drop = FALSE]
    colnames(tmp1) <- feature_name
    tmp1 <- rownames_to_column(tmp1, "drugCombs")
    test_prediction_df <- full_join(test_prediction_df, tmp1, by = "drugCombs")
    rm(tmp1)
    
  }
  

  # Calculate accuracy metrics for training 
  train_prediction_df <- column_to_rownames(train_prediction_df, "drugCombs")
  train_prediction_df$final_predicted_probability <- rowMeans(train_prediction_df)
  
  train_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(train_prediction_df), drugCombs_cat$comb_name)]
  train_prediction_df$class_EffAdv <- factor(train_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))
  
  roc_curve_data <- roc_curve(data = train_prediction_df, truth = class_EffAdv, final_predicted_probability)
  roc_curve_data$Youden_stat <- roc_curve_data$sensitivity + roc_curve_data$specificity - 1
  best_probability_threshold <- roc_curve_data[roc_curve_data$Youden_stat == max(roc_curve_data$Youden_stat), ".threshold", drop = TRUE]
  
  train_prediction_df$final_predicted_category <- ifelse(train_prediction_df$final_predicted_probability > best_probability_threshold, "Eff", "Adv")
  train_prediction_df$final_predicted_category <- factor(train_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))
  
  metrics_list[[fold]][["Train"]] <- rbind( metrics(data = train_prediction_df, truth = class_EffAdv, estimate = final_predicted_category),
                                            data.frame(".metric" = "probability_threshold", ".estimator" = NA, ".estimate" = best_probability_threshold) )
  
  
  # Calculate accuracy metrics for test 
  test_prediction_df <- column_to_rownames(test_prediction_df, "drugCombs")
  test_prediction_df$final_predicted_probability <- rowMeans(test_prediction_df)
  
  test_prediction_df$final_predicted_category <- ifelse(test_prediction_df$final_predicted_probability > best_probability_threshold, "Eff", "Adv")
  test_prediction_df$final_predicted_category <- factor(test_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))
  
  test_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(test_prediction_df), drugCombs_cat$comb_name)]
  test_prediction_df$class_EffAdv <- factor(test_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))
  
  metrics_list[[fold]][["Test"]] <- metrics(data = test_prediction_df, truth = class_EffAdv, estimate = final_predicted_category)
  
}


#####


# Plot the metrics
plot_data <- unlist(metrics_list, recursive = FALSE)

plot_data <- lapply(metrics_list, function(x){ bind_rows(x, .id = "group") })
plot_data <- bind_rows(plot_data, .id = "fold")

plot_data$.metric <- factor(plot_data$.metric, levels = c("probability_threshold", "accuracy", "f_meas", "sensitivity", "specificity", "recall", "precision", "bal_accuracy", "ppv", "npv"))
plot_data$group <- factor(plot_data$group, levels = c("Train", "Test"))



if(!dir.exists("OutputFiles/Predictive_model_2/")){ dir.create("OutputFiles/Predictive_model_2/") }

tiff(paste0("OutputFiles/Predictive_model_2/metrics_PredictiveModel2_", disease, "_", drug_target_type, ".tiff"), 
     width = 8,
     height = 5,
     units = "cm", compression = "lzw", res = 1200)


ggplot(plot_data, aes(x = .metric, y = .estimate, fill = group)) +
  geom_boxplot(width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3, 
               outlier.size = 0.5, 
               outlier.stroke = 0.1) +
  ylim(0, 1) +
  scale_fill_manual(values=c("Test" = "#E69F00", "Train" = "#56B4E9")) +
  scale_x_discrete(labels = function(x) scales::label_wrap(20)(x)) +
  labs(x = "Metric", 
       y = "Score",
       fill = "Data") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        # strip.background = element_rect(color = "black", linewidth = 0.25,),
        # strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
        text = element_text(size = 3), 
        # plot.title = element_text(size = 4, hjust = 0.5, face = "plain"),
        axis.text.x = element_text(size = 3, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()


#####



# Train model with complete data


allData_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select))
model_list <- list()


for(feature_name in selected_features$feature){
  
  # Prepare the data for fitting
  fit_data <- fgsea_result_select[, c(feature_name, "category")]
  colnames(fit_data)[1] <- "feature"
  tmp1 <- c("Eff" = 1, "Adv" = 0)
  fit_data$category <- unname(tmp1[fit_data$category])
  rm(tmp1)
  
  
  # Fit the model
  model <- glm(category ~ feature , data = fit_data, family = "binomial")
  model_list[[feature_name]] <- model
  
  # Predict on the training data
  fit_data$predicted_probability <- predict(model, newdata = fit_data, type = "response")
  tmp1 <- fit_data[, c("predicted_probability"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  allData_prediction_df <- full_join(allData_prediction_df, tmp1, by = "drugCombs")
  rm(tmp1)
  
}


allData_prediction_df <- column_to_rownames(allData_prediction_df, "drugCombs")
allData_prediction_df$final_predicted_probability <- rowMeans(allData_prediction_df)

allData_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(allData_prediction_df), drugCombs_cat$comb_name)]
allData_prediction_df$class_EffAdv <- factor(allData_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))

roc_curve_data <- roc_curve(data = allData_prediction_df, truth = class_EffAdv, final_predicted_probability)
roc_curve_data$Youden_stat <- roc_curve_data$sensitivity + roc_curve_data$specificity - 1
best_probability_threshold <- roc_curve_data[roc_curve_data$Youden_stat == max(roc_curve_data$Youden_stat), ".threshold", drop = TRUE]

allData_prediction_df$final_predicted_category <- ifelse(allData_prediction_df$final_predicted_probability > best_probability_threshold, "Eff", "Adv")
allData_prediction_df$final_predicted_category <- factor(allData_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))


#####


# Save the final model

final_model <- list("models" = model_list, "probability_threshold" = best_probability_threshold)

if(!dir.exists("OutputFiles/Predictive_model_2/")){ dir.create("OutputFiles/Predictive_model_2/") }
saveRDS(final_model, file = paste0("OutputFiles/Predictive_model_2/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


print(warnings())