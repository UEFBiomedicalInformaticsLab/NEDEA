set.seed(5081)



# Use the trained model to predict the labels on the validation data 2



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(yardstick)
library(gridExtra)



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

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}



# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type



cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type))


#####


# Import the model for the selected disease and the drug target type
model <- readRDS(file = paste0("OutputFiles/Predictive_model_2/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


# Read the drug combinations for validation
valid_drugCombs_cat <- readRDS(file = paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))
valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "comb_name", "class_EffAdv", 
                                               "Drug1_indications", "Drug2_indications", 
                                               "DDI_description", "DDI_type")]


# Read the FGSEA results
fgsea_result <- readRDS(file = paste0("OutputFiles/Validation_data_2/Features/fgseaNES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


# Read the important features 
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"))


#####


# Extract the prediction data
predict_data <- fgsea_result[row.names(fgsea_result) %in% selected_features$feature, 
                             colnames(fgsea_result) %in% valid_drugCombs_cat$comb_name, 
                             drop = FALSE]
predict_data <- as.data.frame(t(predict_data))


# Check if all the features to be used are actually in the data
if(!all(names(model$models) %in% colnames(predict_data))){
  stop("Missing features")
}


#####
  
  
# Generate predictions

predictions_df <- data.frame("drugCombs" = row.names(predict_data))

for(feature_name in names(model$models)){
  predict_data_select <- predict_data[, feature_name, drop = FALSE]
  colnames(predict_data_select) <- "feature"
  predict_data_select$predicted_probability <- predict(model$models[[feature_name]], newdata = predict_data_select, type = "response")
  
  tmp1 <- predict_data_select[, c("predicted_probability"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  predictions_df <- full_join(predictions_df, tmp1, by = "drugCombs")
  rm(tmp1)
}

predictions_df <- column_to_rownames(predictions_df, "drugCombs")
predictions_df$final_predicted_probability <- rowMeans(predictions_df)
predictions_df$final_predicted_category <- ifelse(predictions_df$final_predicted_probability > model$probability_threshold, "Eff", "Adv")


# Merge the actual classes
predict_result <- merge(x = valid_drugCombs_cat, 
                        y = predictions_df[, c("final_predicted_probability", "final_predicted_category")], 
                        by.x = "comb_name", 
                        by.y = 0)

if(!dir.exists("OutputFiles/Predictive_model_2/")){ dir.create("OutputFiles/Predictive_model_2/", recursive = TRUE) }
write.csv(predict_result, file = paste0("OutputFiles/Predictive_model_2/val2_predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


#####


# Calculate prediction accuracy for each DDI type

predict_result$class_EffAdv <- factor(predict_result$class_EffAdv, levels = c("Eff", "Adv"))
predict_result$final_predicted_category <- factor(predict_result$final_predicted_category, levels = c("Eff", "Adv"))

predict_metrics <- data.frame("DDI_type" = "all", 
                              "Specificity" = specificity_vec(truth = predict_result$class_EffAdv, estimate = predict_result$final_predicted_category), 
                              "Number_of_combinations" = nrow(predict_result))

for(DDI_type in unique(predict_result$DDI_type)){
  predict_result_select <- predict_result[predict_result$DDI_type %in% DDI_type, ]
  tmp1 <- data.frame("DDI_type" = DDI_type,  
                     "Specificity" = specificity_vec(truth = predict_result_select$class_EffAdv, estimate = predict_result_select$final_predicted_category), 
                     "Number_of_combinations" = nrow(predict_result_select))
  predict_metrics <- rbind(predict_metrics, tmp1)
}


if(!dir.exists("OutputFiles/Predictive_model_2/")){ dir.create("OutputFiles/Predictive_model_2/", recursive = TRUE) }
write.csv(predict_metrics, file = paste0("OutputFiles/Predictive_model_2/val2_predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


print(warnings())