set.seed(5081)



# Use the trained model to predict the labels on the validation data 2a



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
# Extract the feature thresholds
model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
feature_threshold <- model$feature_threshold


# Read the drug combinations for validation
valid_drugCombs_cat <- readRDS(file = paste0("InputFiles/Validation_data_2a/drugCombs_validation2a_", disease, ".rds"))
valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "comb_name", "class_EffAdv", 
                                               "Drug1_indications", "Drug2_indications", 
                                               "DDI_description", "DDI_type")]


# Read the FGSEA results
fgsea_result <- readRDS(file = paste0("OutputFiles/Validation_data_2a/Features/fgseaNES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


# Check if all the features needed for prediction are present in the FGSEA result
if(!all(feature_threshold$feature %in% row.names(fgsea_result))){
  stop("Missing feature in the input FGSEA results", call. = TRUE)
}


#####


# Extract the prediction data
predict_data <- fgsea_result[row.names(fgsea_result) %in% feature_threshold$feature, 
                             colnames(fgsea_result) %in% valid_drugCombs_cat$comb_name, 
                             drop = FALSE]
predict_data <- as.data.frame(t(predict_data))


#####


# Assign category based on each feature
for(feature_name in feature_threshold$feature){
  
  best_feature_threshold <- feature_threshold[feature_threshold$feature == feature_name, "threshold"]
  
  if(grepl("^\\[DISEASE\\]", feature_name)){
    predict_data[, feature_name] <- ifelse(predict_data[, feature_name] >= best_feature_threshold, 1, -1)
  }
  if(grepl("^\\[ADR\\]", feature_name)){
    predict_data[, feature_name] <- ifelse(predict_data[, feature_name] >= best_feature_threshold, -1, 1)
  }
  
}


# Compile and generate final classification
predict_data <- predict_data %>%
  mutate( efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
          safety_score = rowMeans(select(., starts_with("[ADR]"))) ) %>% 
  mutate( final_score = rowMeans(select(., c(efficacy_score, safety_score))) )

predict_data$final_predicted_category <- ifelse(predict_data$final_score > 0, "Eff", "Adv")
predict_data$final_predicted_category <- factor(predict_data$final_predicted_category, levels = c("Eff", "Adv"))


# Merge the actual classes
predict_result <- merge(x = valid_drugCombs_cat, 
                        y = predict_data[, c("efficacy_score", "safety_score", "final_score", "final_predicted_category")], 
                        by.x = "comb_name", 
                        by.y = 0)


#####


# Read the drug info from Drug Bank
DrugBank_drug_info <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_info <- DrugBank_drug_info$drugs$general_information
DrugBank_drug_info <- DrugBank_drug_info[, c("primary_key", "name")]
colnames(DrugBank_drug_info) <- c("DrugBank_drug_ID", "name")


# Add annotations about drugs
predict_result <- predict_result %>%
  left_join(DrugBank_drug_info %>% rename_with(.cols = everything(),
                                               .fn = ~ paste0("Drug1_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_DrugBank_drug_ID")) %>%
  left_join(DrugBank_drug_info %>% rename_with(.cols = everything(),
                                               .fn = ~ paste0("Drug2_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_DrugBank_drug_ID")) 


# Read the ATC codes of the drugs from Drug Bank
DrugBank_drug_ATC <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_drug_ATC$drugs$atc_codes
colnames(DrugBank_drug_ATC) <- gsub("drugbank-id", "DrugBank_drug_ID", colnames(DrugBank_drug_ATC))


# Concatenate the ATC level 1 term for each drug into a single string
ATC_l1 <- DrugBank_drug_ATC[, c("DrugBank_drug_ID", "level_1", "code_1")]
ATC_l1 <- ATC_l1 %>% 
  group_by(DrugBank_drug_ID, level_1)  %>% 
  summarise(code_1 = paste(unique(code_1), collapse = "; "), .groups = "keep") %>% 
  mutate(ATC_level1 = paste(unique(level_1), " (", unique(code_1), ")", sep = "")) %>%
  group_by(DrugBank_drug_ID) %>% 
  summarise(ATC_level1 = paste(unique(ATC_level1), collapse = "; "))

# Add annotations about ATC level 1
predict_result <- predict_result %>%
  left_join(ATC_l1 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug1_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_DrugBank_drug_ID")) %>%
  left_join(ATC_l1 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug2_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_DrugBank_drug_ID")) 

# Concatenate the ATC level 2 term for each drug into a single string
ATC_l2 <- DrugBank_drug_ATC[, c("DrugBank_drug_ID", "level_2", "code_2")]
ATC_l2 <- ATC_l2 %>% 
  group_by(DrugBank_drug_ID, level_2)  %>% 
  summarise(code_2 = paste(unique(code_2), collapse = "; "), .groups = "keep") %>% 
  mutate(ATC_level2 = paste(unique(level_2), " (", unique(code_2), ")", sep = "")) %>%
  group_by(DrugBank_drug_ID) %>% 
  summarise(ATC_level2 = paste(unique(ATC_level2), collapse = "; "))

# Add annotations about ATC level 2
predict_result <- predict_result %>%
  left_join(ATC_l2 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug1_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_DrugBank_drug_ID")) %>%
  left_join(ATC_l2 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug2_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_DrugBank_drug_ID")) 

# Concatenate the ATC level 3 term for each drug into a single string
ATC_l3 <- DrugBank_drug_ATC[, c("DrugBank_drug_ID", "level_3", "code_3")]
ATC_l3 <- ATC_l3 %>% 
  group_by(DrugBank_drug_ID, level_3)  %>% 
  summarise(code_3 = paste(unique(code_3), collapse = "; "), .groups = "keep") %>% 
  mutate(ATC_level3 = paste(unique(level_3), " (", unique(code_3), ")", sep = "")) %>%
  group_by(DrugBank_drug_ID) %>% 
  summarise(ATC_level3 = paste(unique(ATC_level3), collapse = "; "))


# Add annotations about ATC codes
predict_result <- predict_result %>%
  left_join(ATC_l3 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug1_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_DrugBank_drug_ID")) %>%
  left_join(ATC_l3 %>% rename_with(.cols = everything(),
                                   .fn = ~ paste0("Drug2_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_DrugBank_drug_ID")) 


#####


if(!dir.exists("OutputFiles/Validation_data_2a/Predictions/")){ dir.create("OutputFiles/Validation_data_2a/Predictions/", recursive = TRUE) }
write.csv(predict_result, file = paste0("OutputFiles/Validation_data_2a/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


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


if(!dir.exists("OutputFiles/Validation_data_2a/Prediction_metrics/")){ dir.create("OutputFiles/Validation_data_2a/Prediction_metrics/", recursive = TRUE) }
write.csv(predict_metrics, file = paste0("OutputFiles/Validation_data_2a/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)



print(warnings())