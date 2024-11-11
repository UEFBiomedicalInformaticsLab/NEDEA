set.seed(5081)



# Use the trained model to predict the labels on the validation data 2b



# Load libraries
library(unixtools)
library(optparse)
library(ComplexHeatmap)
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
valid_drugCombs_cat <- readRDS(file = paste0("InputFiles/Validation_data_2b/drugCombs_validation2b_", disease, ".rds"))
valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "comb_name", "class_EffAdv", 
                                               "Drug1_indications", "Drug2_indications", 
                                               "DDI_description", "DDI_type")]


# Read the FGSEA results
fgsea_result <- readRDS(file = paste0("OutputFiles/Validation_data_2b/Features/fgseaNES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


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


if(!dir.exists("OutputFiles/Validation_data_2b/Predictions/")){ dir.create("OutputFiles/Validation_data_2b/Predictions/", recursive = TRUE) }
write.csv(predict_result, file = paste0("OutputFiles/Validation_data_2b/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


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


if(!dir.exists("OutputFiles/Validation_data_2b/Prediction_metrics/")){ dir.create("OutputFiles/Validation_data_2b/Prediction_metrics/", recursive = TRUE) }
write.csv(predict_metrics, file = paste0("OutputFiles/Validation_data_2b/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


#####


# Plot heat map of the NES of the drug combinations with the results
plot_data <- t(fgsea_result) 

# Create annotation for the rows
left_annot_color <- list(labelled_class = c("Eff" = "#77DD77", "Adv" = "#FF6961", "Unk" = "#808080"), 
                         predicted_class = c("Eff" = "#77DD77", "Adv" = "#FF6961", "Unk" = "#808080"))

left_annot <- predict_result[, c("comb_name", "class_EffAdv", "final_predicted_category")]
colnames(left_annot) <- c("comb_name", "labelled_class", "predicted_class")
left_annot <- left_annot[match(row.names(plot_data), left_annot$comb_name), ]
row.names(left_annot) <- NULL
left_annot <- column_to_rownames(left_annot, "comb_name")
left_annot <- HeatmapAnnotation(which = "row",
                                df = left_annot,
                                col = left_annot_color,
                                simple_anno_size = unit(0.25, "cm"),
                                annotation_name_gp = gpar(fontsize = 4),
                                annotation_name_rot = 45, 
                                annotation_legend_param = list(title_gp = gpar(fontsize = 4),
                                                               labels_gp = gpar(fontsize = 4),
                                                               grid_height = unit(0.25, "cm"),
                                                               grid_width = unit(0.25, "cm")
                                ))


# Create annotation for the columns
top_annot_color <- list(Feature_type = c("ADR" = "#FF6961", "DISEASE" = "#77DD77"),
                        selected_feature = c("Yes" = "#228B22", "No" = "#FF0000"))

top_annot <- as.data.frame(colnames(plot_data))
colnames(top_annot) <- "Features"
top_annot$selected_feature <- ifelse(top_annot$Features %in% feature_threshold$feature, "Yes", "No")

top_annot$Feature_type <- gsub("^\\[(.*)\\] .+", "\\1", top_annot$Features)
top_annot$Features <- gsub("^\\[(.*)\\] ", "", top_annot$Features)
# top_annot <- top_annot[order(top_annot$Features), ]
row.names(top_annot) <- NULL
top_annot$Features<- str_wrap(top_annot$Features, 30)
top_annot <- column_to_rownames(top_annot, "Features")
top_annot <- HeatmapAnnotation(which = "column",
                               df = top_annot,
                               col = top_annot_color,
                               simple_anno_size = unit(0.25, "cm"),
                               annotation_name_gp = gpar(fontsize = 4),
                               annotation_legend_param = list(title_gp = gpar(fontsize = 4),
                                                              labels_gp = gpar(fontsize = 4),
                                                              grid_height = unit(0.25, "cm"),
                                                              grid_width = unit(0.25, "cm")
                               ))


# Define color function
col_fun <- circlize::colorRamp2(breaks = c(min(plot_data, na.rm = TRUE), max(plot_data, na.rm = TRUE)),
                                colors = c("#CCF9FF", "#0080BF"))

colnames(plot_data) <- gsub("^\\[(.*)\\] ", "", colnames(plot_data))
colnames(plot_data) <- str_wrap(colnames(plot_data), 30)

heatmap <- Heatmap(plot_data,
                   
                   # col = col_fun,
                   cluster_columns = FALSE,
                   
                   row_title = "Drug combinations",
                   column_title = "Features",
                   row_title_side = "left",
                   column_title_side = "bottom",
                   row_title_gp = gpar(fontsize = 5, face = "bold"),
                   column_title_gp = gpar(fontsize = 5, face = "bold"),
                   
                   row_dend_gp = gpar(lwd = 0.5),
                   column_dend_gp = gpar(lwd = 0.5),
                   
                   left_annotation = left_annot,
                   top_annotation = top_annot,
                   
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize = 1),
                   column_names_gp = gpar(fontsize = 4, linebreak = TRUE),
                   column_names_rot = 45, 
                   
                   heatmap_legend_param = list(title = "NES",
                                               title_gp = gpar(fontsize = 4),
                                               labels_gp = gpar(fontsize = 4),
                                               legend_height = unit(4, "cm"),
                                               legend_width = unit(0.1, "cm")
                   ))



if(!dir.exists("OutputFiles/Plots/Validation_data_heatmaps/")){
  dir.create("OutputFiles/Plots/Validation_data_heatmaps/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/Validation_data_heatmaps/plot_validation2b_heatmap_combinedEfficacySafety_", disease, "_", drug_target_type, ".tiff"),
     width = 25, height = 21,
     units = "cm", compression = "lzw", res = 1200)

draw(heatmap)

dev.off()


#####


print(warnings())