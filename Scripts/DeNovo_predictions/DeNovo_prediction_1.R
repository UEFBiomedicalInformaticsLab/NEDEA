set.seed(5081)



# Use the trained model to predict the labels on the de novo data 1


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
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--select_top_combs"), type = "numeric", default = 10,
              help = "The final number of drug combinations to be selected. Default: 10", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call. = FALSE)
}



# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type
select_top_combs <- opt$select_top_combs


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type))
cat(paste0("\nNumber of drug combinations to select: ", select_top_combs))


#####


# Import the model for the selected disease and the drug target type
model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
feature_threshold <- model$feature_threshold


# Read the drug combinations for de novo predictions
denovo_drugCombs <- readRDS(file = paste0("InputFiles/DeNovo_data_1/drugCombs_denovo1_", disease, ".rds"))
denovo_drugCombs$comb_name <- paste(denovo_drugCombs$Drug1_DrugBank_id, denovo_drugCombs$Drug2_DrugBank_id, sep = "_")
denovo_drugCombs <- denovo_drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "comb_name", 
                                         "Drug1_EMA", "Drug1_FDA", "Drug1_WHO", 
                                         "Drug1_Year", "Drug1_Generic", "Drug1_Product", "Drug1_Indications", 
                                         "Drug2_EMA", "Drug2_FDA", "Drug2_WHO", 
                                         "Drug2_Year", "Drug2_Generic", "Drug2_Product", "Drug2_Indications", 
                                         "DDI_description" )]


# Read the FGSEA results
denovo_fgsea_result <- readRDS(file = paste0("OutputFiles/DeNovo_data_1/Features/fgseaNES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


# Check if all the features needed for prediction are present in the FGSEA result
if(!all(feature_threshold$feature %in% row.names(denovo_fgsea_result))){
  stop("Missing feature in the input FGSEA results", call. = TRUE)
}


#####


# Extract the prediction data
predict_data <- denovo_fgsea_result[row.names(denovo_fgsea_result) %in% feature_threshold$feature, 
                             colnames(denovo_fgsea_result) %in% denovo_drugCombs$comb_name, 
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
predict_result <- merge(x = denovo_drugCombs, 
                        y = predict_data[, c("efficacy_score", "safety_score", "final_score", "final_predicted_category")], 
                        by.y = 0, by.x = "comb_name")

if(!dir.exists("OutputFiles/DeNovo_data_1/Predictions/")){ dir.create("OutputFiles/DeNovo_data_1/Predictions/", recursive = TRUE) }
write.csv(predict_result, file = paste0("OutputFiles/DeNovo_data_1/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


#####


# Identify the drugs to prioritize
priority_drugCombs <- as.data.frame(t(denovo_fgsea_result))
priority_drugCombs <- priority_drugCombs[row.names(priority_drugCombs) %in% predict_result[predict_result$final_score >= 1, "comb_name"], ]


# # Select drugs by comparing the mean efficacy score and mean safety score
# priority_drugCombs$mean_efficacy_NES <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, mean)
# priority_drugCombs$mean_safety_NES <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, mean)
# priority_drugCombs$diff_of_means <- priority_drugCombs$mean_efficacy_NES - priority_drugCombs$mean_safety_NES
# tmp1 <- sort(unique(priority_drugCombs$diff_of_means), decreasing = TRUE)[1:select_top_combs]
# priority_drugCombs$isSelected_byDiffOfMean <- ifelse(priority_drugCombs$diff_of_means %in% tmp1, TRUE, FALSE)
# rm(tmp1)


# Select drugs by comparing the max efficacy score and max safety score
priority_drugCombs$max_efficacy_NES <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, max)
priority_drugCombs$max_safety_NES <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, max)
priority_drugCombs$diff_of_maxs <- priority_drugCombs$max_efficacy_NES - priority_drugCombs$max_safety_NES
priority_drugCombs$which_max_efficacy_NES <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
priority_drugCombs$which_max_safety_NES <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
tmp1 <- sort(unique(priority_drugCombs$diff_of_maxs), decreasing = TRUE)[1:select_top_combs]
priority_drugCombs$isSelected_byDiffOfMax <- ifelse(priority_drugCombs$diff_of_maxs %in% tmp1, TRUE, FALSE)
rm(tmp1)


# # Select drugs by comparing the ratio of max efficacy score and max safety score
# priority_drugCombs$max_efficacy_NES <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, max)
# priority_drugCombs$max_safety_NES <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, max)
# priority_drugCombs$ratio_of_maxs <- priority_drugCombs$max_efficacy_NES / priority_drugCombs$max_safety_NES
# priority_drugCombs$which_max_efficacy_NES <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
# priority_drugCombs$which_max_safety_NES <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
# priority_drugCombs$isSelected_byRatioOfMax <- ifelse(priority_drugCombs$ratio_of_maxs >= 2, TRUE, FALSE)


# Save list of priority drug combinations
# priority_drugCombs <- priority_drugCombs[, c("mean_efficacy_NES", "mean_safety_NES", "diff_of_means", 
#                                              "max_efficacy_NES", "max_safety_NES", "diff_of_maxs", "ratio_of_maxs",
#                                              "which_max_efficacy_NES", "which_max_safety_NES", 
#                                              "isSelected_byDiffOfMean", "isSelected_byDiffOfMax", "isSelected_byRatioOfMax")]

priority_drugCombs <- priority_drugCombs[, c("max_efficacy_NES", "max_safety_NES", "diff_of_maxs",
                                             "which_max_efficacy_NES", "which_max_safety_NES", 
                                             "isSelected_byDiffOfMax")]
priority_drugCombs <- merge(x = predict_result,
                            y = priority_drugCombs, 
                            by.x = "comb_name", by.y = 0)
# priority_drugCombs <- priority_drugCombs[priority_drugCombs$isSelected_byDiffOfMean == "TRUE" | priority_drugCombs$isSelected_byDiffOfMax == "TRUE" | priority_drugCombs$isSelected_byRatioOfMax == "TRUE" , ]
priority_drugCombs <- priority_drugCombs[priority_drugCombs$isSelected_byDiffOfMax == "TRUE", ]
priority_drugCombs <- priority_drugCombs[order(priority_drugCombs$diff_of_maxs, decreasing = TRUE), ]


if(!dir.exists("OutputFiles/DeNovo_data_1/Priority_drug_combinations/")){ dir.create("OutputFiles/DeNovo_data_1/Priority_drug_combinations/", recursive = TRUE) }
write.csv(priority_drugCombs, file = paste0("OutputFiles/DeNovo_data_1/Priority_drug_combinations/priorityDrugCombs_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


#####


# Read the important features and select the top to plot
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"))

safety_feature_select <- selected_features[grep("^\\[ADR\\]", selected_features$feature), ]
if(nrow(safety_feature_select) > 0){
  safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)][1]
}else{ safety_feature_select <- c() }

efficacy_feature_select <- selected_features[grep("^\\[DISEASE\\]", selected_features$feature), ]
if(nrow(efficacy_feature_select) > 0){
  efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)][1]
}else{ efficacy_feature_select <- c() }


# Prepare the data for plotting
plot_data <- as.data.frame(t(denovo_fgsea_result))
plot_data <- plot_data[, colnames(plot_data) %in% c(safety_feature_select, efficacy_feature_select)]

colnames(plot_data)[colnames(plot_data) %in% efficacy_feature_select] <- "F1"
colnames(plot_data)[colnames(plot_data) %in% safety_feature_select] <- "F2"

x_axis_label = str_wrap(efficacy_feature_select, 40)
y_axis_label = str_wrap(safety_feature_select, 40)

plot_data$final_score <- predict_result$final_score[match(row.names(plot_data), predict_result$comb_name)]

plot_data <- plot_data %>% 
  rownames_to_column("comb_name") %>% 
  left_join(priority_drugCombs[, c("comb_name", "isSelected_byDiffOfMax")], 
            by = "comb_name")

# plot_data[is.na(plot_data$isSelected_byDiffOfMean), "isSelected_byDiffOfMean"] <- FALSE
plot_data[is.na(plot_data$isSelected_byDiffOfMax), "isSelected_byDiffOfMax"] <- FALSE
# plot_data[is.na(plot_data$isSelected_byRatioOfMax), "isSelected_byRatioOfMax"] <- FALSE


#####


if(!dir.exists("OutputFiles/Plots/DeNovo_predictions/")){
  dir.create("OutputFiles/Plots/DeNovo_predictions/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/DeNovo_predictions/plot_DeNovo_1_predictions_combinedEfficacySafety_", disease, "_", drug_target_type, ".tiff"),
     width = 7, height = 6,
     units = "cm", compression = "lzw", res = 1200)


ggplot() +
  geom_point(data = plot_data[plot_data$isSelected_byDiffOfMax == "FALSE", ], 
             mapping = aes(x = F1, y = F2, color = final_score),  
             size = 0.5, 
             stroke = 0.1,
             shape = 3, 
             alpha = 0.85) +
  # geom_point(data = plot_data[plot_data$isSelected_byDiffOfMean == "TRUE", ], # highlight drugs by mean
  #            mapping = aes(x = F1, y = F2, color = final_score), 
  #            size = 0.5, 
  #            stroke = 0.2,   
  #            shape = 1) +
  geom_point(data = plot_data[plot_data$isSelected_byDiffOfMax == "TRUE", ], # highlight drugs by max
             mapping = aes(x = F1, y = F2, color = final_score), 
             size = 0.6, 
             stroke = 0.2,   
             shape = 2) +
  # geom_point(data = plot_data[plot_data$isSelected_byRatioOfMax == "TRUE", ], # highlight drugs by ratio
  #            mapping = aes(x = F1, y = F2, color = final_score), 
  #            size = 0.6, 
  #            stroke = 0.4,   
  #            shape = 5) +
  geom_hline(yintercept = feature_threshold[feature_threshold$feature %in% safety_feature_select, "threshold"], 
             linetype="dotted",
             color = "gray",
             linewidth = 0.25) +
  geom_vline(xintercept = feature_threshold[feature_threshold$feature %in% efficacy_feature_select, "threshold"], 
             linetype="dotted",
             color = "gray",
             linewidth = 0.25) +
  scale_color_gradient2(low = "#FF2015", mid = "#FFD600", high = "#3ACE3A", midpoint = 0) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "right",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 2.5),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0, "cm")  ) +
  labs(title = disease,
       x = x_axis_label,
       y = y_axis_label,
       color = "Final score") 

dev.off()



print(warnings())