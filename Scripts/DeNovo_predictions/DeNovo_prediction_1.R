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
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
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


# Read the important features 
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"))


#####


# Extract the prediction data
predict_data <- denovo_fgsea_result[row.names(denovo_fgsea_result) %in% selected_features$feature, 
                                    colnames(denovo_fgsea_result) %in% denovo_drugCombs$comb_name]
predict_data <- as.data.frame(t(predict_data))


# Generate predictions
predict_result <- predict(object = model, newdata = predict_data, type = "prob")
colnames(predict_result) <- paste0("predicted_prob", colnames(predict_result))
predict_result$predicted_category <- predict(object = model, newdata = predict_data, type = "raw")


# Merge the actual classes
predict_result <- merge(y = predict_result, 
                        x = denovo_drugCombs, 
                        by.y = 0, by.x = "comb_name")
# predict_result$class_EffAdv <- factor(x = predict_result$class_EffAdv, levels = c("Eff", "Adv"))
predict_result$predicted_category <- factor(x = predict_result$predicted_category, levels = c("Eff", "Adv"))

if(!dir.exists("OutputFiles/DeNovo_data_1/Predictions/")){ dir.create("OutputFiles/DeNovo_data_1/Predictions/", recursive = TRUE) }
write.csv(predict_result, file = paste0("OutputFiles/DeNovo_data_1/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


#####


# Identify the drugs to prioritize
priority_drugCombs <- as.data.frame(t(denovo_fgsea_result))
priority_drugCombs <- priority_drugCombs[row.names(priority_drugCombs) %in% predict_result[predict_result$predicted_category == "Eff", "comb_name"], ]


# Select drugs by comparing the mean efficacy score and mean safety score
priority_drugCombs$mean_efficacy_score <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, mean)
priority_drugCombs$mean_safety_score <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, mean)
priority_drugCombs$diff_of_means <- priority_drugCombs$mean_efficacy_score  - priority_drugCombs$mean_safety_score
tmp1 <- sort((priority_drugCombs$diff_of_means), decreasing = TRUE)[1:select_top_combs]
priority_drugCombs$isSelected_byDiffOfMean <- ifelse(priority_drugCombs$diff_of_means %in% tmp1, TRUE, FALSE)
rm(tmp1)


# Select drugs by comparing the max efficacy score and max safety score
priority_drugCombs$max_efficacy_score <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, max)
priority_drugCombs$max_safety_score <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, max)
priority_drugCombs$diff_of_maxs <- priority_drugCombs$max_efficacy_score  - priority_drugCombs$max_safety_score
priority_drugCombs$which_max_efficacy_score <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
priority_drugCombs$which_max_safety_score <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 1, function(x){ paste(names(which(x == max(x))), collapse = "; ") })
tmp1 <- sort((priority_drugCombs$diff_of_maxs), decreasing = TRUE)[1:select_top_combs]
priority_drugCombs$isSelected_byDiffOfMax <- ifelse(priority_drugCombs$diff_of_maxs %in% tmp1, TRUE, FALSE)
rm(tmp1)


# Save list of priority drug combinations
priority_drugCombs <- priority_drugCombs[, c("mean_efficacy_score", "mean_safety_score", "diff_of_means", 
                                             "max_efficacy_score", "max_safety_score", "diff_of_maxs", 
                                             "which_max_efficacy_score", "which_max_safety_score", 
                                             "isSelected_byDiffOfMean", "isSelected_byDiffOfMax")]
priority_drugCombs <- merge(x = predict_result,
                            y = priority_drugCombs, 
                            by.x = "comb_name", by.y = 0)
priority_drugCombs <- priority_drugCombs[priority_drugCombs$isSelected_byDiffOfMean == "TRUE" | priority_drugCombs$isSelected_byDiffOfMax == "TRUE", ]


if(!dir.exists("OutputFiles/DeNovo_data_1/Priority_drug_combinations/")){ dir.create("OutputFiles/DeNovo_data_1/Priority_drug_combinations/", recursive = TRUE) }
write.csv(priority_drugCombs, file = paste0("OutputFiles/DeNovo_data_1/Priority_drug_combinations/priorityDrugCombs_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)


# tmp1 <- apply(priority_drugCombs[,grep("^\\[DISEASE\\]", colnames(priority_drugCombs))], 2, function(x){ rank(x) })
# tmp2 <- apply(priority_drugCombs[,grep("^\\[ADR\\]", colnames(priority_drugCombs))], 2, function(x){ rank(-x) })
# 
# 
# priority_drugCombs_ranks <- merge(tmp1, tmp2, by = 0)
# priority_drugCombs_ranks <- column_to_rownames(priority_drugCombs_ranks, "Row.names")
# priority_drugCombs_ranks$total_score <- apply(priority_drugCombs_ranks, 1, sum)
# 
# priority_drugCombs_ranks <- priority_drugCombs_ranks[, c("total_score"), drop = FALSE]


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

plot_data$predicted_probEff <- predict_result$predicted_probEff[match(row.names(plot_data), predict_result$comb_name)]

plot_data <- plot_data %>% 
  rownames_to_column("comb_name") %>% 
  left_join(priority_drugCombs[, c("comb_name", "isSelected_byDiffOfMean", "isSelected_byDiffOfMax")], 
            by = "comb_name")

plot_data[is.na(plot_data$isSelected_byDiffOfMean), "isSelected_byDiffOfMean"] <- FALSE
plot_data[is.na(plot_data$isSelected_byDiffOfMax), "isSelected_byDiffOfMax"] <- FALSE


#####


if(!dir.exists("OutputFiles/Plots/DeNovo_predictions/")){
  dir.create("OutputFiles/Plots/DeNovo_predictions/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/DeNovo_predictions/plot_DeNovo_1_predictions_combinedEfficacySafety_", disease, "_", drug_target_type, ".tiff"),
     width = 7, height = 6,
     units = "cm", compression = "lzw", res = 1200)


ggplot() +
  geom_point(data = plot_data[plot_data$isSelected_byDiffOfMean == "FALSE" & plot_data$isSelected_byDiffOfMax == "FALSE", ], 
             mapping = aes(x = F1, y = F2, color = predicted_probEff),  
             size = 0.5, 
             stroke = 0.1,
             shape = 3) +
  geom_point(data = plot_data[plot_data$isSelected_byDiffOfMean == "TRUE", ], # highlight drugs by mean
             mapping = aes(x = F1, y = F2, color = predicted_probEff), 
             size = 0.5, 
             stroke = 0.2,   
             shape = 1) +
  geom_point(data = plot_data[plot_data$isSelected_byDiffOfMax == "TRUE", ], # highlight drugs by mean
             mapping = aes(x = F1, y = F2, color = predicted_probEff), 
             size = 0.6, 
             stroke = 0.2,   
             shape = 5) +
  scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5) +
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
       color = "prob_Eff") 

dev.off()



print(warnings())