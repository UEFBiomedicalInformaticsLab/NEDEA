set.seed(5081)



# Script to plot the applicibility domain for the validation data set 1 for all cancers



# Load libraries
library(tidyverse)
library(ggpubr)



plot_list <- list()
count_df <- data.frame()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the training set drug combination 
  train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]
  
  
  # Read the validation set drug combination 
  valid_drugCombs_cat <- readRDS(file = paste0("InputFiles/Validation_data/Validation1_drugCombs_cat_effVadv_", disease, ".rds"))
  

  # Read the FGSEA results
  fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_known.rds"))
  fgsea_result <- fgsea_result[["combinedEfficacySafety"]]
  
  
  # Read the important features and select the top to plot
  selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_known.csv"))
  
  safety_feature_select <- selected_features[grep("^\\[ADR\\]", selected_features$feature), ]
  if(nrow(safety_feature_select) > 0){
    safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)][1]
  }else{ safety_feature_select <- c() }
  
  
  efficacy_feature_select <- selected_features[grep("^\\[DISEASE\\]", selected_features$feature), ]
  if(nrow(efficacy_feature_select) > 0){
    efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)][1]
  }else{ efficacy_feature_select <- c() }
  
  
  # Get the FGSEA results for the train and test data
  train_fgsea_result <- fgsea_result[row.names(fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                     colnames(fgsea_result) %in% train_drugCombs_cat$comb_name]
  
  valid_fgsea_result <- fgsea_result[row.names(fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                     colnames(fgsea_result) %in% valid_drugCombs_cat$comb_name]
  
  
  train_fgsea_result <- as.data.frame(t(train_fgsea_result))
  train_fgsea_result$category <- train_drugCombs_cat$class_EffAdv[match(row.names(train_fgsea_result), train_drugCombs_cat$comb_name)]
  train_fgsea_result$data_group <- "Training"
  
  valid_fgsea_result <- as.data.frame(t(valid_fgsea_result))
  valid_fgsea_result$category <- valid_drugCombs_cat$class_EffAdv[match(row.names(valid_fgsea_result), valid_drugCombs_cat$comb_name)]
  valid_fgsea_result$data_group <- "Validation"

  
  plot_data <- rbind(train_fgsea_result, valid_fgsea_result)
  
  colnames(plot_data)[colnames(plot_data) %in% efficacy_feature_select] <- "F1"
  colnames(plot_data)[colnames(plot_data) %in% safety_feature_select] <- "F2"
  
  x_axis_label = str_wrap(efficacy_feature_select, 40)
  y_axis_label = str_wrap(safety_feature_select, 40)
  

  plot_list[[disease]] <- ggplot() +
    geom_point(data = plot_data, 
               mapping = aes(x = F1, y = F2, color = category, shape = data_group),  
               size = 0.4, 
               stroke = 0.1) +
    scale_shape_manual(values = c("Training" = 1, "Validation" = 3)) + 
    scale_color_manual(values = c("Eff" = "#0000FF", "Adv" = "#FF0000")) + 
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          text = element_text(size = 3),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 3, margin = margin(2, 1, 1, 1, "pt")),
          plot.margin = margin(1, 2, 1, 2, "pt"), 
          axis.title = element_text(size = 2.5), 
          axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.1),
          legend.position = "right",
          legend.key = element_rect(fill = NA), 
          legend.key.size = unit(0.1, "cm"),
          legend.title = element_text(size = 2.5),
          legend.text = element_text(size = 2),
          legend.margin = margin(1,1,1,1),
          legend.spacing = unit(0.25, "cm"), 
          legend.background =  element_rect(color = "black", linewidth = 0.05)) +
    labs(title = disease,
         x = x_axis_label,
         y = y_axis_label,
         color = "Category", 
         shape = "Group")
  
  
  
  tmp1 <- data.frame("Disease" = disease, 
                     Train_Eff = nrow(plot_data[plot_data$category == "Eff" & plot_data$data_group == "Training", ]), 
                     Train_Adv = nrow(plot_data[plot_data$category == "Adv" & plot_data$data_group == "Training", ]),
                     Valid_Eff = nrow(plot_data[plot_data$category == "Eff" & plot_data$data_group == "Validation", ]),
                     Valid_Adv = nrow(plot_data[plot_data$category == "Adv" & plot_data$data_group == "Validation", ])
                     )
  
  count_df <- rbind(count_df, tmp1)
  
}


if(!dir.exists("OutputFiles/Plots_publication/")){
  dir.create("OutputFiles/Plots_publication/", recursive = TRUE)
}
tiff("OutputFiles/Plots_publication/Validation1_AD_xCancer.tiff",
     width = 8, height = 6,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

dev.off()


write.csv(count_df, "OutputFiles/Plots_publication/Validation1_counts.csv", row.names = FALSE)



print(warnings())