set.seed(5081)



# Script to plot the mean NES from FGSEA for EA drug category



# load libraries
library(unixtools)
library(tidyverse)
library(ggpubr)



# Set the parameters
drug_target_type <- "known"
feature_type <- "combinedEfficacySafety"
drugComb_category_type <- "EA"
proximity_type <- "Separation"


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

# plot for each disease
plot_list <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the network proximity data result
  proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
  proximity_result <- proximity_result[[feature_type]][[proximity_type]]
  proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
  
  
  # Read the drug combination category
  plot_col <- "class_EffAdv"
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
  drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), ]
  
  # Get the plot data
  plot_data <- proximity_result
  plot_data <- plot_data[, colnames(plot_data) %in% drugCombs_cat$comb_name]
  
  
  # Find the features with top variance and select the top three from each feature type
  feature_vars <- apply(plot_data, 1, var)
  feature_vars <- sort(feature_vars, decreasing = TRUE)
  tmp1 <- feature_vars[grep("^\\[ADR\\]", names(feature_vars))][1:3]
  tmp2 <- feature_vars[grep("^\\[DISEASE\\]", names(feature_vars))][1:3]
  feature_vars <- c(tmp1, tmp2)
  plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), ,  drop = FALSE]

  rm(list = c("tmp1", "tmp2"))
  
  
  # Add the drug combination categories to the plot data
  plot_data <- as.data.frame(t(plot_data))
  plot_data <- rownames_to_column(plot_data, "comb_name")
  
  plot_data$category <- drugCombs_cat[,c(plot_col), drop = TRUE][match(plot_data$comb_name, drugCombs_cat$comb_name)]
  plot_data <- plot_data %>% mutate_at("category", ~replace_na(., "unknown"))
  plot_data$category <- as.factor(plot_data$category)
  
  plot_data <- pivot_longer(data = plot_data, 
                            cols = colnames(plot_data)[!colnames(plot_data) %in% c("comb_name", "category")], 
                            cols_vary = "fastest", 
                            names_to = "feature", 
                            values_to = "value")
  
  
  # Calculate the mean and SD for each category for each feature
  plot_data <- plot_data %>% 
    group_by(feature, category) %>% 
    summarise(mean_value = mean(value), 
              se = sd(value) / sqrt(n()), 
              .groups = "drop")
  
  plot_data$feature <- factor(plot_data$feature, levels = names(feature_vars))
  
  
  # Plot
  plot_list[[disease]] <- ggplot(plot_data, aes(x = category, y = mean_value)) +
    geom_bar(width = 0.5, lwd = 0.1, stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_value - se, 
                      ymax = mean_value + se), 
                  width = 0.2, 
                  position = position_dodge(0.9)) +
    facet_wrap(~ feature, nrow = 1, labeller = label_wrap_gen(width = 30, 
                                                              multi_line = FALSE)) +
    labs(title = disease,
         y = "Mean distance +/- SE") +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(size = 6, margin = margin(1,1,1,1)),
          text = element_text(size = 8), 
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 7), 
          axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 5),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
}



# Save to file
if(!dir.exists("OutputFiles/Plots_publication/")){
  dir.create("OutputFiles/Plots_publication/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots_publication/meanProximity", proximity_type, "_xCancer_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".tiff"),
     width = 25,
     height = 25,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list, ncol = 1, nrow = 6)

dev.off()


print(warnings())