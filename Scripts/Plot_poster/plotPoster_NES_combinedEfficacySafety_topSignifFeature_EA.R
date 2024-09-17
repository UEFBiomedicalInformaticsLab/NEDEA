set.seed(5081)



# Script to plot the NES of efficacy and safety gene set



# Load libraries
library(tidyverse)
library(ggpubr)


# Set the feature type to use
feature_type <- "combinedEfficacySafety"


#####


# Plot by disease
plot_list <- list()
signif_feature <- list()

for(drug_target_type in c("known")){
  
  plot_data <- data.frame()
  
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    # Read the FGSEA result
    fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
    fgsea_result <- fgsea_result[[feature_type]]
    
    # Read the drug combination category
    plot_col <- "class_EffAdv"
    drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
    
    stat_data <- fgsea_result
    
    # Add the drug combination categories to the stat data
    stat_data <- as.data.frame(t(stat_data))
    stat_data <- rownames_to_column(stat_data, "comb_name")
    
    stat_data$category <- drugCombs_cat[,c(plot_col), drop = TRUE][match(stat_data$comb_name, drugCombs_cat$comb_name)]
    stat_data <- stat_data %>% mutate_at("category", ~replace_na(., "unknown"))
    stat_data <- stat_data[stat_data$category != "unknown", ]
    stat_data$category <- as.factor(stat_data$category)
    
    stat_data <- pivot_longer(data = stat_data, 
                              cols = colnames(stat_data)[!colnames(stat_data) %in% c("comb_name", "category")], 
                              cols_vary = "fastest", 
                              names_to = "feature", 
                              values_to = "value")
    
    
    # Calculate the statistical difference
    stat_res_final <- data.frame()
    for(lib_name in unique(stat_data$feature)){
      
      stat_data_select <- stat_data[stat_data$feature == lib_name, ]
      
      
      if(grepl("^\\[ADR\\]", lib_name)){
        stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                                y = stat_data_select$value[stat_data_select$category == "Adv"], 
                                alternative = "less")
      }
      
      if(grepl("^\\[DISEASE\\]", lib_name)){
        stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                                y = stat_data_select$value[stat_data_select$category == "Adv"], 
                                alternative = "greater")
      }
      
      
      
      stat_res_final <- rbind(stat_res_final, data.frame("feature" = lib_name,
                                                         "W" = unname(stat_res$statistic),
                                                         "p_val" = stat_res$p.value))
    }
    
    
    stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.001, ]
    stat_res_final <- stat_res_final[order(stat_res_final$p_val, decreasing = FALSE), ]
    
    signif_feature[[drug_target_type]][[disease]] <- stat_res_final
    
    # Select the feature
    safety_feature_select <- stat_res_final[grep("^\\[ADR\\]", stat_res_final$feature), ]
    if(nrow(safety_feature_select) > 0){
      safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)]
    }else{ safety_feature_select <- c() }
    
    efficacy_feature_select <- stat_res_final[grep("^\\[DISEASE\\]", stat_res_final$feature), ]
    if(nrow(efficacy_feature_select) > 0){
      efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)]
    }else{ efficacy_feature_select <- c() }
    
    tmp1 <- stat_data[stat_data$feature %in% c(safety_feature_select, efficacy_feature_select), ]
    tmp1$Disease <- disease
    tmp1$drug_target_type <- drug_target_type
    
    plot_data <- rbind(plot_data, tmp1)
    
    
  }
  
  tmp2 <- ggplot(plot_data, aes(x = feature, y = value, fill = category)) +
    geom_boxplot(width = 0.75, 
                 lwd = 0.25,
                 outlier.shape = 3, 
                 outlier.size = 2) +
    # facet_grid(rows = vars(plot_data$drug_target_type), cols = vars(plot_data$Disease), scales = "free_x") +
    facet_grid(cols = vars(plot_data$Disease), scales = "free_x", labeller = labeller(.cols = function(x){ gsub("Cancer$", " Cancer", x) })) +
    scale_x_discrete(labels = function(x) scales::label_wrap(30)(x)) +
    stat_summary(fun = "mean",
                 geom = "point",
                 color = "blue",
                 size = 2, 
                 position = position_dodge(width = 0.75), 
                 show.legend = FALSE) +
    geom_pwc(method = "wilcox_test",
             group.by = "category",
             label = "p.signif", 
             hide.ns = TRUE, 
             vjust = 0.5, 
             color = "blue",
             label.size = 5) +
    scale_fill_manual(values = c("Adv" = "#FF6961", "Eff" = "#77DD77"), labels = c("Adv" = "Adverse", "Eff" = "Effective")) +
    labs(x = "Feature", 
         y = "NES",
         fill = "Drug combination type:") + 
    theme(plot.margin = margin(l = 0.4, r = 0.1, t = 0.1, b = 0.1, unit = "cm"),
          panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(size = 14, margin = margin(2,1,2,1)),
          text = element_text(size = 14), 
          plot.title = element_text(size = 14, hjust = 0.5, face = "plain"),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.title = element_text(margin = margin(r = 5)),
          legend.key = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.key.spacing.x = unit(0.25, "cm"),
          legend.text = element_text(size = 12, margin = margin(l = 1)),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.25, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25))
  
  if(!dir.exists(paste0("OutputFiles/Plots_poster/NES_", feature_type, "_boxplot_topSignifFeatures_EA/"))){
    dir.create(paste0("OutputFiles/Plots_poster/NES_", feature_type, "_boxplot_topSignifFeatures_EA/"), recursive = TRUE)
  }
  
  tiff(paste0("OutputFiles/Plots_poster/NES_", feature_type, "_boxplot_topSignifFeatures_EA/NES_", feature_type, "_boxplot_topSignifFeatures_EA_", drug_target_type, ".tiff"),
       width = 35,
       height = 12,
       units = "cm", compression = "lzw", res = 1200)
  print(tmp2)
  dev.off()
  
}


print(warnings())