set.seed(5081)



# Script to plot the proximities (separation) to efficacy and safety gene set



# Load libraries
library(tidyverse)
library(ggpubr)



feature_type <- "combinedEfficacySafety"
proximity_type <- "Separation"





# Plot by disease
plot_list <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  plot_data <- data.frame()
  
  
  for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")){
    
    # Read the proximity measures
    proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
    proximity_result <- proximity_result[[feature_type]][[proximity_type]]
    proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
    
    
    # Read the drug combination category
    plot_col <- "class_EffAdv"
    drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
    
    stat_data <- proximity_result
    
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
      
      stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                              y = stat_data_select$value[stat_data_select$category == "Adv"])
      
      
      stat_res_final <- rbind(stat_res_final, data.frame("feature" = lib_name,
                                                         "W" = unname(stat_res$statistic),
                                                         "p_val" = stat_res$p.value))
    }
    
    stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.05, ]
    stat_res_final <- stat_res_final[order(stat_res_final$p_val, decreasing = FALSE), ]
    
    
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
    
    print(paste(disease, drug_target_type, safety_feature_select, efficacy_feature_select, sep = " - "))
    
  }
  
  plot_data$drug_target_type <- factor(x = plot_data$drug_target_type, levels = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")) 
  
  # Create the plots
  tmp2 <- ggplot(plot_data, aes(x = feature, y = value, fill = category)) +
    geom_boxplot(width = 0.5, 
                 lwd = 0.1,
                 outlier.shape = 3, 
                 outlier.size = 1) +
    facet_grid(rows = vars(plot_data$Disease), cols = vars(plot_data$drug_target_type), scales = "free_x") +
    scale_x_discrete(labels = function(x) scales::label_wrap(30)(x)) +
    stat_summary(fun = "mean",
                 geom = "point",
                 color = "red",
                 size = 0.5, 
                 position = position_dodge(width = 0.5), 
                 show.legend = FALSE) +
    # stat_summary(fun.data = "mean_se",
    #              geom = "errorbar",
    #              color = "red",
    #              width = 0.4,
    #              linewidth = 0.1,
    #              position = position_dodge(width = 0.5),
    #              show.legend = FALSE)  +
    geom_pwc(method = "wilcox_test",
             group.by = "category",
             label = "p.signif", 
             hide.ns = TRUE, 
             vjust = 0.5, 
             color = "blue") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    labs(x = "Feature", 
         y = "Separation distance") + 
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
          text = element_text(size = 5), 
          plot.title = element_text(size = 5, hjust = 0.5, face = "plain"),
          axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 5),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
  plot_list[[disease]] <- tmp2
  
  
  if(!dir.exists(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/"))){
    dir.create(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/"))
  }
  
  tiff(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA_", disease, ".tiff"), 
       width = 20,
       height = 7,
       units = "cm", compression = "lzw", res = 1200)
  print(tmp2)
  dev.off()
  
}





# plot by drug target type
plot_list <- list()
signif_feature <- list()

for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")){
  
  plot_data <- data.frame()
  
  
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    # Read the proximity measures
    proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
    proximity_result <- proximity_result[[feature_type]][[proximity_type]]
    proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
    
    
    # Read the drug combination category
    plot_col <- "class_EffAdv"
    drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
    
    stat_data <- proximity_result
    
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
      
      stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                              y = stat_data_select$value[stat_data_select$category == "Adv"])
      
      
      stat_res_final <- rbind(stat_res_final, data.frame("feature" = lib_name,
                                                         "W" = unname(stat_res$statistic),
                                                         "p_val" = stat_res$p.value))
    }
    
    stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.05, ]
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
    
    print(paste(disease, drug_target_type, safety_feature_select, efficacy_feature_select, sep = " - "))
    
  }
  
  plot_data$drug_target_type <- factor(x = plot_data$drug_target_type, levels = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")) 
  
  # Create the plots
  tmp2 <- ggplot(plot_data, aes(x = feature, y = value, fill = category)) +
    geom_boxplot(width = 0.5, 
                 lwd = 0.1,
                 outlier.shape = 3, 
                 outlier.size = 1) +
    facet_grid(rows = vars(plot_data$drug_target_type), cols = vars(plot_data$Disease), scales = "free_x") +
    scale_x_discrete(labels = function(x) scales::label_wrap(30)(x)) +
    stat_summary(fun = "mean",
                 geom = "point",
                 color = "red",
                 size = 0.5, 
                 position = position_dodge(width = 0.5), 
                 show.legend = FALSE) +
    # stat_summary(fun.data = "mean_se",
    #              geom = "errorbar",
    #              color = "red",
    #              width = 0.4,
    #              linewidth = 0.1,
    #              position = position_dodge(width = 0.5),
    #              show.legend = FALSE)  +
    geom_pwc(method = "wilcox_test",
             group.by = "category",
             label = "p.signif", 
             hide.ns = TRUE, 
             vjust = 0.5, 
             color = "blue") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    labs(x = "Feature", 
         y = "Separation distance") + 
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
          text = element_text(size = 5), 
          plot.title = element_text(size = 5, hjust = 0.5, face = "plain"),
          axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 5),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
  plot_list[[disease]] <- tmp2
  
  
  if(!dir.exists(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/"))){
    dir.create(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/"))
  }
  
  tiff(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA_", drug_target_type, ".tiff"), 
       width = 20,
       height = 7,
       units = "cm", compression = "lzw", res = 1200)
  print(tmp2)
  dev.off()
  
}



signif_feature <- unlist(signif_feature, recursive = FALSE)
signif_feature <- bind_rows(signif_feature, .id = "disease")
signif_feature <- separate(signif_feature, col = "disease", into = c("drug_target_type", "disease"), sep = "\\.")


write.csv(signif_feature, paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_EA/Proximity", proximity_type, "_", feature_type, "_WilcoxTest_rest.csv"), row.names = FALSE)


print(warnings())