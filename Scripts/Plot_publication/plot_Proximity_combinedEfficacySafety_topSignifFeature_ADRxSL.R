set.seed(5081)



# Script to plot the proximities (separation) to efficacy and safety gene set



# Load libraries
library(tidyverse)
library(ggpubr)



feature_type <- "combinedEfficacySafety"
proximity_type <- "Separation"



for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  
  plot_list <- list()
  
  for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")){
    
    # Read the proximity measures
    proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
    proximity_result <- proximity_result[[feature_type]][[proximity_type]]
    proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
    
    
    # Read the drug combination category
    plot_col_1 <- "Syn_level"
    drugCombs_cat_1 <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
    drugCombs_cat_1$comb_name <- paste(drugCombs_cat_1$Drug1_DrugBank_id, drugCombs_cat_1$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat_1 <- drugCombs_cat_1[, c("comb_name", plot_col_1)]
    
    plot_col_2 <- paste0("ADR_", disease)
    drugCombs_cat_2 <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
    drugCombs_cat_2$comb_name <- paste(drugCombs_cat_2$Drug1_DrugBank_id, drugCombs_cat_2$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat_2 <- drugCombs_cat_2[, c("comb_name", plot_col_2)]
    
    stat_data <- proximity_result
    
    # Add the drug combination categories to the stat data
    stat_data <- as.data.frame(t(stat_data))
    stat_data <- rownames_to_column(stat_data, "comb_name")
    
    stat_data$category_SL <- drugCombs_cat_1[,c(plot_col_1), drop = TRUE][match(stat_data$comb_name, drugCombs_cat_1$comb_name)]
    stat_data$category_ADR <- drugCombs_cat_2[,c(plot_col_2), drop = TRUE][match(stat_data$comb_name, drugCombs_cat_2$comb_name)]
    
    stat_data <- stat_data %>% mutate_at("category_SL", ~replace_na(., "unknown"))
    stat_data <- stat_data %>% mutate_at("category_ADR", ~replace_na(., "unknown"))
    
    stat_data$category_SL <- as.factor(stat_data$category_SL)
    stat_data$category_ADR <- as.factor(stat_data$category_ADR)
    
    
    stat_data <- pivot_longer(data = stat_data, 
                              cols = colnames(stat_data)[!colnames(stat_data) %in% c("comb_name", "category_SL", "category_ADR")], 
                              cols_vary = "fastest", 
                              names_to = "feature", 
                              values_to = "value")
    
    # Calculate the statistical difference
    stat_res_final <- data.frame()
    for(lib_name in unique(stat_data$feature)){
      for(synergy_level in unique(stat_data$category_SL)){
        
        stat_data_select <- stat_data[stat_data$feature == lib_name & stat_data$category_SL == synergy_level, ]
        
        stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category_ADR == "adr_positive"],
                                y = stat_data_select$value[stat_data_select$category_ADR == "unknown"])
        
        
        stat_res_final <- rbind(stat_res_final, data.frame("feature" = lib_name,
                                                           "synergy_level" = synergy_level,
                                                           "W" = unname(stat_res$statistic),
                                                           "p_val" = stat_res$p.value))
      }
    }
    
    stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.001, ]
    stat_res_final <- stat_res_final[order(stat_res_final$p_val, decreasing = FALSE), ]
    
    
    # Select the feature
    # safety_feature_select <- list()
    # efficacy_feature_select <- list()
    
    
    for(synergy_level in levels(stat_data$category_SL)){
      
      plot_data <- data.frame()
      
      safety_feature_select <- stat_res_final[grep("^\\[ADR\\]", stat_res_final$feature), ]
      safety_feature_select <- safety_feature_select[safety_feature_select$synergy_level == synergy_level, ]
      if(nrow(safety_feature_select) > 0){
        safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)][1] # Used to eliminate situation where two or more terms have same pval
      }else{ safety_feature_select <- c() }
      
      efficacy_feature_select <- stat_res_final[grep("^\\[DISEASE\\]", stat_res_final$feature), ]
      efficacy_feature_select <- efficacy_feature_select[efficacy_feature_select$synergy_level == synergy_level, ]
      
      if(nrow(efficacy_feature_select) > 0){
        efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)][1]
      }else{ efficacy_feature_select <- c() }
      
      
      tmp1 <- stat_data[(stat_data$category_SL == synergy_level) & (stat_data$feature %in% c(safety_feature_select, efficacy_feature_select)), ]
      tmp1$Disease <- disease
      tmp1$Synergy_level <- synergy_level
      tmp1$drug_target_type <- drug_target_type
      
      plot_data <- tmp1
      
      print(paste(disease, synergy_level, drug_target_type, safety_feature_select, efficacy_feature_select, sep = " - "))
      # print(paste(disease, synergy_level, drug_target_type, sep = " - "))
      
      
      if(nrow(plot_data) > 0){
        # Create the plots
        tmp2 <- ggplot(plot_data, aes(x = feature, y = value, fill = category_ADR)) +
          geom_boxplot(width = 0.5, 
                       lwd = 0.1,
                       outlier.shape = 3, 
                       outlier.size = 1) +
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
          labs(title = paste0(drug_target_type, " [", synergy_level, "]")) +
          theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
                panel.grid = element_blank(),
                panel.spacing = unit(0.1, "cm"),
                strip.background = element_rect(color = "black", linewidth = 0.25,),
                strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
                text = element_text(size = 5), 
                plot.title = element_text(size = 5, hjust = 0.5, face = "plain"),
                axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5), 
                axis.ticks = element_line(colour = "black", linewidth = 0.2),
                legend.position = "bottom",
                legend.key = element_blank(),
                legend.key.size = unit(0.1, 'cm'),
                legend.text = element_text(size = 5),
                legend.margin = margin(1,1,1,1),
                legend.box.spacing = unit(0.1, 'cm'),
                legend.box.background = element_rect(colour = "black", linewidth = 0.1))
      }else{
        tmp2 <- ggplot() +
          labs(title = paste0(drug_target_type, " [", synergy_level, "]")) +
          theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
                panel.grid = element_blank(),
                panel.spacing = unit(0.1, "cm"),
                strip.background = element_rect(color = "black", linewidth = 0.25,),
                strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
                text = element_text(size = 5), 
                plot.title = element_text(size = 5, hjust = 0.5, face = "plain"),
                axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5), 
                axis.ticks = element_line(colour = "black", linewidth = 0.2),
                legend.position = "bottom",
                legend.key = element_blank(),
                legend.key.size = unit(0.1, 'cm'),
                legend.text = element_text(size = 5),
                legend.margin = margin(1,1,1,1),
                legend.box.spacing = unit(0.1, 'cm'),
                legend.box.background = element_rect(colour = "black", linewidth = 0.1))
      }
      
      plot_list[[drug_target_type]][[synergy_level]] <- tmp2
      
    }
  }
  
  
  tmp3 <- unlist(plot_list, recursive = FALSE)
  
  if(!dir.exists(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_ADRxSL/"))){
    dir.create(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_ADRxSL/"))
  }
  
  tiff(paste0("OutputFiles/Plots_publication/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_ADRxSL/Proximity", proximity_type, "_", feature_type, "_boxplot_topSignifFeatures_ADRxSL_", disease, ".tiff"), 
       width = 50,
       height = 40,
       units = "cm", compression = "lzw", res = 1200)
  print(ggpubr::ggarrange(plotlist = tmp3, ncol =11, nrow = 7, common.legend = TRUE))
  dev.off()
  
}



print(warnings())