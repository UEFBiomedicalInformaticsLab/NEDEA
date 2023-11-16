
# Load libraries
library(tidyverse)
library(ggpubr)


plot_list <- list()

for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
    # Read the synergy level of the  drug combination
    drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
    drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level")]
    
    
    # Read the DDI data
    DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
    DrugBank_ddi <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", paste0("ADR_", disease))]
    
    
    # Merge the data
    drugCombs_cat <- merge(DrugBank_ddi, drugCombs_data, 
                           by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"),
                           all.y = TRUE)
    

    # Assign the drug combination categories
    drugCombs_cat$class_EffAdv <- NA
    
    drugCombs_cat[drugCombs_cat$Syn_level >= 3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown",]$class_EffAdv <- "Syn_Eff"
    drugCombs_cat[drugCombs_cat$Syn_level >= 3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive",]$class_EffAdv  <- "Syn_Adv"
    
    drugCombs_cat[drugCombs_cat$Syn_level <= -3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown",]$class_EffAdv <- "Ant_Eff"
    drugCombs_cat[drugCombs_cat$Syn_level <= -3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive",]$class_EffAdv  <- "Ant_Adv"
    
    
    drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv),]
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    
    
    
    fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
    fgsea_result <- fgsea_result$combinedEfficacySafety
    fgsea_result <- fgsea_result[, colnames(fgsea_result) %in% drugCombs_cat$comb_name]
    
    plot_data <- fgsea_result
    
    
    
    feature_vars <- apply(plot_data, 1, var)
    feature_vars <- sort(feature_vars, decreasing = TRUE)
    efficacy_feature_select <- names(feature_vars[grep("^\\[DISEASE\\]", names(feature_vars))][1])
    safety_feature_select <- names(feature_vars[grep("^\\[ADR\\]", names(feature_vars))][2])
    
    
    plot_data <- plot_data[row.names(plot_data) %in% c(efficacy_feature_select, safety_feature_select), ,  drop = FALSE]
    plot_data <- as.data.frame(t(plot_data))
    colnames(plot_data) <- c("F1", "F2")
    
    plot_data <- merge(plot_data, drugCombs_cat, by.x = 0, by.y = "comb_name")
    
    plot_list[[drug_target_type]][[disease]] <- ggplot() +
      geom_point(data = plot_data, 
                 mapping = aes(x = F1, y = F2, color = class_EffAdv),  
                 size = 5,
                 shape = 3) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            text = element_text(size = 4),
            plot.title = element_text(hjust = 0.5, size = 4),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.position = "bottom",
            legend.key = element_blank(),
            legend.key.size = unit(0.5, 'cm'),
            legend.title = element_text(size = 3),
            legend.text = element_text(size = 3),
            legend.margin = margin(1,1,1,1),
            legend.box.spacing = unit(0.1, 'cm'),
            legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
      labs(title = paste0("Disease: ", disease,
                          "\nTarget: ", drug_target_type),
           x = str_wrap(efficacy_feature_select, 30),
           y = str_wrap(safety_feature_select, 30),
           color = "Class")
  }
}


