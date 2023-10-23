set.seed(5081)



# Script to plot the NES of top (by variance) features from FGSEA of  SMPDB (drug metabolism pathways) library


# Load libraries
library(tidyverse)
library(ggpubr)


plot_list <- list()


for(drug_target_type in c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("comb_name", "class_EffAdv")]
    
    
    fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_Pathway_", disease, "_", drug_target_type, ".rds"))
    fgsea_result <- fgsea_result$smpdbDrugMet
    fgsea_result <- fgsea_result[, colnames(fgsea_result) %in% drugCombs_cat$comb_name]
    
    plot_data <- fgsea_result
    
    feature_vars <- apply(plot_data, 1, var)
    feature_vars <- sort(feature_vars, decreasing = TRUE)
    feature_vars <- feature_vars[1:2]
    
    plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), ,  drop = FALSE]
    plot_data <- as.data.frame(t(plot_data))
    colnames(plot_data) <- c("F1", "F2")
    
    plot_data <- merge(plot_data, drugCombs_cat, by.x = 0, by.y = "comb_name")
    
    plot_list[[drug_target_type]][[disease]] <- ggplot() +
      geom_point(data = plot_data, 
                 mapping = aes(x = F1, y = F2, color = class_EffAdv),  
                 size = 0.5,
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
           x = str_wrap(names(feature_vars)[1], 30),
           y = str_wrap(names(feature_vars)[2], 30),
           color = "Class")
  }
}


if(!dir.exists("OutputFiles/Plots/top_feature_NES_scatter/")){
  dir.create("OutputFiles/Plots/top_feature_NES_scatter/", recursive = TRUE)
}


tiff("OutputFiles/Plots/top_feature_NES_scatter/plot_NES_smpdbDrugMet_topFeature_byVar.tiff",
     width = 30, height = 25,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = unlist(plot_list, recursive = FALSE), 
          ncol = 6, nrow = 7, 
          common.legend = TRUE, legend = "bottom")

print(warnings())