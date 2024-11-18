set.seed(5081)


# Script to plot/tabulate the number of drug combinations based on the categories and their interaction


# Load libraries
library(tidyverse)
library(ggpubr)


#####


# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")


plot_list <- list()

# Generate the plots 

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  DrugBank_ddi_select <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "ADR_status")]
  
  drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
  drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level", "class_synergyScore")]
  
  drugCombs_cat <- drugCombs_data %>% 
    left_join(DrugBank_ddi_select, by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id")) %>%
    as.data.frame()
  
  # Replace the missing categories as unknown
  drugCombs_cat <- drugCombs_cat %>% 
    mutate_at(c("ADR_status"), 
              ~replace_na(., "unknown"))
  
  
  row.names(drugCombs_cat) <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("Syn_level", "class_synergyScore", "ADR_status")]
  
  category_combinations <- combn(x = colnames(drugCombs_cat), m = 2, simplify = FALSE)
  
  for(i in (category_combinations)){
    plot_data <- as.data.frame(table(drugCombs_cat[, i[1]], drugCombs_cat[,i[2]], useNA = "ifany"))
    
    plot_list[[disease]][[ paste(i, collapse = "__") ]] <- ggplot(plot_data, aes(x = Var1, y = Var2, )) + 
      geom_tile(aes(fill = Freq)) +
      geom_text(aes(label = Freq), size = 1) + 
      labs(title = paste0(disease, " (", nrow(drugCombs_cat), ")"),
           x = i[1],
           y = i[2]) +
      scale_fill_gradient(low = "white", high = "lightblue") +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            text = element_text(size = 4), 
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.position = "none"
      )
  }
  
}


#####


# Plot for all the cancers in single plot
if(!dir.exists("OutputFiles/Plots/drugComb_category_distribution/")){
  dir.create("OutputFiles/Plots/drugComb_category_distribution/", recursive = TRUE)
}
tiff("OutputFiles/Plots/drugComb_category_distribution/drugComb_category_distribution_xCancer.tiff",
     width = 12, height = 15,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = unlist(plot_list, recursive = FALSE), ncol = 3, nrow = 6)

dev.off()


# Plot separately for each cancer
for(disease in names(plot_list)){
  tiff(paste0("OutputFiles/Plots/drugComb_category_distribution/drugComb_category_distribution_", disease, ".tiff"),
       width = 10, height = 3,
       units = "cm", compression = "lzw", res = 1200)
  print(ggarrange(plotlist = plot_list[[disease]], ncol = 3))
  dev.off()
}


#####


print(warnings())