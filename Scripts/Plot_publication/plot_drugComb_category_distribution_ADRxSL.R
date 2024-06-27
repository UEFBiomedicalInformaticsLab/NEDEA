set.seed(5081)


# Script to plot/tabulate the number of drug combinations based on the categories and their interaction (ADR vs SL)


# Load libraries
library(tidyverse)
library(ggpubr)



# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")


plot_list <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  DrugBank_ddi_select <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_therapeuticEfficacy", "class_metabolicEffect", paste0("ADR_", disease))]
  
  drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
  drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level")]
  
  drugCombs_cat <- merge(DrugBank_ddi_select, drugCombs_data, 
                         by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"),
                         all.y = TRUE)
  
  
  # Calculate the number the combinations with positive ADR
  print(paste0(disease, " :: ADR_positive = ", 
               round(( nrow(drugCombs_cat[drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive", ]) / 
                         nrow(drugCombs_cat) ) * 100, 2), 
               "; unknown = ", 
               round(( nrow(drugCombs_cat[drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown", ]) / 
                         nrow(drugCombs_cat) ) * 100, 2)
  ))
  

  # Remove drug combinations with ADR as NA
  # These were the drug combinations for which absolutely no report present
  drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat[, paste0("ADR_", disease)]), ]
  
  
  row.names(drugCombs_cat) <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("Syn_level", paste0("ADR_", disease))]
  
  category_combinations <- combn(x = colnames(drugCombs_cat), m = 2, simplify = FALSE)
  
  

  
  for(i in (category_combinations)){
    plot_data <- as.data.frame(table(drugCombs_cat[, i[1]], drugCombs_cat[,i[2]], useNA = "ifany"))
    
    plot_list[[disease]][[ paste(i, collapse = "__") ]] <- ggplot(plot_data, aes(x = Var1, y = Var2, )) + 
      geom_tile(aes(fill = Freq)) +
      geom_text(aes(label = Freq), size = 0.75) + 
      labs(title = gsub("Cancer$", " Cancer", disease),
           x = "Synergy level",
           y = "ADR status") +
      scale_y_discrete(labels = c("adr_positive" = "ADR positive", "unknown" = "Unknown")) +
      scale_fill_gradient(low = "white", high = "lightblue") +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            text = element_text(size = 2), 
            plot.title = element_text(size = 2.5, hjust = 0.5, face = "bold", margin = margin(2, 1, 1, 1, "pt")),
            plot.margin = margin(1, 2, 1, 2, "pt"),
            axis.title = element_text(size = 2.5), 
            axis.title.y = element_text(margin = margin(r = 0)),
            axis.text.x = element_text(size = 2.5, angle = 0, vjust = 1, hjust = 0.5), 
            axis.text.y = element_text(size = 2.5, hjust = 0.5, margin = margin(l = 1)),
            axis.ticks = element_line(colour = "black", linewidth = 0.1),
            legend.position = "none"
      )
  }
  
}


if(!dir.exists("OutputFiles/Plots_publication/")){
  dir.create("OutputFiles/Plots_publication/", recursive = TRUE)
}
tiff("OutputFiles/Plots_publication/drugComb_category_distribution_xCancer_ADRxSL.tiff",
     width = 8, height = 4,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = unlist(plot_list, recursive = FALSE), ncol = 3, nrow = 2)


dev.off()


print(warnings())