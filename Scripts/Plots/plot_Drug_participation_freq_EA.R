set.seed(5081)



# Script to plot the frequency of participation of the independent drugs in drug combinations (only EA category)


# Load libraries
library(tidyverse)


file_width <- 0
plot_list <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the drug combination categories
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "comb_name", "class_EffAdv")]
  
  # Check the overlap of drugs participating in each category
  drugs_in_combs <- sort(unique(c(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id)))
  file_width <- ifelse(length(drugs_in_combs) > file_width, length(drugs_in_combs), file_width)
  
  drug_participation_mat <- matrix(data = NA, 
                                   nrow = length(drugs_in_combs), 
                                   ncol = length(drugs_in_combs),
                                   dimnames = list(drugs_in_combs, drugs_in_combs)
  )
  
  for(i in 1:nrow(drugCombs_cat)){
    drug_participation_mat[drugCombs_cat[i,]$Drug1_DrugBank_id, drugCombs_cat[i,]$Drug2_DrugBank_id] <- drugCombs_cat[i,]$class_EffAdv
    drug_participation_mat[drugCombs_cat[i,]$Drug2_DrugBank_id, drugCombs_cat[i,]$Drug1_DrugBank_id] <- drugCombs_cat[i,]$class_EffAdv
  }
  
  drug_participation_mat <- apply(drug_participation_mat, 2, function(x){ as.data.frame(table(x)) })
  drug_participation_mat <- bind_rows(drug_participation_mat, .id = "Drugs")
  tmp1 <- drug_participation_mat %>% group_by(Drugs) %>% summarise("Count" = sum(Freq)) %>% arrange(desc(Count))
  drug_participation_mat$Drugs <- factor(drug_participation_mat$Drugs, levels = tmp1$Drugs) 
  drug_participation_mat$x <- factor(drug_participation_mat$x, levels = c("Adv", "Eff")) 
  
  plot_list[[disease]] <- ggplot(drug_participation_mat, 
                                 aes(x = Drugs, y = Freq, fill = x)) +
    geom_bar(stat = "identity") + 
    labs(title = gsub("Cancer$", " Cancer", disease),
         fill = "Drug combination type:") +
    scale_fill_manual(values = c("Adv" = "#FF6961", "Eff" = "#77DD77"), labels = c("Adv" = "Adverse", "Eff" = "Effective")) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          # strip.background = element_rect(color = "black", linewidth = 0.25,),
          # strip.text = element_text(size = 6, margin = margin(1,1,1,1)),
          text = element_text(size = 5), 
          plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.title = element_text(margin = margin(r = 2)),
          legend.key = element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.key.spacing.x = unit(0.1, "cm"),
          legend.text = element_text(size = 5, margin = margin(l = 1)),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
}



tiff("OutputFiles/Plots/Drug_participation_freq_EA.tiff",
     width = file_width/4 + 1, height = 20,
     units = "cm", compression = "lzw", res = 1200)

ggpubr::ggarrange(plotlist = plot_list, nrow = 6, ncol = 1, common.legend = TRUE, legend = "bottom")

dev.off()




print(warnings())