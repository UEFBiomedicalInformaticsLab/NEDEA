set.seed(5081)


# Load libraries
library(tidyverse)


# Script to check the drug class in de novo 1 dataset


# Read the ATC codes of the drugs from Drug Bank
DrugBank_drug_ATC <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_drug_ATC$drugs$atc_codes
colnames(DrugBank_drug_ATC) <- gsub("drugbank-id", "DrugBank_drug_ID", colnames(DrugBank_drug_ATC))


if(!dir.exists("OutputFiles/Plots/ATC_classification/DeNovo_data_1/")){
  dir.create("OutputFiles/Plots/ATC_classification/DeNovo_data_1/", recursive = TRUE)
}

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the drug combination 
  drugCombs <- readRDS(paste0("InputFiles/DeNovo_data_1/drugCombs_denovo1_", disease, ".rds"))
  drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
  
  # Add the ATC codes at level 1
  # Using many-to-many mapping to map all possible ATC codes to a single drug
  # If the level 1
  drugCombs <- drugCombs %>% 
    left_join(DrugBank_drug_ATC %>% 
                select(level_1, DrugBank_drug_ID) %>%  
                rename_with(.cols = everything(), 
                            .fn = ~ paste0("Drug1_ATC_", .)), 
              by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"), 
              relationship = "many-to-many") %>%
    left_join(DrugBank_drug_ATC %>% 
                select(level_1, DrugBank_drug_ID) %>%  
                rename_with(.cols = everything(), 
                            .fn = ~ paste0("Drug2_ATC_", .)), 
              by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"), 
              relationship = "many-to-many") %>% 
    distinct()
  
  drugCombs[is.na(drugCombs$Drug1_ATC_level_1), "Drug1_ATC_level_1"] <- "MISSING"
  drugCombs[is.na(drugCombs$Drug2_ATC_level_1), "Drug2_ATC_level_1"] <- "MISSING"
  
  
  # Get all the listed ATCs
  all_atc <- sort(unique(c(drugCombs$Drug1_ATC_level_1, drugCombs$Drug2_ATC_level_1)))
  
  
  ATC_count_list <- list()
  
  
  # Count the ATC occurance for the combinations
  ATC_count_mat <- matrix(0, 
                          nrow = length(all_atc), 
                          ncol = length(all_atc), 
                          dimnames = list(all_atc, all_atc)
  )
  ATC_count_mat[upper.tri(ATC_count_mat, diag = FALSE)] <- NA
  
  if(nrow(drugCombs) > 0){
    for(i in 1:nrow(drugCombs)){
      
      atc_1 <- drugCombs[i, "Drug1_ATC_level_1", drop = TRUE]
      atc_2 <- drugCombs[i, "Drug2_ATC_level_1", drop = TRUE]
      
      
      if(!is.na(ATC_count_mat[atc_1, atc_2])){
        ATC_count_mat[atc_1, atc_2] <- ATC_count_mat[atc_1, atc_2] + 1
      }else{
        ATC_count_mat[atc_2, atc_1] <- ATC_count_mat[atc_2, atc_1] + 1
        
      }
    }
  }
  
  ATC_count_mat <- as.data.frame(ATC_count_mat) %>% 
    rownames_to_column("ATC1") %>% 
    pivot_longer(-ATC1, 
                 names_to = "ATC2", 
                 values_to = "Count")
  
  
  # Plot the counts of pairs as heatmap
  plot_data <- ATC_count_mat  %>% filter(!(Count == 0 | is.na(Count)))
  plot_data$ATC1 <- factor(x = plot_data$ATC1, levels = sort(unique(plot_data$ATC1), decreasing = FALSE))
  plot_data$ATC2 <- factor(x = plot_data$ATC2, levels = sort(unique(plot_data$ATC2), decreasing = FALSE))
  
  plot <- ggplot(plot_data, aes(x = ATC1, y = ATC2, fill = Count, label = Count)) +
    geom_tile() +
    geom_text(size = 1, na.rm = TRUE) + 
    labs(title = gsub("Cancer$", " Cancer", disease)) +
    scale_fill_distiller(palette = "YlGn", direction = 1, na.value = "white") +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_line(color = "gray", linewidth = 0.05),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
          text = element_text(size = 4),
          plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 4),
          axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 2.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "right",
          legend.key = element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.title = element_text(size = 2, face = "bold", margin = margin(0.5,1,2,1)),
          legend.text = element_text(size = 2),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25))
  
  
  tiff(paste0("OutputFiles/Plots/ATC_classification/DeNovo_data_1/drugCombs_ATCclass_", disease, ".tiff"),
       width = 10,
       height = 8,
       units = "cm", compression = "lzw", res = 1200)
  
  print(plot)
  
  dev.off()
  
}


print(warnings())