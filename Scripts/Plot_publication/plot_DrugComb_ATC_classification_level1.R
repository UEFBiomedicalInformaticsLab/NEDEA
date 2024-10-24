set.seed(5081)


# Script to plot the level_1 ATC codes as heatmap

# Load libraries
library(tidyverse)


#####


# Read the ATC codes of the drugs from Drug Bank
DrugBank_drug_ATC <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_drug_ATC$drugs$atc_codes
colnames(DrugBank_drug_ATC) <- gsub("drugbank-id", "DrugBank_drug_ID", colnames(DrugBank_drug_ATC))


#####


# Read the drug combinations from all the datasets 

drugCombs_master <- list()

for(dataset in c("Training", "Validation1", "Validation2", "Validation3")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    switch(dataset, 
           "Training" = {
             # Read the drug combination categories
             drugCombs <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
             drugCombs$comb_name <- paste(drugCombs$Drug1_DrugBank_id, drugCombs$Drug2_DrugBank_id, sep = "_")
             drugCombs <- drugCombs[!is.na(drugCombs$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             drugCombs_master[[dataset]][[disease]] <- drugCombs
           },
           
           "Validation1" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             drugCombs_master[[dataset]][[disease]] <- drugCombs
           },
           
           "Validation2" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             drugCombs_master[[dataset]][[disease]] <- drugCombs
           },
           
           "Validation3" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             drugCombs_master[[dataset]][[disease]] <- drugCombs
           }
    )
    
  }
}

drugCombs <- drugCombs_master
rm(drugCombs_master)

drugCombs <- unlist(drugCombs, recursive = FALSE)

drugCombs <- bind_rows(drugCombs, .id = "source") %>%
  separate(col = "source",
           into = c("Dataset", "Disease"),
           sep = "\\.")
drugCombs$comb_name <- paste(drugCombs$Drug1_DrugBank_id, drugCombs$Drug2_DrugBank_id, sep = "_")


#####


# Add the ATC codes at level_1
# Using many-to-many mapping to map all possible ATC codes to a single drug
drugCombs <- drugCombs %>%
  left_join(DrugBank_drug_ATC %>%
              select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug1_ATC_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  left_join(DrugBank_drug_ATC %>%
              select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug2_ATC_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  distinct()


drugCombs[is.na(drugCombs$Drug1_ATC_code_1), "Drug1_ATC_code_1"] <- "MISSING"
drugCombs[is.na(drugCombs$Drug2_ATC_code_1), "Drug2_ATC_code_1"] <- "MISSING"


#####


# Get all the listed ATCs
all_atc <- sort(unique(c(drugCombs$Drug1_ATC_code_1, drugCombs$Drug2_ATC_code_1)))


ATC_count_list <- list()

# Count the ATC occurance for the combinations

for(dataset in c("Training", "Validation1", "Validation2", "Validation3")){
  
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    ATC_count_mat <- matrix(0, 
                            nrow = length(all_atc), 
                            ncol = length(all_atc), 
                            dimnames = list(all_atc, all_atc)
    )
    ATC_count_mat[upper.tri(ATC_count_mat, diag = FALSE)] <- NA
    
    drugCombs_select <- drugCombs[drugCombs$Dataset %in% dataset & drugCombs$Disease %in% disease, ]
    
    
    if(nrow(drugCombs_select) > 1){
      for(i in 1:nrow(drugCombs_select)){
        
        atc_1 <- drugCombs_select[i, "Drug1_ATC_code_1"]
        atc_2 <- drugCombs_select[i, "Drug2_ATC_code_1"]
        
        
        if(!is.na(ATC_count_mat[atc_1, atc_2])){
          ATC_count_mat[atc_1, atc_2] <- ATC_count_mat[atc_1, atc_2] + 1
        }else{
          ATC_count_mat[atc_2, atc_1] <- ATC_count_mat[atc_2, atc_1] + 1
          
        }
      }
    }
    
    tmp1 <- as.data.frame(ATC_count_mat) %>% 
      rownames_to_column("ATC1") %>% 
      pivot_longer(-ATC1, 
                   names_to = "ATC2", 
                   values_to = "Count")
    tmp1$Disease <- disease
    tmp1$Dataset <- dataset
    
    ATC_count_list[[dataset]][[disease]] <- tmp1
    rm(tmp1)
  }}



#####


# Plot
plot_data <- unlist(ATC_count_list, recursive = FALSE)
plot_data <- bind_rows(plot_data) 

plot_data$log_count <- log(plot_data$Count, 10)
plot_data[is.infinite(plot_data$log_count), "log_count"] <- NA

tiff(paste0("OutputFiles/Plots_publication/drugCombs_ATCclass_level1.tiff"),
     width = 29,
     height = 20,
     units = "cm", compression = "lzw", res = 1200)


ggplot(plot_data, aes(x = ATC1, y = ATC2, fill = log_count, label = Count)) +
  geom_tile() +
  # geom_text(size = 1) + 
  labs(fill = "log(Count, 1)") + 
  facet_grid(Dataset ~ Disease, 
             labeller = labeller(.cols = function(x){ gsub("Cancer$", " Cancer", x) }) ) +
  scale_fill_continuous(low = "blue", high = "red", na.value = "white") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_line(color = "gray", linewidth = 0.05),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4),
        plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4),
        axis.text.x = element_text(size = 1, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.title = element_text(size = 2, face = "bold", margin = margin(0.5,1,2,1)),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))

dev.off()


#####


print(warnings())