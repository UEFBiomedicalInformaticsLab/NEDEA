set.seed(5081)


# Script to compare the ATC codes (level 3) of the drugs forming various datasets

# Load libraries
library(unixtools)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Read the ATC codes of the drugs from Drug Bank
DrugBank_drug_ATC <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_drug_ATC$drugs$atc_codes
colnames(DrugBank_drug_ATC) <- gsub("drugbank-id", "DrugBank_drug_ID", colnames(DrugBank_drug_ATC))


#####


# Read the drug combinations from all the datasets and extract the ATC codes

ATC_master <- list()

for(dataset in c("Training", "Validation1", "Validation1a", "Validation2", "Validation2a", "Validation3", "Validation3a")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    switch(dataset, 
           "Training" = {
             # Read the drug combination categories
             drugCombs <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
             drugCombs <- drugCombs[!is.na(drugCombs$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation1" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_1/drugCombs_validation1_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation1a" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_1a/drugCombs_validation1a_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation2" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation2a" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_2a/drugCombs_validation2a_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation3" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_3/drugCombs_validation3_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           },
           
           "Validation3a" = {
             drugCombs <- readRDS(paste0("InputFiles/Validation_data_3a/drugCombs_validation3a_", disease, ".rds"))
             drugCombs <- drugCombs[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")]
             
             drugCombs <- drugCombs %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug1_ATC_", .)),
                         by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               left_join(DrugBank_drug_ATC %>%
                           select(code_3, DrugBank_drug_ID) %>%
                           rename_with(.cols = everything(),
                                       .fn = ~ paste0("Drug2_ATC_", .)),
                         by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
                         relationship = "many-to-many") %>%
               distinct()
             
             ATC_master[[disease]][[dataset]] <- sort(na.exclude(unique(c(drugCombs$Drug1_ATC_code_3, drugCombs$Drug2_ATC_code_3))))
           }
           
    )
    
  }
}


#####


# Compare the ATC codes
res <- data.frame()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  for(dataset1 in c("Training")){
    for(dataset2 in c("Validation1", "Validation1a", "Validation2", "Validation2a", "Validation3", "Validation3a")){
      dataset1_atc <- ATC_master[[disease]][[dataset1]]
      dataset2_atc <- ATC_master[[disease]][[dataset2]]
      
      res <- rbind(res, data.frame("Disease" = disease, 
                                   "Dataset1" = dataset1, 
                                   "Dataset2" = dataset2,
                                   "Dataset1_atc_count" = length(dataset1_atc), 
                                   "Dataset2_atc_count" = length(dataset2_atc),
                                   "Dataset1_unique" = length(setdiff(dataset1_atc, dataset2_atc)), 
                                   "Common" =  length(intersect(dataset1_atc, dataset2_atc)), 
                                   "Dataset2_unique" =  length(setdiff(dataset2_atc, dataset1_atc)))
      )
    }
  }
}


#####


# Plot
plot_data <- res %>%
  pivot_longer(- c(Disease, Dataset1, Dataset2, 
                   Dataset1_atc_count, Dataset2_atc_count), 
               names_to = "count_type", 
               values_to = "count")

plot_data$Dataset2 <- factor(x = plot_data$Dataset2, levels = c("Validation1", "Validation1a", "Validation2", "Validation2a", "Validation3", "Validation3a"))
plot_data$count_type <- factor(x = plot_data$count_type, levels = c("Dataset1_unique", "Common", "Dataset2_unique"))


if(!dir.exists("OutputFiles/Plots_publication/ATC_classification/")){dir.create("OutputFiles/Plots_publication/ATC_classification/", recursive = TRUE)}

tiff(paste0("OutputFiles/Plots_publication/ATC_classification/drugCombs_ATCclass_compare_level3.tiff"),
     width = 8,
     height = 6,
     units = "cm", compression = "lzw", res = 1200)


ggplot(plot_data, aes(x = count_type, y = Dataset2, label = count))+
  geom_tile(fill = "white", color = "blue") +
  geom_text(size = 1) +
  labs(x = "", y = "") +
  facet_wrap( ~ Disease,
              labeller = labeller(.cols = function(x){ gsub("Cancer$", " Cancer", x) }) ) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_line(color = "gray", linewidth = 0.05),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4),
        plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 4),
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