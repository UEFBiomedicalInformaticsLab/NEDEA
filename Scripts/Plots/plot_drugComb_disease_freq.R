set.seed(5081)



# Script to plot the frequency of disease participation for the drug combinations


# Load libraries
library(tidyverse)
library(ggpubr)



plot_list <- list()


# Plot the disease participation with initial count
drugCombs_data_list <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
  drugCombs_data$comb_name <- paste(drugCombs_data$Drug1_DrugBank_id, drugCombs_data$Drug2_DrugBank_id, sep = "_")
  drugCombs_data <- drugCombs_data %>% select(c("comb_name", "Drug1_DrugBank_id", "Drug2_DrugBank_id"))
  drugCombs_data_list[[disease]] <- drugCombs_data
}

drugCombs_data_list <- bind_rows(drugCombs_data_list, .id = "Disease")

plot_data <- drugCombs_data_list %>% 
  select(c("Disease", "comb_name")) %>%
  unique() %>%
  group_by(comb_name) %>% 
  summarise(Count = length(Disease))

plot_list[["Before"]] <- ggplot(plot_data) +
  geom_bar(aes(x = Count)) +
  labs(title = "Before",
       x = "No. of diseases", 
       y = "No. of drug combinations") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2)
  )



# Plot the disease participation with final count
drugCombs_data_list <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_data <- drugCombs_data[!is.na(drugCombs_data$class_EffAdv), ]
  drugCombs_data$comb_name <- paste(drugCombs_data$Drug1_DrugBank_id, drugCombs_data$Drug2_DrugBank_id, sep = "_")
  drugCombs_data <- drugCombs_data %>% select(c("comb_name", "Drug1_DrugBank_id", "Drug2_DrugBank_id"))
  drugCombs_data_list[[disease]] <- drugCombs_data
}

drugCombs_data_list <- bind_rows(drugCombs_data_list, .id = "Disease")

plot_data <- drugCombs_data_list %>% 
  select(c("Disease", "comb_name")) %>%
  unique() %>%
  group_by(comb_name) %>% 
  summarise(Count = length(Disease))

plot_list[["After"]] <- ggplot(plot_data) +
  geom_bar(aes(x = Count)) +
  labs(title = "After", 
       x = "No. of diseases", 
       y = "No. of drug combinations") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2)
  )




# Save plot
if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}
tiff("OutputFiles/Plots/drugCombs_disease_freq.tiff",
     width = 8, height = 4,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list)

dev.off()



print(warnings())