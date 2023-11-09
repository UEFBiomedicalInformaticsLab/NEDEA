set.seed(5081)


# Script to plot the number of drug targets and the number of drug combinations per disease


# Load libraries 
library(tidyverse)




# Read the drug targets
drugCombs_targets_list <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
  drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
  drugCombs_targets <- drugCombs_targets %>% select("comb_name", ends_with("count"), ends_with("_tar_cnt"))
  colnames(drugCombs_targets) <- gsub("^ext_|_tar_cnt$", "", colnames(drugCombs_targets))
  colnames(drugCombs_targets) <- gsub("drugTarget_count", "DrugBank", colnames(drugCombs_targets))
  
  drugCombs_targets_list[[disease]] <- drugCombs_targets
}

drugCombs_targets_list <- bind_rows(drugCombs_targets_list, .id = "Disease")



# Plot histogram of the number of drug targets based on different types
plot_data <- pivot_longer(data = drugCombs_targets_list, 
                          cols = colnames(drugCombs_targets_list)[-c(1:2)], 
                          names_to = "Target_source", 
                          values_to = "Target_count")

plot_stats <- plot_data %>% 
  group_by(Disease, Target_source) %>% 
  summarise(
    Min = min(Target_count),
    Max = max(Target_count),
    Mean = round(mean(Target_count), 0),
    Mode = names(which.max(table(Target_count)))
  )

plot_data <- plot_data %>% left_join(plot_stats, by = c("Disease", "Target_source"))


if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}
tiff("OutputFiles/Plots/drugCombs_target_freq.tiff",
     width = 30, height = 20,
     units = "cm", compression = "lzw", res = 1200)

ggplot(data = plot_data) +
  geom_histogram(aes(x = Target_count)) + 
  geom_text(aes(x = Inf, y = Inf, label = paste0("Min: ", Min, 
                                                 "\nMax: ", Max, 
                                                 "\nMean: ", Mean, 
                                                 "\nMode: ", Mode)
  ), 
  hjust = 2, vjust = 1.5, size = 2) +
  facet_grid(rows = vars(Target_source), cols = vars(Disease)) + 
  labs(x = "Target count", y = "No. of drug combinations") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.text = element_text(size = 8), 
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2)
  )

dev.off()





print(warnings())