set.seed(5081)



# Script to plot the predicted probabilities on the validation data 2


# Load libraries
library(tidyverse)
library(geomtextpath)
library(ggpubr)

# Define the parameters to use
drug_target_type <- "known"


#####


# Read the prediction probabilities 
predict_result <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_2/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_result[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_2/Predictions/predictions_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))
  }
}
predict_result <- bind_rows(predict_result, .id = "disease")



#####


# Plot the data
tmp1 <- predict_result[, c("disease", "comb_name", "DDI_type",
                           "predicted_probAdv", "predicted_probEff", "predicted_category")]

tmp2 <- predict_result[, c("disease", "comb_name", "DDI_type",
                           "predicted_probAdv", "predicted_probEff", "predicted_category")]
tmp2$DDI_type <- "all"

plot_data <- rbind(tmp1, tmp2)

plot_list <- list()


# Histogram
plot_list[["histogram"]]  <- ggplot(plot_data, aes(x = predicted_probEff)) +
  geom_histogram(aes(fill = disease), alpha = 0.6) +
  facet_wrap(.~DDI_type, 
             scales = "free_y", 
             nrow = 1, 
             labeller = label_wrap_gen(width = 30, 
                                       multi_line = TRUE)) +
  labs(x = "Predicted probability (Effective)", 
       y = "Count", 
       fill = "Disease") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        # plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4), 
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, "cm"), 
        legend.key = element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))


# Density
plot_list[["density"]] <- ggplot(plot_data, aes(x = predicted_probEff, color = disease, label = disease)) +
  geom_textdensity(size = 1, linewidth = 0.1) +
  facet_wrap(.~DDI_type, 
             scales = "free_y", 
             nrow = 1, 
             labeller = label_wrap_gen(width = 30, 
                                       multi_line = TRUE)) +
  labs(x = "Predicted probability (Effective)", 
       y = "Density", 
       color = "Disease") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        # plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4), 
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, "cm"), 
        legend.key = element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))


# Save
if(!dir.exists("OutputFiles/Plots_publication/Validation/")){
  dir.create("OutputFiles/Plots_publication/Validation/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots_publication/Validation/Validation2_predictionProbabilities_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
     width = 15,
     height = 8,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list, nrow = 2, common.legend = TRUE, legend = "bottom") 

dev.off()



print(warnings())