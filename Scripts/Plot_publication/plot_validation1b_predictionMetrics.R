set.seed(5081)


# Script to plot the accuracy scores from validation 1b


# Load libraries
library(tidyverse)


# Set the parameters
drug_target_type <- "known"


# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_1b/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_1b/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
  }
}
predict_metrics <- bind_rows(predict_metrics, .id = "disease")


# Select the metrices to plot
plot_data <- predict_metrics[predict_metrics$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]


# plot_data$.metric <- gsub("f_meas", "F1 score", plot_data$.metric)
plot_data$.metric <- gsub("sensitivity", "Sensitivity", plot_data$.metric)
plot_data$.metric <- gsub("specificity", "Specificity", plot_data$.metric)
plot_data$.metric <- gsub("bal_accuracy", "Balanced Accuracy", plot_data$.metric)
# plot_data$.metric <- gsub("roc_auc", "ROC-AUC", plot_data$.metric)
# plot_data$.metric <- gsub("pr_auc", "PR-AUC", plot_data$.metric)


#Plot
if(!dir.exists("OutputFiles/Plots_publication/Validation/")){
  dir.create("OutputFiles/Plots_publication/Validation/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots_publication/Validation/Validation1b_predictionMetrics_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
     width = 8,
     height = 3,
     units = "cm", compression = "lzw", res = 1200)

ggplot(plot_data, aes(x = disease, y = .estimate, label = round(.estimate, 2))) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.5, 
           lwd = 0.1) +
  geom_text(size = 1) +
  geom_hline(yintercept = 0.7,
             linetype = "dotted",
             color = "#006400",
             linewidth = 0.1) +
  facet_wrap(~ .metric, ncol = 3) +
  labs(x = "Disease", y = "Score") + 
  scale_x_discrete(labels = function(x){ gsub("Cancer$", " Cancer", x)  } ) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4), 
        axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()


print(warnings())