set.seed(5081)


# Script to plot the accuracy scores from validation 2


# Load libraries
library(tidyverse)


# Set the parameters
drug_target_type <- "known"


# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_2/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_2/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
  }
}
predict_metrics <- bind_rows(predict_metrics, .id = "disease")

predict_metrics <- predict_metrics[order(predict_metrics$Number_of_combinations, decreasing = TRUE), ]
# predict_metrics$DDI_type <- as.factor(predict_metrics$DDI_type)


# Select the metrices to plot
plot_data <- predict_metrics


#Plot
if(!dir.exists("OutputFiles/Plots_publication/Validation/")){
  dir.create("OutputFiles/Plots_publication/Validation/", recursive = TRUE)
}

# tiff(paste0("OutputFiles/Plots_publication/Validation/Validation2_predictionMetrics_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
#      width = 15,
#      height = 8,
#      units = "cm", compression = "lzw", res = 1200)
# 
# ggplot(plot_data, aes(x = DDI_type, y = Specificity)) +
#   geom_bar(stat = "identity",
#            position = "dodge",
#            width = 0.5,
#            lwd = 0.1) +
#   geom_text(aes(label = Number_of_combinations), 
#             position = position_dodge(width = 0.5), 
#             vjust = -0.25,
#             hjust = 0.5,
#             angle = 0, 
#             size = 1) + 
#   geom_hline(yintercept = 0.7,
#              linetype="dashed",
#              color = "green",
#              linewidth = 0.1) +
#   geom_hline(yintercept = 0.5,
#              linetype="dashed",
#              color = "red",
#              linewidth = 0.1) +
#   facet_wrap(.~disease, 
#              nrow = 2, 
#              scales = "free_x") +
#   scale_x_discrete(labels = label_wrap_gen(width = 35))  +
#   labs(x = "Disease", y = "Specificity") +
#   theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
#         panel.grid = element_blank(),
#         panel.spacing = unit(0.1, "cm"),
#         strip.background = element_rect(color = "black", linewidth = 0.25,),
#         strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
#         text = element_text(size = 4),
#         plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
#         axis.title = element_text(size = 4),
#         axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1),
#         axis.ticks = element_line(colour = "black", linewidth = 0.2))


if(!dir.exists("OutputFiles/Plots_publication/Validation/")){
  dir.create("OutputFiles/Plots_publication/Validation/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots_publication/Validation/Validation2_predictionMetrics_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
     width = 8,
     height = 6,
     units = "cm", compression = "lzw", res = 1200)

plot_data$tile_label <- paste0(round(plot_data$Specificity, 2), "\n(", plot_data$Number_of_combinations, ")")
ggplot(plot_data, aes(x = disease, y = DDI_type, fill = Specificity, label = tile_label)) +
  geom_tile() +
  geom_text(size = 1) + 
  scale_y_discrete(labels = label_wrap_gen(width = 30))  +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(x = "Disease", y = "DDI type") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4),
        plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 4),
        axis.text.x = element_text(size = 2.5, angle = 0, vjust = 0, hjust = 0.5),
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


print(warnings())