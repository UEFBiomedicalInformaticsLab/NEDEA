set.seed(5081)


# Script to plot the accuracy scores from validation (for poster)


# Load libraries
library(tidyverse)
library(ggradar)

# Set the parameters
drug_target_type <- "known"


plot_list <- list()


#####


# Training accuracy

final_model_accuracy <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  
  final_model_accuracy[[disease]] <- model$final_model_accuracy
}

final_model_accuracy <- bind_rows(final_model_accuracy, .id = "disease")


final_model_accuracy <- final_model_accuracy[final_model_accuracy$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
final_model_accuracy$.metric <- factor(final_model_accuracy$.metric, 
                                       levels = c("bal_accuracy", "sensitivity", "specificity"),
                                       labels = c("Balanced Accuracy", "Sensitivity", "Specificity"))



plot_data <- pivot_wider(data = final_model_accuracy, id_cols = "disease", names_from = ".metric", values_from = ".estimate", )


plot_list[["Train"]] <- ggradar(plot.data = plot_data)


#####


# Validation data 1


# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_1/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_1/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
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

plot_data <- pivot_wider(data = plot_data, id_cols = "disease", names_from = ".metric", values_from = ".estimate", )


plot_list[["Val1"]] <-  ggradar(plot.data = plot_data)


#####


# Validation data 2


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
plot_data <- predict_metrics[predict_metrics$DDI_type == "all", c("disease", "Specificity")]
# plot_data$`Balanced Accuracy` <- 0
# plot_data$Sensitivity <- 0


#Plot


plot_list[["Val2"]] <- ggradar(plot.data = plot_data)


#####


# Validation data 3


# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  if(file.exists(paste0("OutputFiles/Validation_data_3/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_3/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
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


plot_data[is.na(plot_data$.estimate), ".estimate"] <- 0

plot_data <- pivot_wider(data = plot_data, id_cols = "disease", names_from = ".metric", values_from = ".estimate", )

#Plot

plot_list[["Val3"]] <-  ggradar(plot.data = plot_data)


#####


# Arrange the plots


if(!dir.exists("OutputFiles/Plots_poster/Validation/")){
  dir.create("OutputFiles/Plots_poster/Validation/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots_poster/Validation/ValidationRadar_predictionMetrics_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
     width = 39,
     height = 15,
     units = "cm", compression = "lzw", res = 1200)

library(ggpubr)

ggpubr::ggarrange(plotlist = plot_list, common.legend = TRUE, ncol = 4, nrow = 1, legend = "bottom")

dev.off()
