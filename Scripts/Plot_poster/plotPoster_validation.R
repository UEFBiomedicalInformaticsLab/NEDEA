library(tidyverse)

drug_target_type <- "known"
plot_data <- list()

#####


comb_count_1 <- read.csv("OutputFiles/Tables/Summary_drugComb_count.csv", header = TRUE, check.names = FALSE)
comb_count_1 <- comb_count_1[comb_count_1$category_type == "class_EffAdv", c("disease", "category_type", "Eff", "Adv")]
comb_count_1$Dataset <- "Training"
comb_count_1$label <- paste(comb_count_1$Eff, comb_count_1$Adv, sep = "/")
comb_count_1 <- comb_count_1[, c("Dataset", "disease", "label")]



comb_count_2 <- read.csv("OutputFiles/Tables/Summary_validNdenovo_drugComb_count.csv", 
                         header = TRUE, 
                         check.names = FALSE, row.names = c(1))

comb_count_2 <- rownames_to_column(comb_count_2, "Dataset")



comb_count_2 <- comb_count_2 %>% pivot_longer(-Dataset, 
                                              names_to = "disease", 
                                              values_to = "count") %>% 
  separate(col = "disease", into = c("disease", "comb_type"), sep = "_") %>%
  filter(comb_type %in% c("Eff", "Adv")) %>%
  filter(!Dataset %in% c("DeNovo1", "Validation4")) 

comb_count_2 <- comb_count_2 %>% pivot_wider(names_from = "comb_type", values_from = "count")
comb_count_2$label <- paste(comb_count_2$Eff, comb_count_2$Adv, sep = "/")  
comb_count_2 <- comb_count_2[, c("Dataset", "disease", "label")]

comb_count_2$Dataset <- gsub("Validation", "Validation dataset ", comb_count_2$Dataset)

comb_count <- rbind(comb_count_1, comb_count_2)

#####


# Training

final_model_accuracy<- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  
  final_model_accuracy[[disease]] <- model$final_model_accuracy
}

final_model_accuracy <- bind_rows(final_model_accuracy, .id = "disease")
final_model_accuracy <- final_model_accuracy[final_model_accuracy$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
plot_data[["Training"]] <- final_model_accuracy


#####


# Validation data 1
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_1/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_1/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
  }
}
predict_metrics <- bind_rows(predict_metrics, .id = "disease")
predict_metrics <- predict_metrics[predict_metrics$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
plot_data[["Validation dataset 1"]] <- predict_metrics


#####


# Validation dataset 2

# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  if(file.exists(paste0("OutputFiles/Validation_data_2/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_2/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
  }
}
predict_metrics <- bind_rows(predict_metrics, .id = "disease")

# Select the metrices to plot
predict_metrics <- predict_metrics[predict_metrics$DDI_type == "all", c("disease", "Specificity")]
predict_metrics <- pivot_longer(predict_metrics, cols = "Specificity", names_to = ".metric", values_to = ".estimate")
predict_metrics$.metric <- "specificity"
plot_data[["Validation dataset 2"]] <- predict_metrics


#####


# Read the all the metrics
predict_metrics <- list()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  if(file.exists(paste0("OutputFiles/Validation_data_3/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"))){
    predict_metrics[[disease]] <- read.csv(file = paste0("OutputFiles/Validation_data_3/Prediction_metrics/predictionMetrics_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".csv"), header = TRUE)
  }
}
predict_metrics <- bind_rows(predict_metrics, .id = "disease")
predict_metrics <- predict_metrics[predict_metrics$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
plot_data[["Validation dataset 3"]] <- predict_metrics


#####

plot_data <- bind_rows(plot_data, .id = "Dataset")

# plot_data$.metric <- gsub("sensitivity", "Sensitivity", plot_data$.metric)
# plot_data$.metric <- gsub("specificity", "Specificity", plot_data$.metric)
# plot_data$.metric <- gsub("bal_accuracy", "Balanced Accuracy", plot_data$.metric)

plot_data$.metric <- factor(plot_data$.metric, 
                            levels = c("bal_accuracy", "sensitivity", "specificity"),
                            labels = c("Balanced Accuracy", "Sensitivity", "Specificity"))


plot_data <- plot_data %>% left_join(comb_count, by = c("Dataset", "disease"))


if(!dir.exists("OutputFiles/Plots_poster/Validation/")){
  dir.create("OutputFiles/Plots_poster/Validation/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots_poster/Validation/Validation_predictionMetrics_NES_combinedEfficacySafety_", drug_target_type, ".tiff"),
     width = 39,
     height = 11,
     units = "cm", compression = "lzw", res = 1200)

ggplot(plot_data) +
  geom_bar(aes(x = disease, y = .estimate, fill = .metric),
           stat = "identity",
           position =  position_dodge(preserve = "single"),
           width = 0.5,
           lwd = 0.1, na.rm = TRUE) +
  geom_hline(yintercept = 0.7,
             linetype="dotted",
             color = "#006400",
             linewidth = 0.1) +
  geom_text(aes(x = disease, y = 1, label = label), size = 3, na.rm = TRUE) +
  facet_wrap(~ Dataset, ncol = 4) +
  labs(x = "Disease", y = "Score", fill = "Metric: ") +
  scale_x_discrete(labels = function(x){ gsub("Cancer$", " Cancer", x)  } ) +
  theme(plot.margin = margin(l = 0.1, r = 0.1, t = 0.1, b = 0.1, unit = "cm"),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 14, margin = margin(2,1,2,1)),
        text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "plain"),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 5)),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.key.spacing.x = unit(0.25, "cm"),
        legend.text = element_text(size = 12, margin = margin(l = 1)),
        legend.margin = margin(2,1,2,1),
        legend.box.spacing = unit(0.25, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))

dev.off()


print(warnings())



