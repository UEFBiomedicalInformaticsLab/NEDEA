set.seed(5081)


# Script to plot the characteristics of the predictive model like thresholds and accuracy


# Load libraries
library(tidyverse)
library(ggpubr)


# Specify the drug target type to use
drug_target_type <- "known"


#####


# Plot the threshold selected in during CV and the final threshold
plot_list <- list()
final_model_feature_threshold <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  final_model_feature_threshold[[disease]] <- model$feature_threshold
  cv_feature_threshold <- model$internal_CV_results$selected_threshold
  
  plot_list[[disease]] <- ggplot() +
    geom_boxplot(data = cv_feature_threshold, 
                 mapping = aes(x = feature, y = threshold), 
                 fill = "grey", 
                 width = 0.5, 
                 lwd = 0.1,
                 outlier.shape = 3, 
                 outlier.size = 0.5, 
                 outlier.stroke = 0.1) +
    geom_point(data = final_model_feature_threshold[[disease]], 
               mapping = aes(x = feature, y = threshold),
               size = 1,
               shape = 8,
               stroke = 0.25,
               color = "red") +
    scale_x_discrete(labels = function(x) scales::label_wrap(35)(x)) +
    labs(title = gsub("Cancer$", " Cancer", disease), 
         x = "Feature", 
         y = "Threshold NES") + 
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          # strip.background = element_rect(color = "black", linewidth = 0.25,),
          # strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
          text = element_text(size = 4), 
          plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, "cm"),
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, "cm"),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
}


if(!dir.exists("OutputFiles/Plots_publication/Predictive_model/")){ dir.create("OutputFiles/Plots_publication/Predictive_model/") }

tiff(paste0("OutputFiles/Plots_publication/Predictive_model/model_selected_thresholds_", drug_target_type, ".tiff"), 
     width = 25,
     height = 10,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list)

dev.off()


#####


# Save the thresholds as table
final_model_feature_threshold <- bind_rows(final_model_feature_threshold, .id = "disease")
colnames(final_model_feature_threshold) <- c("Disease", "Feature_name", "Threshold")
final_model_feature_threshold$Feature_type <- gsub("^\\[(.*)\\] .*", "\\1", final_model_feature_threshold$Feature_name)
final_model_feature_threshold$Feature_type <- factor(final_model_feature_threshold$Feature_type, levels = c("DISEASE", "ADR"), labels = c("Efficacy", "Safety"))
final_model_feature_threshold <- final_model_feature_threshold[, c("Disease", "Feature_type", "Feature_name", "Threshold")]
# final_model_feature_threshold$Threshold <- round(final_model_feature_threshold$Threshold , 3)
  
  
if(!dir.exists("OutputFiles/Tables_publication/")){ dir.create("OutputFiles/Tables_publication/", recursive = TRUE) }
write.csv(final_model_feature_threshold, "OutputFiles/Tables_publication/Predictive_model_feature_threshold.csv", row.names = FALSE)


#####


# Plot the final thresholds as barplot 

plot_list <- list()
final_model_feature_threshold <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  final_model_feature_threshold[[disease]] <- model$feature_threshold

  plot_list[[disease]] <- ggplot() +
    geom_bar(data = final_model_feature_threshold[[disease]], 
             mapping = aes(x = feature, y = threshold), 
             stat = "identity",
             width = 0.5, 
             lwd = 0.1) +
    scale_x_discrete(labels = function(x) scales::label_wrap(35)(x)) +
    labs(title = gsub("Cancer$", " Cancer", disease), 
         x = "Feature", 
         y = "Threshold NES") + 
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          # strip.background = element_rect(color = "black", linewidth = 0.25,),
          # strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
          text = element_text(size = 4), 
          plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, "cm"),
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, "cm"),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
}


if(!dir.exists("OutputFiles/Plots_publication/Predictive_model/")){ dir.create("OutputFiles/Plots_publication/Predictive_model/") }

tiff(paste0("OutputFiles/Plots_publication/Predictive_model/model_selected_thresholds_", drug_target_type, "_barplot.tiff"), 
     width = 25,
     height = 10,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list)

dev.off()


#####


# Plot the accuracy of the model 

final_model_accuracy <- list()
model_accuracy <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  
  final_model_accuracy[[disease]] <- model$final_model_accuracy
  
  tmp1 <- model$internal_CV_results$model_accuracy
  tmp1 <- lapply(tmp1, function(x){ bind_rows(x, .id = "group") })
  model_accuracy[[disease]] <- bind_rows(tmp1, .id = "fold")
  
}

final_model_accuracy <- bind_rows(final_model_accuracy, .id = "disease")
model_accuracy <- bind_rows(model_accuracy, .id = "disease")

model_accuracy <- left_join(x = model_accuracy, y = final_model_accuracy, by = c("disease" = "disease", ".metric" = ".metric"), suffix = c("_cv", "_final"))

final_model_accuracy <- final_model_accuracy[final_model_accuracy$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
final_model_accuracy$.metric <- factor(final_model_accuracy$.metric, 
                                       levels = c("bal_accuracy", "sensitivity", "specificity"),
                                       labels = c("Balanced Accuracy", "Sensitivity", "Specificity"))


model_accuracy <- model_accuracy[model_accuracy$.metric %in% c("bal_accuracy", "sensitivity", "specificity"), ]
model_accuracy$.metric <- factor(model_accuracy$.metric, 
                                 levels = c("bal_accuracy", "sensitivity", "specificity"),
                                 labels = c("Balanced Accuracy", "Sensitivity", "Specificity"))
model_accuracy$group <- factor(model_accuracy$group, levels = c("Train", "Test"))

if(!dir.exists("OutputFiles/Plots_publication/Predictive_model/")){ dir.create("OutputFiles/Plots_publication/Predictive_model/") }

tiff(paste0("OutputFiles/Plots_publication/Predictive_model/model_accuracy_", drug_target_type, ".tiff"), 
     width = 8,
     height = 6,
     units = "cm", compression = "lzw", res = 1200)

ggplot(data = model_accuracy) +
  geom_boxplot(mapping = aes(x = disease, y = .estimate_cv, fill = group),
               width = 0.5, 
               lwd = 0.05,
               outlier.shape = 3, 
               outlier.size = 0.5, 
               outlier.stroke = 0.1) +
  geom_hline(yintercept = 0.75,
             linewidth = 0.1,
             linetype = "dotted",
             color = "#006400") +
  geom_point(mapping = aes(x = disease, y = .estimate_final),
             size = 0.5,
             shape = 8,
             stroke = 0.1,
             color = "red") +
  facet_grid(.metric ~ .) + 
  scale_fill_manual(values=c("Test" = "#FF9F00", "Train" = "#007FFF")) +
  scale_x_discrete(labels = function(x){ gsub("Cancer$", " Cancer", x)  } ) +
  labs(x = "Disease", 
       y = "Score",
       fill = "Data: ") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4.5, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        # plot.title = element_text(size = 4, hjust = 0.5, face = "plain"),
        axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 2)),
        legend.key = element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.key.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(size = 4, margin = margin(l = 1)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()


#####


print(warnings())