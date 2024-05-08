
library(tidyverse)
library(ggpubr)

drug_target_type <- "known"



#####


# Plot the threshold selected in during CV and the final threshold
plot_list <- list()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
  final_model_feature_threshold <- model$feature_threshold
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
    geom_point(data = final_model_feature_threshold, 
               mapping = aes(x = feature, y = threshold),
               size = 1,
               shape = 8,
               stroke = 0.25,
               color = "red") +
    scale_x_discrete(labels = function(x) scales::label_wrap(35)(x)) +
    labs(title = disease, 
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
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
}


if(!dir.exists("OutputFiles/Plots_publication/Predictive_model/")){ dir.create("OutputFiles/Plots_publication/Predictive_model/") }

tiff(paste0("OutputFiles/Plots_publication/Predictive_model/model_selected_thresholds_", drug_target_type, ".tiff"), 
     width = 20,
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

final_model_accuracy <- final_model_accuracy[final_model_accuracy$.metric %in% c("f_meas", "sensitivity", "specificity"), ]
final_model_accuracy$.metric <- factor(final_model_accuracy$.metric, 
                                       levels = c("f_meas", "sensitivity", "specificity"),
                                       labels = c("F1", "Sensitivity", "Specificity"))


model_accuracy <- model_accuracy[model_accuracy$.metric %in% c("f_meas", "sensitivity", "specificity"), ]
model_accuracy$.metric <- factor(model_accuracy$.metric, 
                                 levels = c("f_meas", "sensitivity", "specificity"),
                                 labels = c("F1", "Sensitivity", "Specificity"))
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
  geom_point(mapping = aes(x = disease, y = .estimate_final),
             size = 0.5,
             shape = 8,
             stroke = 0.1,
             color = "red") +
  facet_grid(rows = vars(model_accuracy$.metric)) +
  scale_fill_manual(values=c("Test" = "#FF9F00", "Train" = "#007FFF")) +
  scale_x_discrete(labels = function(x) scales::label_wrap(20)(x)) +
  labs(x = "Disease", 
       y = "Score",
       fill = "Data") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        # plot.title = element_text(size = 4, hjust = 0.5, face = "plain"),
        axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()


#####

































# plot_list <- list()
# 
# for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
#   
#   model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
#   
#   final_model_feature_level_accuracy <- model$final_model_feature_level_accuracy
#   final_model_feature_level_accuracy <- pivot_longer(data = final_model_feature_level_accuracy, 
#                                                      cols = c("specificity", "sensitivity"), 
#                                                      names_to = "metric", 
#                                                      values_to = "score")
#   
#   
#   cv_feature_level_accuracy <- model$internal_CV_results$feature_level_accuracy
#   cv_feature_level_accuracy <- pivot_longer(data = cv_feature_level_accuracy, 
#                                             cols = c("specificity", "sensitivity"), 
#                                             names_to = "metric", 
#                                             values_to = "score")
#   
#   cv_feature_level_accuracy <- left_join(x = cv_feature_level_accuracy, 
#                                          y = final_model_feature_level_accuracy, 
#                                          by = c("feature" = "feature", "metric" = "metric"), 
#                                          suffix = c("_cv", "_final"))
#   
#   
#   plot_list[[disease]] = ggplot(data = cv_feature_level_accuracy) +
#     geom_boxplot(mapping = aes(x = feature, y = score_cv, fill = group),
#                  width = 0.5, 
#                  lwd = 0.05,
#                  outlier.shape = 3, 
#                  outlier.size = 0.5, 
#                  outlier.stroke = 0.1) +
#     geom_point(mapping = aes(x = feature, y = score_final),
#                size = 0.5,
#                shape = 8,
#                stroke = 0.1,
#                color = "red") +
#     facet_grid(rows = vars(cv_feature_level_accuracy$metric)) +
#     scale_fill_manual(values  =c("Test" = "#E69F00", "Train" = "#56B4E9")) +
#     scale_x_discrete(labels = function(x) scales::label_wrap(20)(x)) +
#     labs(title = disease,
#          x = "Disease", 
#          y = "Score",
#          fill = "Data") +
#     theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
#           panel.grid = element_blank(),
#           panel.spacing = unit(0.1, "cm"),
#           strip.background = element_rect(color = "black", linewidth = 0.25,),
#           strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
#           text = element_text(size = 4), 
#           plot.title = element_text(size = 4, hjust = 0.5, face = "bold"),
#           axis.text.x = element_text(size = 2.5, angle = 45, vjust = 1, hjust = 1), 
#           axis.ticks = element_line(colour = "black", linewidth = 0.2),
#           legend.position = "bottom",
#           legend.key = element_blank(),
#           legend.key.size = unit(0.1, 'cm'),
#           legend.text = element_text(size = 4),
#           legend.margin = margin(1,1,1,1),
#           legend.box.spacing = unit(0.1, 'cm'),
#           legend.box.background = element_rect(colour = "black", linewidth = 0.1))
#   
#   # print( plot_list[[disease]] )
#   
# }
# 
# 
# 
# # if(!dir.exists("OutputFiles/Plots_publication/Predictive_model/")){ dir.create("OutputFiles/Plots_publication/Predictive_model/") }
# # 
# # tiff(paste0("OutputFiles/Plots_publication/Predictive_model/model_feature_level_accuracy_", drug_target_type, ".tiff"),
# #      width = 30,
# #      height = 20,
# #      units = "cm", compression = "lzw", res = 1200)
# # 
# # ggarrange(plotlist = plot_list, common.legend = FALSE)
# # 
# # 
# # dev.off()
# 
# 
# 
# 
# 
# 
# 
# # Combine plots using ggarrange
# ggarrange(plotlist = plot_list, common.legend = FALSE,
#           nrow = 2, ncol = 6, align = "v")
