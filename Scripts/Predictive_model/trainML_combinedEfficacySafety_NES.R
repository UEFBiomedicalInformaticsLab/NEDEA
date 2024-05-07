set.seed(5081)


# Script to train predictive model



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(yardstick)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known",
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type



cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


#####


# Read the FGSEA result
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
fgsea_result <- fgsea_result[["combinedEfficacySafety"]]


# Read the drug combination category
drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
drugCombs_cat <- drugCombs_cat[, c("comb_name", "class_EffAdv")]
drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), ]


# Read the list of selected features
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"))
# selected_features <- selected_features[order(selected_features$feature), ]

# Sub-set the FGSEA results for only the selected drug combinations and features
fgsea_result_select <- fgsea_result[row.names(fgsea_result) %in% selected_features$feature, 
                                    colnames(fgsea_result) %in% drugCombs_cat$comb_name]
fgsea_result_select <- as.data.frame(t(fgsea_result_select))
fgsea_result_select$category <-  drugCombs_cat$class_EffAdv[match(row.names(fgsea_result_select), drugCombs_cat$comb_name)]


#####


# Create different folds of data to verify the stability
data_folds <- list()
for(i in seq(1:10)){
  set.seed(5081 + i*2)
  data_folds[[paste0("Repeat", i)]] <- caret::createFolds(y = fgsea_result_select$category,
                                                          k = 3,
                                                          list = TRUE,
                                                          returnTrain = TRUE)
}

data_folds <- unlist(data_folds, recursive = FALSE)
names(data_folds) <- gsub("\\.", "_", names(data_folds))
set.seed(5081)


#####


metrics_list <- list()
feature_level_accuracy <- data.frame()
metrics <- metric_set(accuracy, f_meas, sensitivity, specificity, recall, precision, bal_accuracy)
selected_threshold <- data.frame()

for(fold in names(data_folds)){
  
  train_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select[data_folds[[fold]], ]))
  test_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select[-data_folds[[fold]], ]))
  
  
  for(feature_name in selected_features$feature){
    
    # Prepare the data for fitting
    fit_data <- fgsea_result_select[data_folds[[fold]], c(feature_name, "category")]
    colnames(fit_data)[1] <- "feature"
    tmp1 <- c("Eff" = 1, "Adv" = 0)
    fit_data$category <- unname(tmp1[fit_data$category])
    rm(tmp1)
    
    
    # Fit the model
    model <- glm(category ~ feature , data = fit_data, family = "binomial")
    fit_data$predicted_probability <- predict(model, newdata = fit_data, type = "response")
    
    
    # Calculate the probability threshold with best separation
    fit_data$category <- factor(fit_data$category, levels = c("1", "0"))
    roc_curve_data <- roc_curve(data = fit_data, truth = category, predicted_probability)
    roc_curve_data$Youden_stat <- roc_curve_data$sensitivity + roc_curve_data$specificity - 1
    roc_curve_data <- roc_curve_data[roc_curve_data$Youden_stat == max(roc_curve_data$Youden_stat), ]
    if(nrow(roc_curve_data) > 1){
      roc_curve_data <- roc_curve_data[roc_curve_data$specificity == max(roc_curve_data$specificity), ]
    }
    best_probability_threshold <- roc_curve_data$.threshold
    
    
    # Get the actual feature threshold using model coefficients and the selected probability
    best_feature_threshold <- (log(best_probability_threshold / (1 - best_probability_threshold)) - as.numeric(coef(model)[1])) / as.numeric(coef(model)[2]) 
    selected_threshold <- rbind(selected_threshold, 
                                data.frame("fold" = fold, 
                                           "feature" = feature_name, 
                                           "threshold" = best_feature_threshold))
    
    
    # Generate predictions for training
    tmp1 <- fit_data
    
    if(grepl("^\\[DISEASE\\]", feature_name)){
      tmp1$predicted_category <- ifelse(fit_data$feature >= best_feature_threshold, 1, -1)
    }
    if(grepl("^\\[ADR\\]", feature_name)){
      tmp1$predicted_category <- ifelse(fit_data$feature >= best_feature_threshold, -1, 1)
    } 
    
    
    tmp2 <- tmp1[, c("category", "predicted_category")]
    tmp2$category <- factor(tmp2$category, levels = c("1", "0"), labels = c("Eff", "Adv"))
    tmp2$predicted_category <- factor(tmp2$predicted_category, levels = c("1", "-1"), labels = c("Eff", "Adv"))
    feature_level_accuracy <- rbind(feature_level_accuracy, 
                                    data.frame("fold" = fold, 
                                               "feature" = feature_name, 
                                               "group" = "Train", 
                                               "specificity" = specificity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category), 
                                               "sensitivity" = sensitivity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category)))
    
    
    tmp1 <- tmp1[, c("predicted_category"), drop = FALSE]
    colnames(tmp1) <- feature_name
    tmp1 <- rownames_to_column(tmp1, "drugCombs")
    train_prediction_df <- full_join(train_prediction_df, tmp1, by = "drugCombs")
    rm(tmp1)
    rm(tmp2)
    
    
    # Generate predictions for test
    predict_data <- fgsea_result_select[-data_folds[[fold]], c(feature_name, "category")]
    colnames(predict_data)[1] <- "feature"
    
    if(grepl("^\\[DISEASE\\]", feature_name)){
      predict_data$predicted_category <- ifelse(predict_data$feature >= best_feature_threshold, 1, -1)
    }
    if(grepl("^\\[ADR\\]", feature_name)){
      predict_data$predicted_category <- ifelse(predict_data$feature >= best_feature_threshold, -1, 1)
    } 
    
    
    tmp2 <- predict_data[, c("category", "predicted_category")]
    tmp2$category <- factor(tmp2$category, levels = c("Eff", "Adv"), labels = c("Eff", "Adv"))
    tmp2$predicted_category <- factor(tmp2$predicted_category, levels = c("1", "-1"), labels = c("Eff", "Adv"))
    feature_level_accuracy <- rbind(feature_level_accuracy, 
                                    data.frame("fold" = fold, 
                                               "feature" = feature_name, 
                                               "group" = "Test", 
                                               "specificity" = specificity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category), 
                                               "sensitivity" = sensitivity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category)))
    
    
    tmp1 <- predict_data[, c("predicted_category"), drop = FALSE]
    colnames(tmp1) <- feature_name
    tmp1 <- rownames_to_column(tmp1, "drugCombs")
    test_prediction_df <- full_join(test_prediction_df, tmp1, by = "drugCombs")
    rm(tmp1)
    rm(tmp2)
    
  }
  
  
  # Calculate accuracy metrics for training 
  train_prediction_df <- column_to_rownames(train_prediction_df, "drugCombs")

  train_prediction_df <- train_prediction_df %>%
    mutate( efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
            safety_score = rowMeans(select(., starts_with("[ADR]"))) ) %>% 
    mutate( final_score = rowMeans(select(., c(efficacy_score, safety_score))) )
  

  train_prediction_df$final_predicted_category <- ifelse(train_prediction_df$final_score > 0, "Eff", "Adv")
  train_prediction_df$final_predicted_category <- factor(train_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))
  
  train_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(train_prediction_df), drugCombs_cat$comb_name)]
  train_prediction_df$class_EffAdv <- factor(train_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))
  
  metrics_list[[fold]][["Train"]] <- metrics(data = train_prediction_df, truth = class_EffAdv, estimate = final_predicted_category)
  
  
  # Calculate accuracy metrics for test 
  test_prediction_df <- column_to_rownames(test_prediction_df, "drugCombs")

  test_prediction_df <- test_prediction_df %>%
    mutate( efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
            safety_score = rowMeans(select(., starts_with("[ADR]"))) ) %>% 
    mutate( final_score = rowMeans(select(., c(efficacy_score, safety_score))) )
  
  
  test_prediction_df$final_predicted_category <- ifelse(test_prediction_df$final_score > 0, "Eff", "Adv")
  test_prediction_df$final_predicted_category <- factor(test_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))
  
  test_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(test_prediction_df), drugCombs_cat$comb_name)]
  test_prediction_df$class_EffAdv <- factor(test_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))
  
  metrics_list[[fold]][["Test"]] <- metrics(data = test_prediction_df, truth = class_EffAdv, estimate = final_predicted_category)
  
}


#####


plot_list <- list()

# Plot the threshold
plot_list[["thresholds"]] <- ggplot(selected_threshold, aes(x = feature, y = threshold)) +
  geom_boxplot(fill = "grey", 
               width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3, 
               outlier.size = 0.5, 
               outlier.stroke = 0.1) +
  scale_x_discrete(labels = function(x) scales::label_wrap(30)(x)) +
  labs(x = "Feature", 
       y = "Threshold") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4), 
        axis.text.x = element_text(size = 3, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))


# Plot the metrics
plot_data <- unlist(metrics_list, recursive = FALSE)

plot_data <- lapply(metrics_list, function(x){ bind_rows(x, .id = "group") })
plot_data <- bind_rows(plot_data, .id = "fold")

plot_data$.metric <- factor(plot_data$.metric, levels = c("accuracy", "f_meas", "sensitivity", "specificity", "recall", "precision", "bal_accuracy", "ppv", "npv"))
plot_data$group <- factor(plot_data$group, levels = c("Train", "Test"))


plot_list[["metrics"]] <- ggplot(plot_data, aes(x = .metric, y = .estimate, fill = group)) +
  geom_boxplot(width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3, 
               outlier.size = 0.5, 
               outlier.stroke = 0.1) +
  ylim(0, 1) +
  scale_fill_manual(values=c("Test" = "#E69F00", "Train" = "#56B4E9")) +
  scale_x_discrete(labels = function(x) scales::label_wrap(20)(x)) +
  labs(x = "Metric", 
       y = "Score",
       fill = "Data") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4), 
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))


# Plot the per feature accuracy
plot_data <- feature_level_accuracy
plot_data <- pivot_longer(data = plot_data, 
                          cols = c("specificity", "sensitivity"), 
                          names_to = "metric", 
                          values_to = "score")
plot_data$group <- factor(plot_data$group, levels = c("Train", "Test"))


plot_list[["feature_level_accuracy"]] <- ggplot(plot_data, aes(x = feature, y = score, fill = group)) +
  geom_boxplot(width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3, 
               outlier.size = 0.5, 
               outlier.stroke = 0.1) +
  scale_fill_manual(values=c("Test" = "#E69F00", "Train" = "#56B4E9")) +
  facet_grid(rows = vars(plot_data$metric)) +
  scale_x_discrete(labels = function(x) scales::label_wrap(30)(x)) +
  labs(x = "Feature",
       y = "Score",
       fill = "Data") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        axis.text.x = element_text(size = 3, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))



if(!dir.exists("OutputFiles/Plots/Predictive_model/")){ dir.create("OutputFiles/Plots/Predictive_model/") }

tiff(paste0("OutputFiles/Plots/Predictive_model/model_characterization_", disease, "_", drug_target_type, ".tiff"), 
     width = 20,
     height = 10,
     units = "cm", compression = "lzw", res = 1200)

ggpubr::ggarrange(plot_list$feature_level_accuracy, 
                  ggpubr::ggarrange(plot_list$thresholds, plot_list$metrics, nrow = 2),
                  ncol = 2)

dev.off()
rm(list = c("plot_data", "plot_list"))


#####


# Train model with complete data
allData_prediction_df <- data.frame("drugCombs" = row.names(fgsea_result_select))
model_list <- list()
final_model_selected_threshold <- data.frame()
final_model_feature_level_accuracy <- data.frame()
plot_list <- list()

for(feature_name in selected_features$feature){
  
  # Prepare the data for fitting
  fit_data <- fgsea_result_select[, c(feature_name, "category")]
  colnames(fit_data)[1] <- "feature"
  tmp1 <- c("Eff" = 1, "Adv" = 0)
  fit_data$category <- unname(tmp1[fit_data$category])
  rm(tmp1)
  
  
  # Fit the model
  model <- glm(category ~ feature , data = fit_data, family = "binomial")
  model_list[[feature_name]] <- model
  
  
  # Prediction on the full data 
  fit_data$predicted_probability <- predict(model, newdata = fit_data, type = "response")
  
  
  # Calculate the probability threshold with best separation
  fit_data$category <- factor(fit_data$category, levels = c("1", "0"))
  roc_curve_data <- roc_curve(data = fit_data, truth = category, predicted_probability)
  roc_curve_data$Youden_stat <- roc_curve_data$sensitivity + roc_curve_data$specificity - 1
  roc_curve_data <- roc_curve_data[roc_curve_data$Youden_stat == max(roc_curve_data$Youden_stat), ]
  if(nrow(roc_curve_data) > 1){
    roc_curve_data <- roc_curve_data[roc_curve_data$specificity == max(roc_curve_data$specificity), ]
  }
  best_probability_threshold <- roc_curve_data$.threshold
  
  
  # Get the actual feature threshold using model coefficients and the selected probability
  best_feature_threshold <- (log(best_probability_threshold / (1 - best_probability_threshold)) - as.numeric(coef(model)[1])) / as.numeric(coef(model)[2]) ### RECHECK
  final_model_selected_threshold <- rbind(final_model_selected_threshold, 
                                          data.frame("feature" = feature_name, 
                                                     "threshold" = best_feature_threshold))
  
  
  # Generate predictions for all data
  tmp1 <- fit_data
  
  if(grepl("^\\[DISEASE\\]", feature_name)){
    tmp1$predicted_category <- ifelse(fit_data$feature >= best_feature_threshold, 1, -1)
  }
  if(grepl("^\\[ADR\\]", feature_name)){
    tmp1$predicted_category <- ifelse(fit_data$feature >= best_feature_threshold, -1, 1)
  } 
  
  
  # Calculate the feature level accuracy
  tmp2 <- tmp1[, c("category", "predicted_category")]
  tmp2$category <- factor(tmp2$category, levels = c("1", "0"), labels = c("Eff", "Adv"))
  tmp2$predicted_category <- factor(tmp2$predicted_category, levels = c(1, -1), labels = c("Eff", "Adv"))
  final_model_feature_level_accuracy <- rbind(final_model_feature_level_accuracy, 
                                              data.frame("feature" = feature_name, 
                                                         "specificity" = specificity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category), 
                                                         "sensitivity" = sensitivity_vec(data = tmp2, truth = tmp2$category, estimate = tmp2$predicted_category)))
  rm(tmp2)
  
  
  # Plot the partition
  plot_data <- tmp1
  plot_data$feature_name <- feature_name
  
  plot_data$category <- factor(plot_data$category, levels = c("1", "0"), labels = c("Eff", "Adv"))
  plot_data$predicted_category <- factor(plot_data$predicted_category, levels = c("1", "-1"), labels = c("Eff", "Adv"))
  
  
  plot_list[[feature_name]] <- ggplot(plot_data, aes(x = feature_name, y = feature, group = category)) +
    geom_boxplot(fill = "grey",
                 width = 0.5,
                 lwd = 0.1,
                 outliers = FALSE) +
    geom_jitter(aes(shape = category,
                    color = predicted_category),
                size = 1) +
    geom_hline(yintercept = best_feature_threshold,
               linewidth = 0.1,
               linetype = "dashed",
               color = "red") +
    scale_x_discrete(labels = function(x) scales::label_wrap(40)(x)) +
    scale_color_manual(values = c("Adv" = "#E69F00", "Eff" = "#56B4E9")) +
    scale_shape_manual(values = c("Adv" = 1, "Eff" = 3)) +
    labs(x = "Feature",
         y = "NES",
         shape = "Acual cat.",
         color = "Predicted cat.") +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          # strip.background = element_rect(color = "black", linewidth = 0.25,),
          # strip.text = element_text(size = 5, margin = margin(1,1,1,1)),
          text = element_text(size = 5),
          # plot.title = element_text(size = 4, hjust = 0.5, face = "plain"),
          axis.text.x = element_text(size = 5, angle = 0, vjust = 0, hjust = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.1))
  
  
  tmp1 <- tmp1[, c("predicted_category"), drop = FALSE]
  colnames(tmp1) <- feature_name
  tmp1 <- rownames_to_column(tmp1, "drugCombs")
  allData_prediction_df <- full_join(allData_prediction_df, tmp1, by = "drugCombs")
  rm(tmp1)
  
}


# Calculate accuracy metrics on the complete data
allData_prediction_df <- column_to_rownames(allData_prediction_df, "drugCombs")

allData_prediction_df <- allData_prediction_df %>%
  mutate( efficacy_score = rowMeans(select(., starts_with("[DISEASE]"))), 
          safety_score = rowMeans(select(., starts_with("[ADR]"))) ) %>% 
  mutate( final_score = rowMeans(select(., c(efficacy_score, safety_score))) )
  

allData_prediction_df$final_predicted_category <- ifelse(allData_prediction_df$final_score > 0, "Eff", "Adv")
allData_prediction_df$final_predicted_category <- factor(allData_prediction_df$final_predicted_category, levels = c("Eff", "Adv"))

allData_prediction_df$class_EffAdv <-  drugCombs_cat$class_EffAdv[match(row.names(allData_prediction_df), drugCombs_cat$comb_name)]
allData_prediction_df$class_EffAdv <- factor(allData_prediction_df$class_EffAdv, levels = c("Eff", "Adv"))

final_model_accuracy <- metrics(data = allData_prediction_df, truth = class_EffAdv, estimate = final_predicted_category)


# Save the data partition plot
if(!dir.exists("OutputFiles/Plots/Predictive_model/")){ dir.create("OutputFiles/Plots/Predictive_model/") }

tiff(paste0("OutputFiles/Plots/Predictive_model/training_data_partition_", disease, "_", drug_target_type, ".tiff"), 
     width = 20,
     height = 15,
     units = "cm", compression = "lzw", res = 1200)

ggpubr::ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "bottom")

dev.off()


#####


# Save the final model

final_model <- list("models" = model_list, 
                    "feature_threshold" = final_model_selected_threshold, 
                    "final_model_accuracy" = final_model_accuracy, 
                    "final_model_feature_level_accuracy" = final_model_feature_level_accuracy,
                    "internal_CV_results" = list("selected_threshold" = selected_threshold, 
                                                 "feature_level_accuracy" = feature_level_accuracy, 
                                                 "model_accuracy" = metrics_list)
)


if(!dir.exists("OutputFiles/Predictive_model/")){ dir.create("OutputFiles/Predictive_model/", recursive = TRUE) }
saveRDS(final_model, file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))



print(warnings())