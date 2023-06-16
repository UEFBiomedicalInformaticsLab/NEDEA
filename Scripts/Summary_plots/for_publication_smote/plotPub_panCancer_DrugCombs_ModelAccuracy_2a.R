# Script for plotting the model accuracy (for publication: comparision of combined-efficacy-safety vs Barabasi separation measure)



# Load libraries
library(openxlsx)
library(tidyverse)
library(svglite)





model_stats <- data.frame()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
    cat(paste0("\nReading data for: ", disease, "\n"))

    # Get the statistics for unbalanced models
    files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                        pattern = "^models_none_[a-zA-Z_]+.xlsx", 
                        ignore.case = TRUE, full.names = TRUE)

    none_model_stats <- list()

    for(file in files){
      sheet_names_none <- getSheetNames(file)

      ## Read files
      for(name in sheet_names_none){
        tmp <- strsplit(x = file, split = "\\/")[[1]][4]
        tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
        none_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
      }
      if(grepl(pattern = "BarabasiProx", x = file)){
        prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
        prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
        names(none_model_stats[[tmp]]) <- paste(prox_comp, names(none_model_stats[[tmp]]), sep = "_")
      }
    }
    rm(tmp)

    none_model_stats <- unlist(none_model_stats, recursive = FALSE)
    none_model_stats <- bind_rows(none_model_stats, .id = "model")
    none_model_stats <- separate(none_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

    none_model_stats <- none_model_stats[, c("featureType", "model", "Fold", 
                                               "PRAUC_train", "PRAUC_test",
                                                "BalancedAccuracy_train", "BalancedAccuracy_test",
                                                "Precision_train", "Precision_test",
                                                "Recall_train", "Recall_test",
                                                "F1_train", "F1_test")]
    none_model_stats$imbalance <- "none"





    # Get the statistics for SMOTE models
    files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                        pattern = "^models_SMOTE_[a-zA-Z_]+.xlsx", 
                        ignore.case = TRUE, full.names = TRUE)

    smote_model_stats <- list()

    for(file in files){
      sheet_names_smote <- getSheetNames(file)

      ## Read files
      for(name in sheet_names_smote){
        tmp <- strsplit(x = file, split = "\\/")[[1]][4]
        tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
        smote_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
      }
      if(grepl(pattern = "BarabasiProx", x = file)){
        prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
        prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
        names(smote_model_stats[[tmp]]) <- paste(prox_comp, names(smote_model_stats[[tmp]]), sep = "_")
      }
    }
    rm(tmp)

    smote_model_stats <- unlist(smote_model_stats, recursive = FALSE)
    smote_model_stats <- bind_rows(smote_model_stats, .id = "model")
    smote_model_stats <- separate(smote_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

    smote_model_stats <- smote_model_stats[, c("featureType", "model", "Fold", 
                                               "PRAUC_train", "PRAUC_test",
                                                "BalancedAccuracy_train", "BalancedAccuracy_test",
                                                "Precision_train", "Precision_test",
                                                "Recall_train", "Recall_test",
                                                "F1_train", "F1_test")]
    smote_model_stats$imbalance <- "SMOTE"





    # Get the statistics for upSample models
    files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                        pattern = "^models_upSample_[a-zA-Z_]+.xlsx", 
                        ignore.case = TRUE, full.names = TRUE)

    upSample_model_stats <- list()

    for(file in files){
      sheet_names_upSample <- getSheetNames(file)

      ## Read files
      for(name in sheet_names_upSample){
        tmp <- strsplit(x = file, split = "\\/")[[1]][4]
        tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
        upSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
      }
      if(grepl(pattern = "BarabasiProx", x = file)){
        prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
        prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
        names(upSample_model_stats[[tmp]]) <- paste(prox_comp, names(upSample_model_stats[[tmp]]), sep = "_")
      }
    }
    rm(tmp)

    upSample_model_stats <- unlist(upSample_model_stats, recursive = FALSE)
    upSample_model_stats <- bind_rows(upSample_model_stats, .id = "model")
    upSample_model_stats <- separate(upSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

    upSample_model_stats <- upSample_model_stats[, c("featureType", "model", "Fold", 
                                               "PRAUC_train", "PRAUC_test",
                                                "BalancedAccuracy_train", "BalancedAccuracy_test",
                                                "Precision_train", "Precision_test",
                                                "Recall_train", "Recall_test",
                                                "F1_train", "F1_test")]
    upSample_model_stats$imbalance <- "upSample"





    # Get the statistics for downSample models
    files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                        pattern = "^models_downSample_[a-zA-Z_]+.xlsx", 
                        ignore.case = TRUE, full.names = TRUE)

    downSample_model_stats <- list()

    for(file in files){
      sheet_names_downSample <- getSheetNames(file)

      ## Read files
      for(name in sheet_names_downSample){
        tmp <- strsplit(x = file, split = "\\/")[[1]][4]
        tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
        downSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
      }
      if(grepl(pattern = "BarabasiProx", x = file)){
        prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
        prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
        names(downSample_model_stats[[tmp]]) <- paste(prox_comp, names(downSample_model_stats[[tmp]]), sep = "_")
      }
    }
    rm(tmp)

    downSample_model_stats <- unlist(downSample_model_stats, recursive = FALSE)
    downSample_model_stats <- bind_rows(downSample_model_stats, .id = "model")
    downSample_model_stats <- separate(downSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

    downSample_model_stats <- downSample_model_stats[, c("featureType", "model", "Fold", 
                                               "PRAUC_train", "PRAUC_test",
                                                "BalancedAccuracy_train", "BalancedAccuracy_test",
                                                "Precision_train", "Precision_test",
                                                "Recall_train", "Recall_test",
                                                "F1_train", "F1_test")]
    downSample_model_stats$imbalance <- "downSample"





    # Merge all model stats and rearrange for plotting
    tmp <- rbind(none_model_stats, smote_model_stats, upSample_model_stats, downSample_model_stats) 
    tmp$disease <- disease
    model_stats <- rbind(model_stats, tmp)
  
}
saveRDS(model_stats, "Scripts/Summary_plots/for_publication/model_stats_2.rds")
# model_stats <- readRDS("Scripts/Summary_plots/for_publication/model_stats_2.rds")


model_stats <- reshape(model_stats, direction = "long",
                       varying = c("PRAUC_train", "BalancedAccuracy_train", "Precision_train", "Recall_train", "F1_train",
                                    "PRAUC_test", "BalancedAccuracy_test", "Precision_test", "Recall_test", "F1_test"),
                       v.names = "value",
                       timevar = "scoreType",
                       times = c("PRAUC_train", "BalancedAccuracy_train", "Precision_train", "Recall_train", "F1_tn",
                                    "PRAUC_test", "BalancedAccuracy_test", "Precision_test", "Recall_test", "F1_test"))

model_stats <- model_stats[model_stats$imbalance == "SMOTE", ]

rownames(model_stats) <- NULL
model_stats$value <- as.numeric(model_stats$value)

model_stats <- separate(model_stats, col = "scoreType", into = c("scoreType", "scoreType_class"), sep = "_")




# Plot the model test accuracy for selected features

if(!dir.exists("OutputFiles/Plots/Publication")){
  dir.create("OutputFiles/Plots/Publication", recursive = TRUE)
}



features_to_plot <- unique(model_stats$featureType)
features_to_plot <- features_to_plot[grep("_BbsiProx_separation$", features_to_plot)]
features_to_plot <- c("CombDisAdr2Gene", features_to_plot)
select_model_stats <- model_stats[(model_stats$scoreType_class == "test"),]
select_model_stats <- select_model_stats[(select_model_stats$value > 0.5 ),] # removing insignificant scores
select_model_stats <- select_model_stats[(select_model_stats$featureType %in% features_to_plot),]


select_model_stats$featureType <- gsub("CombDisAdr2Gene", "Combined-efficacy-safety", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugAdr_BbsiProx_separation", "Drug-ADR", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugDisease_BbsiProx_separation", "Drug-Disease", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugDrug_BbsiProx_separation", "Drug-Drug", select_model_stats$featureType)
select_model_stats$disease <- gsub("Cancer$", "", select_model_stats$disease)


select_model_stats$featureType <- factor(x = select_model_stats$featureType,
                                            levels = c("Combined-efficacy-safety", "Drug-Drug", "Drug-Disease", "Drug-ADR"))
select_model_stats$imbalance <- factor(x = select_model_stats$imbalance,
                                            levels = c("none", "SMOTE", "upSample", "downSample"))

select_model_stats <- na.exclude(select_model_stats) # Some scores are NA if while calculating ratio the denominator is 0





# svglite("OutputFiles/Plots/Publication/panCancer_ModelAccuracy_Test.svg", width = 6, height = length(unique(model_stats$disease)))
tiff("OutputFiles/Plots/Publication/panCancer_ModelAccuracy_Test_2a.tiff", 
     width = 15, height = 10, 
     units = "cm", compression = "lzw", res = 1200)

ggplot(select_model_stats, aes(x = disease, y = value, fill = model)) + 
  geom_boxplot(lwd = 0.1, outlier.shape = NA, position = position_dodge2(preserve = "single")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", size = 0.1,),
        strip.text = element_text(margin = margin(1,1,1,1)),
        text = element_text(size = 3), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", size = 0.1)) +
  scale_fill_manual(values = c("#008080", "#ffa500", "#00ff7f", "#00bfff", "#deb887")) + 
  # ggtitle(paste0("Model test accuracy for drug combinations")) +
  xlab("Cancers") + ylab("Accuracy scores") + 
  labs(colour = "Models : ") +
  facet_grid(rows = vars(scoreType), cols = vars(featureType)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "red", size = 0.05)
dev.off()



print(warnings())