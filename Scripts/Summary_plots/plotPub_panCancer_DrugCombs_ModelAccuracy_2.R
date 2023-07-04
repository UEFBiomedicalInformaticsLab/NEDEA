set.seed(5081)



# Script for plotting the model accuracy (for publication: comparision of combined-efficacy-safety vs Barabasi separation measure)



# Load libraries
library(openxlsx)
library(tidyverse)
library(svglite)
library(RColorBrewer)



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
    if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
      prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
      prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
      prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
      names(none_model_stats[[tmp]]) <- gsub(pattern = "\\.", replacement = paste0("_", prox_comp, "."), x = names(none_model_stats[[tmp]]), )
    }
  }   
  rm(tmp)
  
  none_model_stats <- unlist(none_model_stats, recursive = FALSE)
  none_model_stats <- bind_rows(none_model_stats, .id = "model")
  none_model_stats <- separate(none_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")
  
  none_model_stats <- none_model_stats[, c("featureType", "model", "Fold", 
                                           "PRAUC_train", "PRAUC_test",
                                           "F1_train", "F1_test")]
  none_model_stats <- none_model_stats[none_model_stats$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]
  none_model_stats$imbalance <- "none"
  

  
  # Merge all model stats and rearrange for plotting
  tmp <- none_model_stats 
  tmp$disease <- disease
  model_stats <- rbind(model_stats, tmp)
  
}

if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}

saveRDS(model_stats, "OutputFiles/Plots/model_stats_2.rds")
# model_stats <- readRDS("OutputFiles/Plots/model_stats_2.rds")

model_stats <- model_stats[, !colnames(model_stats) %in% c("BalancedAccuracy_train", "BalancedAccuracy_test", 
                                                           "Precision_train", "Precision_test", 
                                                           "Recall_train", "Recall_test")]
model_stats <- reshape(model_stats, direction = "long",
                       varying = c("PRAUC_train", "F1_train",
                                   "PRAUC_test", "F1_test"),
                       v.names = "value",
                       timevar = "scoreType",
                       times = c("PRAUC_train", "F1_train",
                                 "PRAUC_test", "F1_test"))


rownames(model_stats) <- NULL
model_stats$value <- as.numeric(model_stats$value)

model_stats <- separate(model_stats, col = "scoreType", into = c("scoreType", "scoreType_class"), sep = "_")




# Plot the model test accuracy for selected features
features_to_plot <- unique(model_stats$featureType)
features_to_plot <- features_to_plot[grep("_BbsiProx_separation$", features_to_plot)]
features_to_plot <- features_to_plot[grep("^DrugDrug_", features_to_plot, invert = TRUE)]
features_to_plot <- c("CombDisAdr2Gene", features_to_plot)
select_model_stats <- model_stats[(model_stats$scoreType_class == "test"),]
# select_model_stats <- select_model_stats[(select_model_stats$value > 0.5 ),] # removing insignificant scores
select_model_stats <- select_model_stats[(select_model_stats$featureType %in% features_to_plot),]


select_model_stats$featureType <- gsub("CombDisAdr2Gene", "Combined-efficacy-safety", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugAdr_BbsiProx_separation", "Drug-ADR", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugDisease_BbsiProx_separation", "Drug-Disease", select_model_stats$featureType)
# select_model_stats$featureType <- gsub("DrugDrug_BbsiProx_separation", "Drug-Drug", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrgDisAdr_BbsiProx_separation", "Drug-Disease-ADR", select_model_stats$featureType)
select_model_stats$disease <- gsub("Cancer$", "", select_model_stats$disease)


select_model_stats$featureType <- factor(x = select_model_stats$featureType,
                                         levels = c("Combined-efficacy-safety", "Drug-Drug", "Drug-Disease", "Drug-ADR", "Drug-Disease-ADR"))
select_model_stats$imbalance <- factor(x = select_model_stats$imbalance,
                                       levels = c("none"))

select_model_stats <- na.exclude(select_model_stats) # Some scores are NA if while calculating ratio the denominator is 0





tiff("OutputFiles/Plots/panCancer_ModelAccuracy_Test_2.tiff", 
     width = 17, height = length(features_to_plot) * 4, 
     units = "cm", compression = "lzw", res = 1200)

ggplot(select_model_stats, aes(x = disease, y = value, fill = model)) + 
  geom_boxplot(lwd = 0.2, outlier.shape = NA, position = position_dodge(preserve = "single")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", size = 0.25,),
        strip.text = element_text(margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", size = 0.25)) +
  scale_fill_manual(values = brewer.pal(4, "Set2")) + 
  # ggtitle(paste0("Model test accuracy for drug combinations")) +
  xlab("Cancers") + ylab("Accuracy scores") + 
  labs(colour = "Models : ") +
  facet_grid(cols = vars(scoreType), rows = vars(featureType)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "red", size = 0.1)
dev.off()



print(warnings())