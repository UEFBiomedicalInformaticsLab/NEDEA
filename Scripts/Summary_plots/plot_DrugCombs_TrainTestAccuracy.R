set.seed(5081)
rm(list = ls())



# Script for plotting the model train vs test accuracy 



# Load libraries
library(optparse)
library(openxlsx)
library(tidyverse)
library(svglite)






# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
    make_option(c("--metric"), type = "character", default = "F1", 
              help = "The accuracy metric to be plotted. Options: BalancedAccuracy, F1, PRAUC. Default: F1.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$metric %in% c("BalancedAccuracy", "F1", "PRAUC")){
  print_help(opt_parser)
  stop("--metric should be either of BalancedAccuracy, F1, PRAUC.", call.=FALSE)
}

# Define global options for this script 
disease <- opt$disease
metric <- opt$metric

cat(paste0("\n\nPlotting accuracy for: ", disease, "\n\n"))





# Get the statistics for unbalanced models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "models_none_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

none_model_stats <- list()

for(file in files){
  sheet_names_none <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_none){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    none_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
}
rm(tmp)

none_model_stats <- unlist(none_model_stats, recursive = FALSE)
none_model_stats <- bind_rows(none_model_stats, .id = "model")
none_model_stats <- separate(none_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

none_model_stats <- none_model_stats[, c("featureType", "model", "Resample", 
                                           "PR_AUC_train", "PR_AUC_test",
                                            "F1_train", "F1_test", 
                                            "BalancedAccuracy_train", "BalancedAccuracy_test")]
none_model_stats$imbalance <- "none"





# Get the statistics for SMOTE models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "models_SMOTE_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

smote_model_stats <- list()

for(file in files){
  sheet_names_smote <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_smote){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    smote_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
}
rm(tmp)

smote_model_stats <- unlist(smote_model_stats, recursive = FALSE)
smote_model_stats <- bind_rows(smote_model_stats, .id = "model")
smote_model_stats <- separate(smote_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

smote_model_stats <- smote_model_stats[, c("featureType", "model", "Resample", 
                                           "PR_AUC_train", "PR_AUC_test",
                                            "F1_train", "F1_test", 
                                            "BalancedAccuracy_train", "BalancedAccuracy_test")]
smote_model_stats$imbalance <- "SMOTE"





# Get the statistics for upSample models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "models_upSample_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

upSample_model_stats <- list()

for(file in files){
  sheet_names_upSample <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_upSample){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    upSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
}
rm(tmp)

upSample_model_stats <- unlist(upSample_model_stats, recursive = FALSE)
upSample_model_stats <- bind_rows(upSample_model_stats, .id = "model")
upSample_model_stats <- separate(upSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

upSample_model_stats <- upSample_model_stats[, c("featureType", "model", "Resample", 
                                           "PR_AUC_train", "PR_AUC_test",
                                            "F1_train", "F1_test", 
                                            "BalancedAccuracy_train", "BalancedAccuracy_test")]
upSample_model_stats$imbalance <- "upSample"





# Get the statistics for downSample models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "models_downSample_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

downSample_model_stats <- list()

for(file in files){
  sheet_names_downSample <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_downSample){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    downSample_model_stats[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
}
rm(tmp)

downSample_model_stats <- unlist(downSample_model_stats, recursive = FALSE)
downSample_model_stats <- bind_rows(downSample_model_stats, .id = "model")
downSample_model_stats <- separate(downSample_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

downSample_model_stats <- downSample_model_stats[, c("featureType", "model", "Resample", 
                                           "PR_AUC_train", "PR_AUC_test",
                                            "F1_train", "F1_test", 
                                            "BalancedAccuracy_train", "BalancedAccuracy_test")]
downSample_model_stats$imbalance <- "downSample"





# Merge all model stats and rearrange for plotting
model_stats <- rbind(none_model_stats, smote_model_stats, upSample_model_stats, downSample_model_stats) 

model_stats <- reshape(model_stats, direction = "long",
                       varying = c("PR_AUC_train", "F1_train", "BalancedAccuracy_train",
                                    "PR_AUC_test", "F1_test", "BalancedAccuracy_test"),
                       v.names = "value",
                       timevar = "scoreType",
                       times = c("PR_AUC_train", "F1_train", "BalancedAccuracy_train",
                                 "PR_AUC_test", "F1_test", "BalancedAccuracy_test"))

rownames(model_stats) <- NULL
model_stats$value <- as.numeric(model_stats$value)

model_stats$scoreType <- gsub("^PR_AUC", "PRAUC", model_stats$scoreType) # To be removed once all scripts are fixed

model_stats <- separate(model_stats, col = "scoreType", into = c("scoreType", "scoreType_class"), sep = "_")





# Plot train vs test accuracy
svglite(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/ModelAccuracy_TrainVsTest_",metric, "_", disease, ".svg"), width = 4, height = 4)
features_to_plot <- c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene",
                        "BbsiProx_separation", "SteinerTopol", 
                        "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet")
select_model_stats <- model_stats[(model_stats$scoreType == metric),]
select_model_stats <- select_model_stats[(select_model_stats$featureType %in% features_to_plot),]
select_model_stats$featureType <- factor(x = select_model_stats$featureType,
                                            levels = c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene",
                                                        "BbsiProx_separation", "SteinerTopol", 
                                                        "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet"))
select_model_stats$imbalance <- factor(x = select_model_stats$imbalance,
                                            levels = c("none", "SMOTE", "upSample", "downSample"))
select_model_stats$scoreType_class <- factor(select_model_stats$scoreType_class, levels = c("train", "test"))

select_model_stats <- na.exclude(select_model_stats) # Some scores are NA if while calculating ratio the denominator is 0

ggplot(select_model_stats, aes(x = model, y = value, fill = scoreType_class)) + 
  geom_boxplot(width = 0.5, lwd = 0.1, outlier.shape = NA) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        strip.background = element_rect(color = "black", size = 0.1,),
        text = element_text(size = 5), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.ticks = element_line(colour = "black", size = 0.1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 5),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", size = 0.1)) +
  scale_fill_manual(values = c("#ffe119", "#4363d8")) +
  ggtitle(paste0("Train vs Test ", metric, " for ", disease, " drug combinations")) +
  xlab("Models") + ylab("Accuracy measures") + 
  labs(colour = "Score type : ") +
  facet_grid(rows = vars(imbalance), cols = vars(featureType))
dev.off()