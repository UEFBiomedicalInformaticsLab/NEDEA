set.seed(5081)



# Script for plotting the model accuracy 



# Load libraries
library(optparse)
library(openxlsx)
library(tidyverse)
library(svglite)



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

# Define global options for this script 
disease <- opt$disease

cat(paste0("\n\nPlotting accuracy for: ", disease, "\n\n"))





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
    names(none_model_stats[[tmp]]) <- gsub(pattern = "\\.", 
                                           replacement = paste0("_", prox_comp, "."), 
                                           x = names(none_model_stats[[tmp]]), )
  }
}
rm(tmp)

none_model_stats <- unlist(none_model_stats, recursive = FALSE)
none_model_stats <- bind_rows(none_model_stats, .id = "model")
none_model_stats <- separate(none_model_stats, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

none_model_stats <- none_model_stats[, c("featureType", "model", "Fold", 
                                         "PRAUC_train", "PRAUC_test",
                                         "F1_train", "F1_test", 
                                         "BalancedAccuracy_train", "BalancedAccuracy_test")]
none_model_stats$imbalance <- "none"



# Merge all model stats and rearrange for plotting
model_stats <- none_model_stats 

model_stats <- reshape(model_stats, direction = "long",
                       varying = c("PRAUC_train", "F1_train", "BalancedAccuracy_train",
                                   "PRAUC_test", "F1_test", "BalancedAccuracy_test"),
                       v.names = "value",
                       timevar = "scoreType",
                       times = c("PRAUC_train", "F1_train", "BalancedAccuracy_train",
                                 "PRAUC_test", "F1_test", "BalancedAccuracy_test"))

rownames(model_stats) <- NULL
model_stats$value <- as.numeric(model_stats$value)

model_stats <- separate(model_stats, col = "scoreType", into = c("scoreType", "scoreType_class"), sep = "_")





# Plot the model test accuracy for selected features

features_to_plot <- c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene",
                      "DrugDrug_BbsiProx_separation", "DrugDisease_BbsiProx_separation",
                      "DrugAdr_BbsiProx_separation", "DrgDisAdr_BbsiProx_separation",
                      "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet", "miscGeneSet")
select_model_stats <- model_stats[(model_stats$scoreType_class == "test"),]
select_model_stats <- select_model_stats[(select_model_stats$featureType %in% features_to_plot),]

select_model_stats$featureType <- gsub("Dis2Gene", "Efficacy", select_model_stats$featureType)
select_model_stats$featureType <- gsub("WdrlAdr2Gene", "Safety", select_model_stats$featureType)
select_model_stats$featureType <- gsub("CombDisAdr2Gene", "Combined-efficacy-safety", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugDrug_BbsiProx_separation", "Drug-Drug (Sep.)", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugAdr_BbsiProx_separation", "Drug-ADR (Sep.)", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrugDisease_BbsiProx_separation", "Drug-Disease (Sep.)", select_model_stats$featureType)
select_model_stats$featureType <- gsub("DrgDisAdr_BbsiProx_separation", "Drug-Disease-ADR (Sep.)", select_model_stats$featureType)
select_model_stats$featureType <- gsub("keggPath", "KEGG pathways", select_model_stats$featureType)
select_model_stats$featureType <- gsub("SMPDbPath_DrugAction", "Drug action pathways", select_model_stats$featureType)
select_model_stats$featureType <- gsub("SMPDbPath_DrugMet", "Drug metabolism pathways", select_model_stats$featureType)
select_model_stats$featureType <- gsub("miscGeneSet", "Misc. gene sets", select_model_stats$featureType)

select_model_stats$featureType <- factor(x = select_model_stats$featureType,
                                         levels = c("Efficacy", "Safety", "Combined-efficacy-safety", 
                                                    "Drug-Drug (Sep.)", "Drug-Disease (Sep.)", "Drug-ADR (Sep.)", "Drug-Disease-ADR (Sep.)", 
                                                    "KEGG pathways", "Drug action pathways", "Drug metabolism pathways", 
                                                    "Misc. gene sets"))
select_model_stats$imbalance <- factor(x = select_model_stats$imbalance,
                                       levels = c("none"))

select_model_stats <- na.exclude(select_model_stats) # Some scores are NA if while calculating ratio the denominator is 0


if(!dir.exists(paste0("OutputFiles/Plots/", disease))){
  dir.create(paste0("OutputFiles/Plots/", disease), recursive = TRUE)
}

svglite(paste0("OutputFiles/Plots/", disease, "/ModelAccuracy_Test_", disease, ".svg"), 
        width = 12, 
        height = 5)

ggplot(select_model_stats, aes(x = model, y = value, fill = imbalance)) +
  geom_boxplot(width = 0.5, lwd = 0.1, outlier.shape = NA) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        strip.background = element_rect(color = "black", size = 0.1,),
        text = element_text(size = 5), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.ticks = element_line(colour = "black", size = 0.1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 5),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", size = 0.1)) +
  scale_fill_manual(values = c("#ff7f50", "#06b8b9", "#9a6324", "#911eb4")) +
  ggtitle(paste0("Model test accuracy for ", disease, " drug combinations")) +
  xlab("Models") + ylab("Accuracy measures") +
  labs(colour = "Imbalance : ") +
  facet_grid(rows = vars(scoreType), cols = vars(featureType))+
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "red", size = 0.1)
dev.off()



print(warnings())