rm(list = ls())



# Script for plotting the model parameters across resampling



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

cat(paste0("\n\nPlotting parameters for: ", disease, "\n\n"))





# Get the parameters for unbalanced models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "^models_none_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

none_model_param <- list()

for(file in files){
  sheet_names_none <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_none){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    none_model_param[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
  if(grepl(pattern = "BarabasiProx", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][5]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
    names(none_model_param[[tmp]]) <- paste(prox_comp, names(none_model_param[[tmp]]), sep = "_")
  }
}
rm(tmp)

none_model_param <- unlist(none_model_param, recursive = FALSE)
none_model_param <- bind_rows(none_model_param, .id = "model")
none_model_param <- separate(none_model_param, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

keep <- colnames(none_model_param)[grep("BestTune_", colnames(none_model_param))]
none_model_param <- none_model_param[, c("featureType", "model", "Fold", keep)]
none_model_param$imbalance <- "none"





# Get the parameters for SMOTE models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "^models_SMOTE_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

smote_model_param <- list()

for(file in files){
  sheet_names_smote <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_smote){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    smote_model_param[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
  if(grepl(pattern = "BarabasiProx", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][5]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
    names(smote_model_param[[tmp]]) <- paste(prox_comp, names(smote_model_param[[tmp]]), sep = "_")
  }
}
rm(tmp)

smote_model_param <- unlist(smote_model_param, recursive = FALSE)
smote_model_param <- bind_rows(smote_model_param, .id = "model")
smote_model_param <- separate(smote_model_param, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

keep <- colnames(smote_model_param)[grep("BestTune_", colnames(smote_model_param))]
smote_model_param <- smote_model_param[, c("featureType", "model", "Fold", keep)]
smote_model_param$imbalance <- "SMOTE"





# Get the parameters for upSample models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "^models_upSample_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

upSample_model_param <- list()

for(file in files){
  sheet_names_upSample <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_upSample){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    upSample_model_param[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
  if(grepl(pattern = "BarabasiProx", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][5]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
    names(upSample_model_param[[tmp]]) <- paste(prox_comp, names(upSample_model_param[[tmp]]), sep = "_")
  }
}
rm(tmp)

upSample_model_param <- unlist(upSample_model_param, recursive = FALSE)
upSample_model_param <- bind_rows(upSample_model_param, .id = "model")
upSample_model_param <- separate(upSample_model_param, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

keep <- colnames(upSample_model_param)[grep("BestTune_", colnames(upSample_model_param))]
upSample_model_param <- upSample_model_param[, c("featureType", "model", "Fold", keep)]
upSample_model_param$imbalance <- "upSample"





# Get the parameters for downSample models
files <- list.files(path = paste0("Analysis/STRING/DrugCombs_v5/", disease),
                    pattern = "^models_downSample_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

downSample_model_param <- list()

for(file in files){
  sheet_names_downSample <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_downSample){
    tmp <- strsplit(x = file, split = "\\/")[[1]][5]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    downSample_model_param[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
  if(grepl(pattern = "BarabasiProx", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][5]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
    names(downSample_model_param[[tmp]]) <- paste(prox_comp, names(downSample_model_param[[tmp]]), sep = "_")
  }
}
rm(tmp)

downSample_model_param <- unlist(downSample_model_param, recursive = FALSE)
downSample_model_param <- bind_rows(downSample_model_param, .id = "model")
downSample_model_param <- separate(downSample_model_param, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

keep <- colnames(downSample_model_param)[grep("BestTune_", colnames(downSample_model_param))]
downSample_model_param <- downSample_model_param[, c("featureType", "model", "Fold", keep)]
downSample_model_param$imbalance <- "downSample"





# Merge all model stats and rearrange for plotting
model_param <- rbind(none_model_param, smote_model_param, upSample_model_param, downSample_model_param) 
keep <- colnames(model_param)[grep("BestTune_", colnames(model_param))]

model_param <- reshape(model_param, direction = "long",
                       varying = keep,
                       v.names = "value",
                       timevar = "parameter",
                       times = keep)

rownames(model_param) <- NULL
model_param <- model_param[!model_param$parameter %in% c("BestTune_usekernel", "BestTune_kernel"), ]
model_param$value <- as.numeric(model_param$value)
model_param <- na.exclude(model_param)

model_param$model <- factor(model_param$model, levels = c("glmnet", "rf", "svmRadial", "nb", "knn"))
model_param$parameter <- factor(model_param$parameter, levels = c("BestTune_alpha", "BestTune_lambda", 
                                                                  "BestTune_mtry", 
                                                                  "BestTune_C", "BestTune_sigma",
                                                                  "BestTune_adjust", "BestTune_laplace", "BestTune_usekernel",
                                                                  "BestTune_kmax", "BestTune_distance", "BestTune_kernel"))
model_param$featureType <- factor(x = model_param$featureType,
                                            levels = c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene",
                                                        "BbsiProx_separation", "SteinerTopol", 
                                                        "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet"))




# Plot the results
svglite(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/ModelParameters_", disease, ".svg"), width = 15, height = 3)
labels <- unique(model_param[, c("model", "parameter")])
labels$plot_label <- paste0(labels$parameter, " (", labels$model, ")")
plot_label <- labels$plot_label
names(plot_label) <- labels$parameter
rm(labels)
ggplot(model_param, aes(x = featureType, y = value, fill = imbalance)) + 
  geom_boxplot(width = 0.6, lwd = 0.1, outlier.shape = NA) +
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
  ggtitle(paste0("Model parameters for ", disease, " drug combinations")) +
  xlab("Feature type") + ylab("Parameters") + 
  labs(colour = "Imbalance : ") +
  facet_wrap(vars(parameter), nrow = 1, scales = "free_y", labeller = as_labeller(plot_label))
dev.off()

print(warnings())