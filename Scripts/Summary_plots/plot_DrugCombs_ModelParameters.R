set.seed(5081)



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
files <- list.files(path = paste0("OutputFiles/Model_train/", disease),
                    pattern = "^models_none_[a-zA-Z_]+.xlsx", 
                    ignore.case = TRUE, full.names = TRUE)

none_model_param <- list()

for(file in files){
  sheet_names_none <- getSheetNames(file)
  
  ## Read files
  for(name in sheet_names_none){
    tmp <- strsplit(x = file, split = "\\/")[[1]][4]
    tmp <- strsplit(x = tmp, split = "\\.")[[1]][1]
    none_model_param[[tmp]][[name]] <- read.xlsx(file, sheet = name)
  }
  if(grepl(pattern = "BarabasiProx", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][4]
    names(none_model_param[[tmp]]) <- paste(prox_comp, names(none_model_param[[tmp]]), sep = "_")
  }
  if(grepl(pattern = "BarabasiProx_DrgDisAdr", x = file)){
    prox_comp <- strsplit(x = file, split = "\\/")[[1]][4]
    prox_comp <- strsplit(x = prox_comp, split = "\\_")[[1]][6]
    prox_comp <- strsplit(x = prox_comp, split = "\\.")[[1]][1]
    names(none_model_param[[tmp]]) <- gsub(pattern = "\\.", 
                                           replacement = paste0("_", prox_comp, "."), 
                                           x = names(none_model_param[[tmp]]), )
  }
}
rm(tmp)

none_model_param <- unlist(none_model_param, recursive = FALSE)
none_model_param <- bind_rows(none_model_param, .id = "model")
none_model_param <- separate(none_model_param, col = "model", into = c("file", "featureType", "model"), sep = "\\.")

keep <- colnames(none_model_param)[grep("BestTune_", colnames(none_model_param))]
none_model_param <- none_model_param[, c("featureType", "model", "Fold", keep)]
none_model_param$imbalance <- "none"



# Merge all model stats and rearrange for plotting
model_param <- none_model_param 
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

model_param$model <- factor(model_param$model, levels = c("glmnet", "rf", "svmRadial", "nb"))
model_param$parameter <- factor(model_param$parameter, levels = c("BestTune_alpha", "BestTune_lambda", 
                                                                  "BestTune_mtry", 
                                                                  "BestTune_C", "BestTune_sigma",
                                                                  "BestTune_adjust", "BestTune_laplace"))

model_param$featureType <- factor(x = model_param$featureType,
                                  levels = unique(model_param$featureType))




# Plot the results

features_to_plot <- c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene",
                      "DrugDrug_BbsiProx_separation", "DrugDisease_BbsiProx_separation",
                      "DrugAdr_BbsiProx_separation", "DrgDisAdr_BbsiProx_separation",
                      "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet", "miscGeneSet")

select_model_param <- model_param[(model_param$featureType %in% features_to_plot),]

select_model_param$featureType <- gsub("Dis2Gene", "Efficacy", select_model_param$featureType)
select_model_param$featureType <- gsub("WdrlAdr2Gene", "Safety", select_model_param$featureType)
select_model_param$featureType <- gsub("CombDisAdr2Gene", "Combined-efficacy-safety", select_model_param$featureType)
select_model_param$featureType <- gsub("DrugDrug_BbsiProx_separation", "Drug-Drug", select_model_param$featureType)
select_model_param$featureType <- gsub("DrugAdr_BbsiProx_separation", "Drug-ADR", select_model_param$featureType)
select_model_param$featureType <- gsub("DrugDisease_BbsiProx_separation", "Drug-Disease", select_model_param$featureType)
select_model_param$featureType <- gsub("DrgDisAdr_BbsiProx_separation", "Drug-Disease-ADR", select_model_param$featureType)
select_model_param$featureType <- gsub("keggPath", "KEGG pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("SMPDbPath_DrugAction", "Drug action pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("SMPDbPath_DrugMet", "Drug metabolism pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("miscGeneSet", "Misc. gene sets", select_model_param$featureType)

labels <- unique(select_model_param[, c("model", "parameter")])
labels$plot_label <- paste0(gsub("^BestTune_", "", labels$parameter), " (", labels$model, ")")
plot_label <- labels$plot_label
names(plot_label) <- labels$parameter
rm(labels)


if(!dir.exists(paste0("OutputFiles/Plots/", disease))){
  dir.create(paste0("OutputFiles/Plots/", disease), recursive = TRUE)
}

svglite(paste0("OutputFiles/Plots/", disease, "/ModelParameters_", disease, ".svg"), 
        width = length(unique(model_param$featureType)), height = 3)

ggplot(select_model_param, aes(x = featureType, y = value, fill = imbalance)) + 
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