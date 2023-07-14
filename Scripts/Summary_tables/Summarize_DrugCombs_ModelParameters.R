set.seed(5081)



# Script for tabulate the model parameters across different folds



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
keep_cols <- colnames(model_param)[grep("BestTune_", colnames(model_param))]
keep_rows <- c("Dis2Gene", "WdrlAdr2Gene", "CombDisAdr2Gene", 
               "DrugDisease_BbsiProx_separation", "DrugAdr_BbsiProx_separation", "DrgDisAdr_BbsiProx_separation", 
               "keggPath", "SMPDbPath_DrugMet", "SMPDbPath_DrugAction", "miscGeneSet")

model_param <- model_param[model_param$featureType %in% keep_rows, ]
model_param <- model_param[model_param$imbalance == "none", ]
model_param$featureType <- factor(model_param$featureType, levels = keep_rows)
row.names(model_param) <- NULL
model_param <- split(model_param, f = model_param$featureType)

if(!dir.exists(paste0("OutputFiles/Tables/", disease, "/"))){
  dir.create(paste0("OutputFiles/Tables/", disease, "/"), recursive = TRUE)
}
write.xlsx(model_param, paste0("OutputFiles/Tables/", disease, "/ModelParameters_", disease, ".xlsx"), overwrite = TRUE)



print(warnings())