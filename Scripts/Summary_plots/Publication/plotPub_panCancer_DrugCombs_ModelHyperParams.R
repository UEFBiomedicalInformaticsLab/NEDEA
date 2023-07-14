set.seed(5081) 



# Script for plotting the model parameters (for publication: comparision of combined-efficacy-safety vs pthways vs misc gene set)



# Load libraries
library(openxlsx)
library(tidyverse)
library(svglite)
library(RColorBrewer)



model_param <- data.frame()
for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  cat(paste0("\nReading data for: ", disease, "\n"))
  
  # Get the statistics for unbalanced models
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
  
  tmp <- names(none_model_param)[grep("^BestTune_", names(none_model_param))]
  none_model_param <- none_model_param[, c("featureType", "model", "Fold", tmp)]
  none_model_param <- none_model_param[none_model_param$model %in% c("glmnet", "nb", "rf", "svmRadial"), ]
  none_model_param$imbalance <- "none"
  rm(tmp)
  
  
  
  # Merge all model stats and rearrange for plotting
  tmp <- none_model_param
  tmp$disease <- disease
  model_param <- rbind(model_param, tmp)
  
}


if(!dir.exists("OutputFiles/Plots/Publication/")){
  dir.create("OutputFiles/Plots/Publication/", recursive = TRUE)
}

saveRDS(model_param, "OutputFiles/Plots/Publication/model_param.rds")
# model_param <- readRDS("OutputFiles/Plots/Publication/model_param.rds")

model_param <- model_param[, !colnames(model_param) %in% c("BestTune_usekernel")]
tmp <- names(model_param)[grep("^BestTune_", names(model_param))]
model_param <- reshape(model_param, direction = "long",
                       varying = tmp,
                       v.names = "value",
                       timevar = "HyperParam",
                       times = tmp)

rownames(model_param) <- NULL
model_param$value <- as.numeric(model_param$value)
model_param <- na.exclude(model_param)
model_param$model <- factor(model_param$model, levels = c("glmnet", "rf", "svmRadial", "nb"))
model_param$HyperParam <- factor(model_param$HyperParam, levels = c("BestTune_alpha", "BestTune_lambda", 
                                                                  "BestTune_mtry", 
                                                                  "BestTune_C", "BestTune_sigma",
                                                                  "BestTune_adjust", "BestTune_laplace"))
rm(tmp)



# Plot the model test accuracy for selected features
features_to_plot <- unique(model_param$featureType)
features_to_plot <- features_to_plot[grep("_BbsiProx_separation$", features_to_plot)]
features_to_plot <- features_to_plot[grep("^DrugDrug_", features_to_plot, invert = TRUE)]
features_to_plot <- c("CombDisAdr2Gene", features_to_plot)
features_to_plot <- c("Dis2Gene", "WdrlAdr2Gene", features_to_plot, "keggPath", "SMPDbPath_DrugAction", "SMPDbPath_DrugMet", "miscGeneSet")

select_model_param <- model_param[(model_param$featureType %in% features_to_plot),]

select_model_param$featureType <- gsub("Dis2Gene", "Efficacy", select_model_param$featureType)
select_model_param$featureType <- gsub("WdrlAdr2Gene", "Safety", select_model_param$featureType)
select_model_param$featureType <- gsub("CombDisAdr2Gene", "Combined-efficacy-safety", select_model_param$featureType)
select_model_param$featureType <- gsub("DrugAdr_BbsiProx_separation", "Drug-ADR", select_model_param$featureType)
select_model_param$featureType <- gsub("DrugDisease_BbsiProx_separation", "Drug-Disease", select_model_param$featureType)
select_model_param$featureType <- gsub("DrgDisAdr_BbsiProx_separation", "Drug-Disease-ADR", select_model_param$featureType)
select_model_param$featureType <- gsub("keggPath", "KEGG pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("SMPDbPath_DrugAction", "Drug action pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("SMPDbPath_DrugMet", "Drug metabolism pathways", select_model_param$featureType)
select_model_param$featureType <- gsub("miscGeneSet", "Misc. gene sets", select_model_param$featureType)

select_model_param$disease <- gsub("Cancer$", "", select_model_param$disease)


select_model_param$featureType <- factor(x = select_model_param$featureType,
                                         levels = c("Efficacy", "Safety", "Combined-efficacy-safety", 
                                                    "Drug-Disease", "Drug-ADR", "Drug-Disease-ADR", 
                                                    "KEGG pathways", "Drug action pathways", "Drug metabolism pathways", 
                                                    "Misc. gene sets"))
select_model_param$imbalance <- factor(x = select_model_param$imbalance,
                                       levels = c("none"))

select_model_param$HyperParam_label <- paste0(gsub("^BestTune_", "", select_model_param$HyperParam), " (", select_model_param$model, ")")

# select_model_param <- na.exclude(select_model_param) # Some scores are NA if while calculating ratio the denominator is 0







tiff("OutputFiles/Plots/Publication/panCancer_ModelHyperparameters.tiff", 
     width = length(unique(select_model_param$featureType)) * 4, 
     height = length(unique(select_model_param$HyperParam)) * 4, 
     units = "cm", compression = "lzw", res = 1200)


ggplot(select_model_param, aes(x = disease, y = value)) + #
  geom_violin() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", size = 0.25,),
        strip.text = element_text(margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", size = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 4),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", size = 0.25)) +
  scale_fill_manual(values = brewer.pal(6, "Set2")) + 
  xlab("Cancers") + ylab("Accuracy scores") + #
  labs(colour = "Models : ") +
  facet_grid(rows = vars(HyperParam_label), cols = vars(featureType), scales = "free_y")
dev.off()



print(warnings())