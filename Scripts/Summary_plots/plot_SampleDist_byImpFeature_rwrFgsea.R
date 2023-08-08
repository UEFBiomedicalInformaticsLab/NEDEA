set.seed(5081)



# Script for plotting the distribution of samples based on the top two features selected by the model


# Load libraries
library(optparse)
library(tidyverse)
library(ggfortify)
library(ggpubr)


# Get arguments
option_list = list(
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: none. Default: none.", metavar = "character"),
  make_option(c("--feature_type"), type = "character", default = NULL, 
              help = "The feature type to use for modelling. Possible values: Disease2Gene, WithdrawalAdr2Gene, CombinedDisAdr2Gene, keggPath, SMPDbPath_DrugMet, SMPDbPath_DrugAction, miscGeneSet. Default: NULL", metavar = "character"),
  make_option(c("--model"), type = "character", default = NULL, 
              help = "The modelling technique to use. Possible values: glmnet, nb, rf, svmRadial. Default: NULL", metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if(!opt$data_balance_method %in% c("none")){
  stop("--data_balance_method: No data balancing methods currently included. Default: none")
}

if(is.null(opt$model)){
  print_help(opt_parser)
  stop("--model argument needed", call.=FALSE)
}


if(!opt$model %in% c("glmnet", "nb", "rf", "svmRadial")){
  stop("--model must be glmnet, nb, rf, svmRadial")
}


if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
}


# Define global options for this script 
data_balance_method <- opt$data_balance_method
model <- opt$model
feature_type <- opt$feature_type

# Read the feature importance file
feature_type_tmp <- feature_type
model_tmp <- model
if(feature_type == "Disease2Gene"){feature_type_tmp <- "Dis2Gene"}
if(feature_type == "WithdrawalAdr2Gene"){feature_type_tmp <- "WdrlAdr2Gene"}
if(feature_type == "CombinedDisAdr2Gene"){feature_type_tmp <- "CombDisAdr2Gene"}
if(model == "svmRadial"){model_tmp <- "svmRd"}

featureImp_xDis <- read.table(paste0("OutputFiles/Tables/featureImportance_xDis/", feature_type_tmp, "_", model_tmp, "_", data_balance_method, ".tsv"), sep = "\t", header = TRUE)


plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  # Select the top two important features
  featureImp_select <- featureImp_xDis[, grep(disease, colnames(featureImp_xDis))]
  featureImp_select <- featureImp_select[1:2, 1]
  
  
  # Read the features
  switch(feature_type,
         "Disease2Gene" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_Disease2Gene
         }, 
         "WithdrawalAdr2Gene" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_WithdrawalAdr2Gene
         }, 
         "CombinedDisAdr2Gene" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
         }, 
         "keggPath" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_keggPath
         }, 
         "SMPDbPath_DrugMet" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_SMPDbPath_DrugMet
         }, 
         "SMPDbPath_DrugAction" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_SMPDbPath_DrugAction
         }, 
         "miscGeneSet" = {
           feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
           feature_matrix <- feature_matrix$NES_miscGeneSet
         }
  )
  
  feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
  feature_matrix <- as.data.frame(t(feature_matrix))
  feature_matrix <- feature_matrix[,featureImp_select, drop = FALSE]
  
  
  # Extract data for PCA
  pca_data <- feature_matrix
  
  # Remove columns with zero variance
  keep <- names(which(apply(pca_data, 2, var) != 0))
  pca_data <- pca_data[, keep]
  
  # Perform PCA on training data
  pca <- prcomp(x = pca_data, center = TRUE, scale. = TRUE)
  
  # Transform training and test set
  pca_res <- predict(pca, pca_data)
  
  plot_data <- data.frame(pca_res)
  plot_data$Class <- substr(row.names(plot_data), 1, 3)
  
  # Plot PCA scatter plot
  pca_scatter <- ggplot() +
    geom_point(data = plot_data, 
               aes(x = PC1, y = PC2, color = Class), 
               size = 0.5,
               shape = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(margin = margin(1,1,1,1)),
          text = element_text(size = 8), 
          plot.title = element_text(hjust = 0.5, size = 4),
          axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
    ggtitle(paste0(c(disease, featureImp_select), collapse = "\n")) +
    xlab(paste0("PC1 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[2], "%)"))
  
  plot_list[[disease]] <- pca_scatter
}
  

# Save the plot
if(!dir.exists("OutputFiles/Plots/SampleDist_byImpFeature/")){
  dir.create("OutputFiles/Plots/SampleDist_byImpFeature/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/SampleDist_byImpFeature/SampleDist_byImp_", feature_type, "_", model, "_", data_balance_method, ".tiff"),
     width = 20, height = 25,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()



print(warnings())