set.seed(5081)



# Applicability domain analysis using the combined efficacy and safety estimates from the RWR-FGSEA approach


# Load libraries 
library(ggfortify)
library(ggpubr)
library(tidyverse)


featureType <- "CombDisAdr2Gene"
model <- "nb"
featureImp_xDis <- read.table(paste0("OutputFiles/Tables/featureImportance_xDis/", featureType, "_", model, "_none.tsv"), sep = "\t", header = TRUE)


# For external drug combinations from C-DCDB (RX/OTC)
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Select the top two important features
  featureImp_select <- featureImp_xDis[, grep(disease, colnames(featureImp_xDis))]
  featureImp_select <- featureImp_select[1:2, 1]
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  trainData_feature_matrix <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB/features/CDCDB_OrangeBook_rxotc__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  extData_feature_matrix <- extData_feature_matrix[,featureImp_select, drop = FALSE]
  
  # Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Remove columns with zero variance
  keep <- names(which(apply(pca_data_train, 2, var) != 0))
  pca_data_train <- pca_data_train[, keep, drop = FALSE]
  pca_data_ext <- pca_data_ext[, keep, drop = FALSE]
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = TRUE, scale. = TRUE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Plot
  plot_data_train <- data.frame(pca_res_train)
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  
  plot_data_ext <- data.frame(pca_res_ext)
  plot_data_ext$Class <- substr(row.names(plot_data_ext), 1, 3)
  
  pca_scatter <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = PC1, y = PC2, color = Class), 
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext,
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
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
    ggtitle(paste0(featureImp_select, collapse = "\n")) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
  
  plot_list[[disease]] <- pca_scatter
}
if(!dir.exists("OutputFiles/Plots/Applicibility_domain_impFeature/")){
  dir.create("OutputFiles/Plots/Applicibility_domain_impFeature/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain_impFeature/AD_NES_CombinedDisAdr2Gene_", model, "__CDCDB_OrangeBook_rxotc.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()



# For external drug combinations from C-DCDB (AACT)
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Select the top two important features
  featureImp_select <- featureImp_xDis[, grep(disease, colnames(featureImp_xDis))]
  featureImp_select <- featureImp_select[1:2, 1]
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  trainData_feature_matrix <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB_AACT/features/CDCDB_AACT__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  extData_feature_matrix <- extData_feature_matrix[,featureImp_select, drop = FALSE]
  
  #  Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Remove columns with zero variance
  keep <- names(which(apply(pca_data_train, 2, var) != 0))
  pca_data_train <- pca_data_train[, keep, drop = FALSE]
  pca_data_ext <- pca_data_ext[, keep, drop = FALSE]
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = TRUE, scale. = TRUE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Plot
  plot_data_train <- data.frame(pca_res_train)
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  
  plot_data_ext <- data.frame(pca_res_ext)
  plot_data_ext$Class <- substr(row.names(plot_data_ext), 1, 3)
  plot_data_ext$Class <- gsub("^NCT", "Unk", plot_data_ext$Class)
  
  pca_scatter <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = PC1, y = PC2, color = Class), 
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext,
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
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
    ggtitle(paste0(featureImp_select, collapse = "\n")) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
  
  plot_list[[disease]] <- pca_scatter
}

if(!dir.exists("OutputFiles/Plots/Applicibility_domain_impFeature/")){
  dir.create("OutputFiles/Plots/Applicibility_domain_impFeature/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain_impFeature/AD_NES_CombinedDisAdr2Gene_", model, "__CDCDB_AACT.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()




# For external drug combinations from CDD
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Select the top two important features
  featureImp_select <- featureImp_xDis[, grep(disease, colnames(featureImp_xDis))]
  featureImp_select <- featureImp_select[1:2, 1]
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  trainData_feature_matrix <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDD/features/CDD_drugCombs__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  extData_feature_matrix <- extData_feature_matrix[,featureImp_select, drop = FALSE]
  
  # Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Remove columns with zero variance
  keep <- names(which(apply(pca_data_train, 2, var) != 0))
  pca_data_train <- pca_data_train[, keep, drop = FALSE]
  pca_data_ext <- pca_data_ext[, keep, drop = FALSE]
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = TRUE, scale. = TRUE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Plot
  plot_data_train <- data.frame(pca_res_train)
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  
  plot_data_ext <- data.frame(pca_res_ext)
  plot_data_ext$Class <- substr(row.names(plot_data_ext), 1, 3)
  
  pca_scatter <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = PC1, y = PC2, color = Class), 
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext,
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
          legend.text = element_text(size = 4),
          legend.margin = margin(1,1,1,1),
          legend.box.spacing = unit(0.1, 'cm'),
          legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
    ggtitle(paste0(featureImp_select, collapse = "\n")) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
  
  plot_list[[disease]] <- pca_scatter
}

if(!dir.exists("OutputFiles/Plots/Applicibility_domain_impFeature/")){
  dir.create("OutputFiles/Plots/Applicibility_domain_impFeature/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain_impFeature/AD_NES_CombinedDisAdr2Gene_", model, "__CDD.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()


print(warnings())