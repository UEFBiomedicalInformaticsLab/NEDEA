set.seed(5081)



# Applicability domain analysis using the combined efficacy and safety estimates from the RWR-FGSEA approach
# Checks for outliers


# Load libraries 
library(ggfortify)
library(ggpubr)
library(tidyverse)
library(FNN)



# For external drug combinations from C-DCDB (RX/OTC)
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB/features/CDCDB_OrangeBook_rxotc__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  
  # Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = FALSE, scale. = FALSE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Compute k-NN distances for training set
  k <- 10 # you can adjust this value
  knn_train <- get.knn(pca_res_train[,1:2], k = k)
  
  # Compute distribution of average k-NN distances for training set
  knn_train_avg <- rowMeans(knn_train$nn.dist)
  quants <- quantile(knn_train_avg, probs = c(0.25, 0.75))
  iqr <- IQR(knn_train_avg)
  
  # Combine train_set and test_set
  combined_set <- rbind(pca_res_train[,1:2], pca_res_ext[,1:2])
  
  # Compute k-NN distances for combined set
  knn_combined <- get.knn(combined_set, k = k)
  
  # Get k-NN distances for the test set from the combined set k-NN distances
  knn_ext <- knn_combined$nn.dist[(nrow(pca_res_train) + 1):nrow(combined_set),]
  
  # Compute average k-NN distances for test set
  knn_ext_avg <- rowMeans(knn_ext)
  
  # Check if test instances fall within the AD
  ad_ext <- (knn_ext_avg >= (quants[1] - 1.5 * iqr)) & (knn_ext_avg <= (quants[2] + 1.5 * iqr))
  
  # Add applicability domain info to the test set
  pca_res_ext <- data.frame(pca_res_ext, ad = ad_ext)
  
  # Visualize the training and test data in the reduced space
  
  plot_list[[disease]] <- ggplot() +
    geom_point(data = data.frame(pca_res_train), 
               aes(x = PC1, y = PC2), 
               color = "blue",
               size = 0.5 ,
               shape = 3,) +
    geom_point(data = pca_res_ext, 
               aes(x = PC1, y = PC2, color = ad), 
               size = 0.5 ,
               shape = 3) +
    scale_color_manual(values = c("red", "green"), labels = c("Outside AD", "Inside AD")) +
    labs(color = "Applicability Domain") +
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
    ggtitle(disease) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
}
if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/ADcheck_NES_CombinedDisAdr2Gene__CDCDB_OrangeBook_rxotc.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()




# For external drug combinations from C-DCDB (AACT)
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB_AACT/features/CDCDB_AACT__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  
  #  Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = FALSE, scale. = FALSE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Compute k-NN distances for training set
  k <- 10 # you can adjust this value
  knn_train <- get.knn(pca_res_train[,1:2], k = k)
  
  # Compute distribution of average k-NN distances for training set
  knn_train_avg <- rowMeans(knn_train$nn.dist)
  quants <- quantile(knn_train_avg, probs = c(0.25, 0.75))
  iqr <- IQR(knn_train_avg)
  
  # Combine train_set and test_set
  combined_set <- rbind(pca_res_train[,1:2], pca_res_ext[,1:2])
  
  # Compute k-NN distances for combined set
  knn_combined <- get.knn(combined_set, k = k)
  
  # Get k-NN distances for the test set from the combined set k-NN distances
  knn_ext <- knn_combined$nn.dist[(nrow(pca_res_train) + 1):nrow(combined_set),]
  
  # Compute average k-NN distances for test set
  knn_ext_avg <- rowMeans(knn_ext)
  
  # Check if test instances fall within the AD
  ad_ext <- (knn_ext_avg >= (quants[1] - 1.5 * iqr)) & (knn_ext_avg <= (quants[2] + 1.5 * iqr))
  
  # Add applicability domain info to the test set
  pca_res_ext <- data.frame(pca_res_ext, ad = ad_ext)
  
  # Visualize the training and test data in the reduced space
  
  plot_list[[disease]] <- ggplot() +
    geom_point(data = data.frame(pca_res_train), 
               aes(x = PC1, y = PC2), 
               color = "blue",
               size = 0.5 ,
               shape = 3,) +
    geom_point(data = pca_res_ext, 
               aes(x = PC1, y = PC2, color = ad), 
               size = 0.5 ,
               shape = 3) +
    scale_color_manual(values = c("red", "green"), labels = c("Outside AD", "Inside AD")) +
    labs(color = "Applicability Domain") +
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
    ggtitle(disease) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
}
if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/ADcheck_NES_CombinedDisAdr2Gene__CDCDB_AACT.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()




# For external drug combinations from CDD
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
  trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDD/features/CDD_drugCombs__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))
  
  # Extract data for PCA
  pca_data_train <- trainData_feature_matrix
  pca_data_ext <- extData_feature_matrix
  
  # Perform PCA on training data
  pca_train <- prcomp(x = pca_data_train, center = FALSE, scale. = FALSE)
  
  # Transform training and test set
  pca_res_train <- predict(pca_train, pca_data_train)
  pca_res_ext <- predict(pca_train, pca_data_ext)
  
  # Compute k-NN distances for training set
  k <- 10 # you can adjust this value
  knn_train <- get.knn(pca_res_train[,1:2], k = k)
  
  # Compute distribution of average k-NN distances for training set
  knn_train_avg <- rowMeans(knn_train$nn.dist)
  quants <- quantile(knn_train_avg, probs = c(0.25, 0.75))
  iqr <- IQR(knn_train_avg)
  
  # Combine train_set and test_set
  combined_set <- rbind(pca_res_train[,1:2], pca_res_ext[,1:2])
  
  # Compute k-NN distances for combined set
  knn_combined <- get.knn(combined_set, k = k)
  
  # Get k-NN distances for the test set from the combined set k-NN distances
  knn_ext <- knn_combined$nn.dist[(nrow(pca_res_train) + 1):nrow(combined_set),]
  
  # Compute average k-NN distances for test set
  knn_ext_avg <- rowMeans(knn_ext)
  
  # Check if test instances fall within the AD
  ad_ext <- (knn_ext_avg >= (quants[1] - 1.5 * iqr)) & (knn_ext_avg <= (quants[2] + 1.5 * iqr))
  
  # Add applicability domain info to the test set
  pca_res_ext <- data.frame(pca_res_ext, ad = ad_ext)
  
  # Visualize the training and test data in the reduced space
  
  plot_list[[disease]] <- ggplot() +
    geom_point(data = data.frame(pca_res_train), 
               aes(x = PC1, y = PC2), 
               color = "blue",
               size = 0.5 ,
               shape = 3) +
    geom_point(data = pca_res_ext, 
               aes(x = PC1, y = PC2, color = ad), 
               size = 0.5 ,
               shape = 3) +
    scale_color_manual(values = c("red", "green"), labels = c("Outside AD", "Inside AD")) +
    labs(color = "Applicability Domain") +
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
    ggtitle(disease) +
    xlab(paste0("PC1 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[1], "%)")) + 
    ylab(paste0("PC2 (",   round(pca_train$sdev^2 / sum(pca_train$sdev^2) * 100, 2)[2], "%)"))
}
if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/ADcheck_NES_CombinedDisAdr2Gene__CDD.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()



print(warnings())