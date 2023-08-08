set.seed(5081)



# Script to check the distribution of features from training set data (RWR-FGSEA)



# Load libraries
library(ggfortify)
library(ggpubr)
library(tidyverse)
library(patchwork)



plot_list <- list()


for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
 
  
  # For efficacy and safety features from RWR-FGSEA
  fgsea_result <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  
  for(featureType in names(fgsea_result)){
    
    # Read the training data feature matrix
    feature_matrix <- fgsea_result[[featureType]]
    feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
    feature_matrix <- as.data.frame(t(feature_matrix))
    
    
    #  Extract data for PCA
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
            plot.title = element_text(hjust = 0.5),
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
      xlab(paste0("PC1 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[1], "%)")) + 
      ylab(paste0("PC2 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[2], "%)"))
    
    
    variance_explained <- as.data.frame(summary(pca)$importance)
    variance_explained <- round(variance_explained, 2)
    variance_explained <- as.data.frame(t(variance_explained))
    variance_explained <- rownames_to_column(variance_explained, "PC")
    colnames(variance_explained) <- gsub(" ", "_", colnames(variance_explained))
    variance_explained$PC <- factor(x = variance_explained$PC, levels = variance_explained$PC)
    if(length(variance_explained$PC) > 10){variance_explained <- variance_explained[1:10,]}
    
    # Plot the variance explained by each PCs
    pca_bar <- ggplot(variance_explained, aes(x = PC, y = Proportion_of_Variance)) +
      geom_bar(stat = "identity") +
      ggtitle(disease) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 8), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2)) +
      ylab("Proportion of variance") 
    
    # Plot the cumulative variance explained by the PCs
    pca_step <- ggplot(variance_explained, aes(x = PC, y = Cumulative_Proportion)) +
      geom_step(aes(group = 1)) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 8), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2)) +
      ylab("Cumulative variance") 
    
    
    plot_list[[featureType]][[disease]] <- pca_scatter + pca_bar + pca_step
    
  }
  
  
  # For common libraries
  fgsea_result <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))
  
  for(featureType in names(fgsea_result)){
    
    # Read the training data feature matrix
    feature_matrix <- fgsea_result[[featureType]]
    feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
    feature_matrix <- as.data.frame(t(feature_matrix))
    
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
            plot.title = element_text(hjust = 0.5),
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
      xlab(paste0("PC1 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[1], "%)")) + 
      ylab(paste0("PC2 (",   round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)[2], "%)"))
    
    
    variance_explained <- as.data.frame(summary(pca)$importance)
    variance_explained <- round(variance_explained, 2)
    variance_explained <- as.data.frame(t(variance_explained))
    variance_explained <- rownames_to_column(variance_explained, "PC")
    colnames(variance_explained) <- gsub(" ", "_", colnames(variance_explained))
    variance_explained$PC <- factor(x = variance_explained$PC, levels = variance_explained$PC)
    if(length(variance_explained$PC) > 10){variance_explained <- variance_explained[1:10,]}
    
    # Plot the variance explained by each PCs
    pca_bar <- ggplot(variance_explained, aes(x = PC, y = Proportion_of_Variance)) +
      geom_bar(stat = "identity") +
      ggtitle(disease) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 8), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2)) +
      ylab("Proportion of variance") 
    
    # Plot the variance explained by the PCs
    pca_step <- ggplot(variance_explained, aes(x = PC, y = Cumulative_Proportion)) +
      geom_step(aes(group = 1)) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 8), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2)) +
      ylab("Cumulative variance") 
    
    
    plot_list[[featureType]][[disease]] <- pca_scatter + pca_bar + pca_step
    
  }
  
  
}



# Save to file
if(!dir.exists("OutputFiles/Plots/Feature_variance/")){
  dir.create("OutputFiles/Plots/Feature_variance/", recursive = TRUE)
}  
for(featureType in names(plot_list)){
  tiff(paste0("OutputFiles/Plots/Feature_variance/", featureType, ".tiff"),
       width = 40, height = 20,
       units = "cm", compression = "lzw", res = 1200)
  plot <- ggarrange(plotlist = plot_list[[featureType]], 
            nrow = 3, ncol = 2, 
            common.legend = TRUE, legend = "bottom")
  print(plot)
  dev.off()
}


print(warnings())