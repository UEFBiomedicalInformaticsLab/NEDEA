set.seed(5081)



# Applicability domain analysis using the combined efficacy and safety estimates from the RWR-FGSEA approach


# Load libraries 
library(ggfortify)
library(ggpubr)
library(tidyverse)


# For external drug combinations from C-DCDB (RX/OTC)
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  # Read the training data feature matrix
  trainData_feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix <- trainData_feature_matrix$NES_CombinedDisAdr2Gene
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB/features/CDCDB_OrangeBook_rxotc__", disease, "_CombinedDisAdr2Gene.rds"))
  
  # Merge the data
  merged_feature_matrix <- merge(trainData_feature_matrix, extData_feature_matrix, by.x = "Term", by.y = 0, all = TRUE)
  merged_feature_matrix <- column_to_rownames(merged_feature_matrix, "Term")
  merged_feature_matrix <- as.data.frame(t(merged_feature_matrix))
  merged_feature_matrix$Class <- factor(substr(row.names(merged_feature_matrix), 1, 3))
  
  # PCA
  pca_data <- merged_feature_matrix[, !colnames(merged_feature_matrix) %in% "Class"]
  pca_res <- prcomp(x = pca_data, center = FALSE, scale. = FALSE)
  pca_scatter <- autoplot(pca_res, 
                          data = merged_feature_matrix, 
                          colour = "Class", 
                          size = 0.5 ,
                          shape = 3,
                          main = disease)  + 
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
          legend.box.background = element_rect(colour = "black", size = 0.25))
  plot_list[[disease]] <- pca_scatter
}
if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/AD_NES_CombinedDisAdr2Gene__CDCDB_OrangeBook_rxotc.tiff"),
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
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB_AACT/features/CDCDB_AACT__", disease, "_CombinedDisAdr2Gene.rds"))
  
  # Merge the data
  merged_feature_matrix <- merge(trainData_feature_matrix, extData_feature_matrix, by.x = "Term", by.y = 0, all = TRUE)
  merged_feature_matrix <- column_to_rownames(merged_feature_matrix, "Term")
  merged_feature_matrix <- as.data.frame(t(merged_feature_matrix))
  merged_feature_matrix$Class <- factor(substr(row.names(merged_feature_matrix), 1, 3))
  merged_feature_matrix$Class <- gsub("^NCT", "Unk", merged_feature_matrix$Class)
  
  # PCA
  pca_data <- merged_feature_matrix[, !colnames(merged_feature_matrix) %in% "Class"]
  pca_res <- prcomp(x = pca_data, center = FALSE, scale. = FALSE)
  pca_scatter <- autoplot(pca_res, 
                          data = merged_feature_matrix, 
                          colour = "Class", 
                          size = 0.5 ,
                          shape = 3,
                          main = disease)  + 
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
          legend.box.background = element_rect(colour = "black", size = 0.25))
  plot_list[[disease]] <- pca_scatter
}

if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/AD_NES_CombinedDisAdr2Gene__CDCDB_AACT.tiff"),
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
  
  # Read the external data feature matrix
  extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDD/features/CDD_drugCombs__", disease, "_CombinedDisAdr2Gene.rds"))
  
  # Merge the data
  merged_feature_matrix <- merge(trainData_feature_matrix, extData_feature_matrix, by.x = "Term", by.y = 0, all = TRUE)
  merged_feature_matrix <- column_to_rownames(merged_feature_matrix, "Term")
  merged_feature_matrix <- as.data.frame(t(merged_feature_matrix))
  merged_feature_matrix$Class <- factor(substr(row.names(merged_feature_matrix), 1, 3))
  merged_feature_matrix$Class <- gsub("^NCT", "Unk", merged_feature_matrix$Class)
  
  # PCA
  pca_data <- merged_feature_matrix[, !colnames(merged_feature_matrix) %in% "Class"]
  pca_res <- prcomp(x = pca_data, center = FALSE, scale. = FALSE)
  pca_scatter <- autoplot(pca_res, 
                          data = merged_feature_matrix, 
                          colour = "Class", 
                          size = 0.5 ,
                          shape = 3,
                          main = disease)  + 
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
          legend.box.background = element_rect(colour = "black", size = 0.25))
  plot_list[[disease]] <- pca_scatter
}

if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/AD_NES_CombinedDisAdr2Gene__CDD.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()


print(warnings())