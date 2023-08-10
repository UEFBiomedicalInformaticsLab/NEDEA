
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(patchwork)

plot_list <- list()


# disease <- "BreastCancer"
data_balance_method <- "none"


for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  print(disease)
  
  featureImp_xDis <- list()
  
  for(model in c("glmnet", "rf", "nb", "svmRadial")){
    model_tmp <- model
    if(model == "svmRadial"){model_tmp <- "svmRd"}
    
    
    featureImp_xDis[[model_tmp]] <- read.table(paste0("OutputFiles/Tables/featureImportance_xDis/CombDisAdr2Gene_", model_tmp, "_", data_balance_method, ".tsv"), 
                                               sep = "\t", 
                                               header = TRUE)
  }
  
  
  featureImp_xDis <- lapply(featureImp_xDis, function(x){x[, grep(disease, colnames(x))]})
  
  
  featureImp_xDis_Dis2gene <- lapply(featureImp_xDis, function(x){x[grep("\\[DISEASE\\]", x[, c(1)]),]})
  featureImp_xDis_Dis2gene <- lapply(featureImp_xDis_Dis2gene, function(x){x[1:3, grep("TopK_median", colnames(x))]})
  featureImp_xDis_Dis2gene <- sort(table(unlist(featureImp_xDis_Dis2gene)), decreasing = TRUE)
  print(featureImp_xDis_Dis2gene)
  featureImp_xDis_Dis2gene <- names(featureImp_xDis_Dis2gene)[1:2]
  
  featureImp_xDis_WdrlAdr2Gene <- lapply(featureImp_xDis, function(x){x[grep("\\[ADR\\]", x[, c(1)]),]})
  featureImp_xDis_WdrlAdr2Gene <- lapply(featureImp_xDis_WdrlAdr2Gene, function(x){x[1:3, grep("TopK_median", colnames(x))]})
  featureImp_xDis_WdrlAdr2Gene <- sort(table(unlist(featureImp_xDis_WdrlAdr2Gene)), decreasing = TRUE)
  print(featureImp_xDis_WdrlAdr2Gene)
  featureImp_xDis_WdrlAdr2Gene <- names(featureImp_xDis_WdrlAdr2Gene)[1:2]
  
  featureImp_xDis_WdrlAdr2Gene <- c("[ADR] Chemical and Drug Induced Liver Injury (MESH:D056486) [CTD]", "[ADR] Drug toxicity (C0013221) [DisGeNET]") # Hard coding to test
  
  
  # Plot for efficacy
  featureImp_select  <- featureImp_xDis_Dis2gene

  feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
  feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
  feature_matrix <- as.data.frame(t(feature_matrix))
  feature_matrix <- feature_matrix[,featureImp_select, drop = FALSE]
  
  
  plot_data <- as.data.frame(feature_matrix)
  x_axis_name <- colnames(plot_data)[1]
  y_axis_name <- colnames(plot_data)[2]
  colnames(plot_data) <- c("F1", "F2")
  plot_data$Class <- substr(row.names(plot_data), 1, 3)
  
  plot_scatter_efficacy <- ggplot() +
    geom_point(data = plot_data, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(margin = margin(1,1,1,1)),
          text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 8),
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
    ggtitle(disease) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  # Plot for safety
  featureImp_select  <- featureImp_xDis_WdrlAdr2Gene

  feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
  feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
  feature_matrix <- as.data.frame(t(feature_matrix))
  feature_matrix <- feature_matrix[,featureImp_select, drop = FALSE]
  
  
  plot_data <- as.data.frame(feature_matrix)
  x_axis_name <- colnames(plot_data)[1]
  y_axis_name <- colnames(plot_data)[2]
  colnames(plot_data) <- c("F1", "F2")
  plot_data$Class <- substr(row.names(plot_data), 1, 3)
  
  plot_scatter_safety <- ggplot() +
    geom_point(data = plot_data, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(margin = margin(1,1,1,1)),
          text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 8),
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
    ggtitle(disease) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  
  # Plot for efficacy and safety
  featureImp_select  <- c(featureImp_xDis_Dis2gene[1], featureImp_xDis_WdrlAdr2Gene[1])
  
  feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
  feature_matrix <- column_to_rownames(feature_matrix, colnames(feature_matrix)[1])
  feature_matrix <- as.data.frame(t(feature_matrix))
  feature_matrix <- feature_matrix[,featureImp_select, drop = FALSE]
  
  
  plot_data <- as.data.frame(feature_matrix)
  x_axis_name <- colnames(plot_data)[1]
  y_axis_name <- colnames(plot_data)[2]
  colnames(plot_data) <- c("F1", "F2")
  plot_data$Class <- substr(row.names(plot_data), 1, 3)
  
  plot_scatter_efficacySafety <- ggplot() +
    geom_point(data = plot_data, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_rect(color = "black", linewidth = 0.25,),
          strip.text = element_text(margin = margin(1,1,1,1)),
          text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 8),
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
    ggtitle(disease) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  plot_list[[disease]] <- plot_scatter_efficacy + plot_scatter_safety + plot_scatter_efficacySafety
}

tiff(paste0("OutputFiles/Plots/SampleDist_NES_topFeature_byDisease.tiff"),
     width = 50, height = 30,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list,
          nrow = 3, ncol = 2)
dev.off()