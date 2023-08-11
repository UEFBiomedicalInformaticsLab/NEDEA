set.seed(5081)



# Script for plotting the distribution of NES of external data over training data for important features


# Load libraries
library(optparse)
library(tidyverse)
library(openxlsx)
library(ggfortify)
library(ggpubr)
library(patchwork)


# Get arguments
option_list = list(make_option(c("--disease"), type = "character", default = NULL, 
                               help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
                   make_option(c("--data_balance_method"), type = "character", default = "none",
                               help = "The method to be used to balance imbalanced data. Possible values: none. Default: none.", metavar = "character"),
                   make_option(c("--feature_type"), type = "character", default = NULL,
                               help = "The feature type to use for modelling. Possible values: Disease2Gene, WithdrawalAdr2Gene, CombinedDisAdr2Gene, keggPath, SMPDbPath_DrugMet, SMPDbPath_DrugAction, miscGeneSet. Default: NULL", metavar = "character"))


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}


if(!opt$data_balance_method %in% c("none")){
  stop("--data_balance_method: No data balancing methods currently included. Default: none")
}


if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
}


# Define global options for this script
data_balance_method <- opt$data_balance_method
disease <- opt$disease
feature_type <- opt$feature_type




# Read the training data feature matrix
switch(feature_type,
       "CombinedDisAdr2Gene" = {
         feature_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
         feature_matrix <- feature_matrix$NES_CombinedDisAdr2Gene
       }
)
trainData_feature_matrix <- feature_matrix
trainData_feature_matrix <- column_to_rownames(trainData_feature_matrix, "Term")
trainData_feature_matrix <- as.data.frame(t(trainData_feature_matrix))

# Read the external data feature matrix
extData_feature_matrix <- readRDS(paste0("OutputFiles/External_predictions/CDCDB_AACT/features/CDCDB_AACT__", disease, "_", feature_type, ".rds"))
extData_feature_matrix <- as.data.frame(t(extData_feature_matrix))



featureImp_xDis <- list()

for(model in c("glmnet", "rf", "nb", "svmRadial")){
  
  feature_type_tmp <- feature_type
  model_tmp <- model
  
  if(feature_type == "Disease2Gene"){feature_type_tmp <- "Dis2Gene"}
  if(feature_type == "WithdrawalAdr2Gene"){feature_type_tmp <- "WdrlAdr2Gene"}
  if(feature_type == "CombinedDisAdr2Gene"){feature_type_tmp <- "CombDisAdr2Gene"}
  if(model == "svmRadial"){model_tmp <- "svmRd"}
  
  
  featureImp_xDis[[model_tmp]] <- read.table(paste0("OutputFiles/Tables/featureImportance_xDis/", feature_type_tmp, "_", model_tmp, "_", data_balance_method, ".tsv"), 
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



plot_list <- list()

for(model in c("glmnet", "rf", "nb", "svmRadial")){
  
  print(model)
  
  predicted_class <- read.xlsx(paste0("OutputFiles/External_predictions/CDCDB_AACT/ExtPred_CDCDB_AACT__", disease, "_", model, "_", data_balance_method, "_", feature_type, ".xlsx"))
  extData_feature_matrix$Class <- predicted_class$predicted_class[match(row.names(extData_feature_matrix), predicted_class$drugCombs)]
  
  print(table(extData_feature_matrix$Class))
  
  
  # Plot with efficacy features
  featureImp_select  <- featureImp_xDis_Dis2gene
  
  plot_data_train <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  x_axis_name <- colnames(plot_data_train)[1]
  y_axis_name <- colnames(plot_data_train)[2]
  colnames(plot_data_train) <- c("F1", "F2")
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  plot_data_train$Class <- paste("Train", plot_data_train$Class, sep = "_")
  
  plot_data_ext <- extData_feature_matrix[,c(featureImp_select, "Class"), drop = FALSE]
  colnames(plot_data_ext) <- c("F1", "F2", "Class")
  plot_data_ext$Class <- paste("Ext", plot_data_ext$Class, sep = "_")
  
  
  plot_scatter_efficacy <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    scale_colour_manual(values = c("Train_Eff" = "#00ff00", "Train_Adv" = "#ffa500", "Ext_Eff" = "#0000ff", "Ext_Adv" = "#ff1493")) + 
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
    ggtitle(paste0("Efficacy + ", model)) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  
  
  
  # Plot with safety features
  featureImp_select  <- featureImp_xDis_WdrlAdr2Gene
  
  plot_data_train <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  x_axis_name <- colnames(plot_data_train)[1]
  y_axis_name <- colnames(plot_data_train)[2]
  colnames(plot_data_train) <- c("F1", "F2")
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  plot_data_train$Class <- paste("Train", plot_data_train$Class, sep = "_")
  
  plot_data_ext <- extData_feature_matrix[,c(featureImp_select, "Class"), drop = FALSE]
  colnames(plot_data_ext) <- c("F1", "F2", "Class")
  plot_data_ext$Class <- paste("Ext", plot_data_ext$Class, sep = "_")
  
  
  plot_scatter_Safety <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    scale_colour_manual(values = c("Train_Eff" = "#00ff00", "Train_Adv" = "#ffa500", "Ext_Eff" = "#0000ff", "Ext_Adv" = "#ff1493")) + 
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
    ggtitle(paste0("Safety  + ", model)) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  
  
  # Plot for efficacy and safety
  featureImp_select  <- c(featureImp_xDis_Dis2gene[1], featureImp_xDis_WdrlAdr2Gene[1])
  
  plot_data_train <- trainData_feature_matrix[,featureImp_select, drop = FALSE]
  x_axis_name <- colnames(plot_data_train)[1]
  y_axis_name <- colnames(plot_data_train)[2]
  colnames(plot_data_train) <- c("F1", "F2")
  plot_data_train$Class <- substr(row.names(plot_data_train), 1, 3)
  plot_data_train$Class <- paste("Train", plot_data_train$Class, sep = "_")
  
  plot_data_ext <- extData_feature_matrix[,c(featureImp_select, "Class"), drop = FALSE]
  colnames(plot_data_ext) <- c("F1", "F2", "Class")
  plot_data_ext$Class <- paste("Ext", plot_data_ext$Class, sep = "_")
  
  
  plot_scatter_efficacySafety <- ggplot() +
    geom_point(data = plot_data_train, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    geom_point(data = plot_data_ext, 
               aes(x = F1, y = F2, color = Class),
               size = 0.5,
               shape = 3) +
    scale_colour_manual(values = c("Train_Eff" = "#00ff00", "Train_Adv" = "#ffa500", "Ext_Eff" = "#0000ff", "Ext_Adv" = "#ff1493")) + 
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
    ggtitle(paste0("Efficacy + Safety + ", model)) +
    xlab(str_wrap(x_axis_name, 30)) +
    ylab(str_wrap(y_axis_name, 30))
  
  
  
  plot_list[[model]] <- plot_scatter_efficacy + plot_scatter_Safety + plot_scatter_efficacySafety + plot_layout(guides = "collect") & theme(legend.position = "none")
}



# Save the plot
if(!dir.exists("OutputFiles/Plots/ExtPred_NES/")){
  dir.create("OutputFiles/Plots/ExtPred_NES/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/ExtPred_NES/ExtPred_NES_", disease, "_", feature_type, "_", data_balance_method, ".tiff"),
     width = 30, height = 40,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = plot_list, 
          nrow = 4, ncol = 1, 
          common.legend = TRUE, 
          legend = "bottom")

dev.off()