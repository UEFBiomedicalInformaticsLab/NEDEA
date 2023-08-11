set.seed(5081)



# Script for plotting the distribution of enrichment scores


# Load libraries
library(optparse)
library(tidyverse)
library(patchwork)
library(ggpubr)


# Get arguments
option_list = list(make_option(c("--disease"), type = "character", default = NULL, 
                               help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
                   make_option(c("--feature_type"), type = "character", default = NULL,
                               help = "The feature type to use for modelling. Possible values: Disease2Gene, WithdrawalAdr2Gene, CombinedDisAdr2Gene, keggPath, SMPDbPath_DrugMet, SMPDbPath_DrugAction, miscGeneSet. Default: NULL", metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
}


# Define global options for this script
disease <- opt$disease
feature_type <- opt$feature_type



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
select_cols <- colnames(feature_matrix)
plot_data <- rownames_to_column(feature_matrix, "drugComb")


plot_data <- reshape(plot_data,
                     varying = select_cols,
                     v.names = "NES",
                     timevar = "Terms",
                     times = select_cols,
                     direction = "long")
row.names(plot_data) <- NULL
plot_data$Class <- substr(plot_data$drugComb, 1,3)


# Plot the complete NES
plot_box1 <- ggboxplot(data = plot_data,
                       x = "Terms", y = "NES", 
                       width = 0.5, 
                       lwd = 0.1,
                       fill = "grey",
                       outlier.shape = NA) +
					   stat_summary(fun = "mean", 
					   color = "red", 
					   size = 0.2) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        strip.background = element_rect(color = "black", linewidth = 0.1,),
        text = element_text(size = 4), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) + 
  coord_flip() 

# Separate the two classes and plot
plot_box2 <- ggboxplot(data = plot_data,
                       x = "Terms", y = "NES", 
                       fill = "Class",
                       width = 0.5, 
                       lwd = 0.1,
                       outlier.shape = NA) +
  stat_compare_means(aes(group = Class), 
                     method = "t.test", 
                     label = "p.signif", 
					 size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        strip.background = element_rect(color = "black", linewidth = 0.1,),
        text = element_text(size = 4), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_blank()) + 
  ylab("") + 
  coord_flip() 


# Save the plot
if(!dir.exists("OutputFiles/Plots/NES_distribution/")){
  dir.create("OutputFiles/Plots/NES_distribution/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/NES_distribution/NES_distribution_", feature_type, "_", disease, ".tiff"),
     width = 30, 
     height = ifelse(ncol(feature_matrix) < 10, 10, ncol(feature_matrix)),
     units = "cm", compression = "lzw", res = 1200)

plot_box1 + plot_box2

dev.off()



print(warnings())