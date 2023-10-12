set.seed(5081)



# Script to plot PCA

# Load libraries
library(FactoMineR)
library(factoextra )

disease <- "BreastCancer"
drug_target_type <- "all"


# Read the FGSEA result
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
fgsea_result <- fgsea_result$combined_efficacy_safety


# Read the extended drug target data
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_targets/Drug_targets_extended_", disease, ".rds"))
drugCombs_targets$combination_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")




# pca_res <- PCA(t(fgsea_result[-which(apply(fgsea_result,1,var) < 0.01),]), 
#                scale.unit = FALSE, 
#                ncp = 5, 
#                graph = FALSE)
# 
# fviz_pca_ind(pca_res, axes = c(1, 2), col.ind = drugCombs_targets$drug_class, label = "none")
# fviz_pca_ind(pca_res, axes = c(1, 2), col.ind = as.factor(drugCombs_targets$totSyn), label = "none")




pca_data <- t(fgsea_result[-which(apply(fgsea_result,1,var) < 0.01),])

pca <- prcomp(x = pca_data, center = FALSE, scale. = FALSE)
pca_res <- predict(pca, pca_data)

plot_data <- data.frame(pca_res)
plot_data$Class <- drugCombs_targets$drug_class[match(row.names(plot_data), drugCombs_targets$combination_name)]
plot_data$totSyn <- drugCombs_targets$totSyn[match(row.names(plot_data), drugCombs_targets$combination_name)]
plot_data$totSyn <- as.factor(plot_data$totSyn)






pca_scatter_1 <- ggplot() +
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



pca_scatter_2 <- ggplot() +
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, color = totSyn), 
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






print(warnings())