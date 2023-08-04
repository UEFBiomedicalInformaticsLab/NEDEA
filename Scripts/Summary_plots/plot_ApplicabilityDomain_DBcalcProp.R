set.seed(5081)



# Script to plot the applicability domain of molecules
# Notes:
# (a) Details on descriptors: https://egonw.github.io/cdkbook/appmoldescs.html


library(tidyverse)
library(Rcpi)
library(ggfortify)
library(ggpubr)
library(foreach)
library(doParallel)


cl <- makeCluster(2)
registerDoParallel(cl) 


# disease <- "BreastCancer"
plot_list <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  print(disease)
  # Read calculated properties for the drugs
  DrugBank_calc_prop <- read.csv("Databases/DrugBank/drug_calculated_properties.csv")
  DrugBank_calc_prop <- DrugBank_calc_prop[, !colnames(DrugBank_calc_prop) %in% c("source")]
  DrugBank_calc_prop <- reshape(data = DrugBank_calc_prop, direction = "wide", idvar = c("parent_key"), v.names = "value", timevar = c("kind"))
  DrugBank_calc_prop <- DrugBank_calc_prop[, c("parent_key", "value.SMILES")]
  
  
  # Read the training data combinations
  trainData_drugCombs <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_drugCombs <- colnames(trainData_drugCombs$NES_CombinedDisAdr2Gene)
  trainData_drugCombs <- trainData_drugCombs[-1]
  
  # Read the external data combinations
  extData_drugCombs <- readRDS(paste0("OutputFiles/External_predictions/CDCDB/features/CDCDB_OrangeBook_rxotc__", disease, "_CombinedDisAdr2Gene.rds"))
  extData_drugCombs <- colnames(extData_drugCombs)
  
  
  # Extract the properties for the drug combinations
  # drugCombs_prop <- data.frame()
  drugCombs_prop <- foreach(drugComb=c(trainData_drugCombs, extData_drugCombs),
                            .combine = rbind,
                            .errorhandling = "pass",
                            .verbose = TRUE,
                            .packages = c("ggfortify", "Rcpi")) %dopar% {
   # for(drugComb in c(trainData_drugCombs, extData_drugCombs)){     
       print(drugComb)
                              drug1 <- strsplit(drugComb, "__")[[1]][2]
                              drug2 <- strsplit(drugComb, "__")[[1]][3]
                              
                              drug1_smiles <- DrugBank_calc_prop[DrugBank_calc_prop$parent_key == drug1, "value.SMILES"]
                              drug2_smiles <- DrugBank_calc_prop[DrugBank_calc_prop$parent_key == drug2, "value.SMILES"]
                              
                              write(drug1_smiles, paste0("tmp_dir/", drug1, "_smiles.smi"))
                              write(drug2_smiles, paste0("tmp_dir/", drug2, "_smiles.smi"))
                              
                              drug1_mol <- readMolFromSmi(smifile = paste0("tmp_dir/", drug1, "_smiles.smi"), type = "mol")
                              drug2_mol <- readMolFromSmi(smifile = paste0("tmp_dir/", drug2, "_smiles.smi"), type = "mol")
                              
                              drug1_prop <- extractDrugAIO(drug1_mol)
                              colnames(drug1_prop) <- paste0("[Drug1] ", colnames(drug1_prop))
                              
                              drug2_prop <- extractDrugAIO(drug2_mol)
                              colnames(drug2_prop) <- paste0("[Drug2] ", colnames(drug2_prop))
                              
                              tmp1 <- cbind(drug1_prop, drug2_prop)
                              tmp1$drugComb <- drugComb
                              tmp1 <- tmp1[, colnames(tmp1)[grep("parent_key", colnames(tmp1), invert = TRUE)]]
                              tmp1
                              # drugCombs_prop <- rbind(drugCombs_prop, tmp1)
                            }
  row.names(drugCombs_prop) <- NULL
  drugCombs_prop <- column_to_rownames(drugCombs_prop, "drugComb")
  
  
  # Several properties could not be calculated as £D coordinates not found
  # drugCombs_prop <- drugCombs_prop[, colSums(is.na(drugCombs_prop)) != nrow(drugCombs_prop)]
  drugCombs_prop <- drugCombs_prop[, colSums(is.na(drugCombs_prop)) == 0]
  archive <- drugCombs_prop
  # Remove descriptors with all zero count  
  drugCombs_prop <- drugCombs_prop[, colSums(drugCombs_prop) != 0]

  drugCombs_prop$Class <- factor(substr(row.names(drugCombs_prop), 1, 3))
  
  
  

  # PCA
  pca_data <- drugCombs_prop[, !colnames(drugCombs_prop) %in% "Class"]
  pca_res <- prcomp(x = pca_data, center = FALSE, scale. = TRUE)
  pca_scatter <- autoplot(pca_res, 
                          data = drugCombs_prop, 
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

stopCluster(cl)



if(!dir.exists("OutputFiles/Plots/Applicibility_domain/")){
  dir.create("OutputFiles/Plots/Applicibility_domain/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/Applicibility_domain/AD_molDesc_CombinedDisAdr2Gene__CDCDB_OrangeBook_rxotc.tiff"),
     width = 20, height = 20,
     units = "cm", compression = "lzw", res = 1200)
ggarrange(plotlist = plot_list, 
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
dev.off()