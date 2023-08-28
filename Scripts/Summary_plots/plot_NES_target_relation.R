set.seed(5081)



# Plot the distribution of NES of combined efficacy and safety library agaist number of targets



# Load libraries
library(tidyverse)
library(ggfortify)
library(ggpubr)


for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Load RWR results
  load(file = paste0("OutputFiles/Model_train/", disease, "/dNetRWR050_", disease, ".rda"))
  
  
  effectiveComb_rwr <- drugCombs_rwr_res_final$effectiveCombinations
  adverseComb_rwr <- drugCombs_rwr_res_final$adverseCombinations
  rwr_res <- c(effectiveComb_rwr, adverseComb_rwr)
  rwr_res <- as.data.frame(t(as.data.frame(lapply(rwr_res, function(x){sum(x$is_union_seed)}))))
  colnames(rwr_res) <- "Targets"
  
  trainData_feature_matrix  <- readRDS(file = paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))
  trainData_feature_matrix  <- trainData_feature_matrix $NES_CombinedDisAdr2Gene
  trainData_feature_matrix  <- column_to_rownames(trainData_feature_matrix , "Term")
  trainData_feature_matrix  <- as.data.frame(t(trainData_feature_matrix ))
  
  enrich_lib_terms <- colnames(trainData_feature_matrix)
  
  trainData_feature_matrix $Targets <- rwr_res$Targets[match(row.names(trainData_feature_matrix ), row.names(rwr_res))]
  trainData_feature_matrix $Class <- substr(row.names(trainData_feature_matrix), 1, 3)
  
  plot_list <- list()
  for(term in enrich_lib_terms){
    plot_data <- trainData_feature_matrix[, c(term, "Targets", "Class")]
    colnames(plot_data) <- c("Term", "Targets", "Class")
    
    plot_list[[term]] <- ggplot() +
      geom_point(data = plot_data, 
                 aes(x = Term, y = Targets, color = Class),
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
      xlab(str_wrap(term, 30)) +
      ylab(str_wrap("Number of targets", 30)) 
  }
  
  
  # Save the plot
  if(!dir.exists("OutputFiles/Plots/NES_target_relation/")){
    dir.create("OutputFiles/Plots/NES_target_relation/", recursive = TRUE)
  }
  
  tiff(paste0("OutputFiles/Plots/NES_target_relation/NES_target_relation_", disease, ".tiff"),
       width = 40, height = 30,
       units = "cm", compression = "lzw", res = 1200)
  
  plot <- ggarrange(plotlist = plot_list, 
                    common.legend = TRUE, 
                    legend = "bottom")
  
  print(plot)
  
  dev.off()
  
}



print(warnings())