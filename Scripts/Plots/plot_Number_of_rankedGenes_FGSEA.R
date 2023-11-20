set.seed(5081)



# Script to check the number of genes selected genes for FGSEA



# Load libraries
library(foreach)
library(doParallel)
library(tidyverse)
source("Scripts/Functions/Functions_parallelprocesses.R")
source("Scripts/Functions/Functions_RWR.R")


cl <- makeCluster(10)
registerDoParallel(cl) 




summary_df <- foreach(drug_target_type=c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR", "all"), 
                      .packages = c("tidyverse")) %:%
  foreach(disease=c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer"), 
          .packages = c("tidyverse")) %dopar% {
            
            tryCatch({
              
              # Read the RWR results
              rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
              if(!file.exists(rwr_result_file_path)){
                stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
              }
              rwr_result <- readRDS(file = rwr_result_file_path)
              
              
              # Extract the thresholds for RWR probabilities
              rwr_threshold <- sapply(apply(rwr_result, 2, func_RWR_threshold), function(x) x$ELB)
              
              
              
              # Extract the number of seected genes
              tmp1 <- data.frame()
              for(drugComb in colnames(rwr_result)){
                rwr_data_select <- rwr_result[, drugComb]
                rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold[drugComb]], decreasing = TRUE)
                
                tmp1 <- rbind(tmp1, 
                              data.frame("Disease" = disease,
                                         "Drug_target_type" = drug_target_type,
                                         "comb_name" = drugComb, 
                                         "ranked_gene_length" = length(names(rankedGeneList))))
                
              }
              
              row.names(tmp1) <- NULL
              tmp1
              
            }, 
            error = function(e) {
              # Code to handle the error
              message(paste0("Error in iteration ", drug_target_type, ", ", disease, ": ", e$message))
            })
          }

stopCluster(cl)
unregister_dopar()



summary_df <- unlist(summary_df, recursive = FALSE)
summary_df <- do.call(rbind, summary_df)



# Plot the data
tiff("OutputFiles/Plots/rankedGeneList_lengths.tiff",
     width = 25, height = 15,
     units = "cm", compression = "lzw", res = 1200)



plot_data <- summary_df 
plot_data$Drug_target_type <- factor(plot_data$Drug_target_type, levels = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR", "all"))
ggplot(plot_data, aes(x = Disease, y = ranked_gene_length)) +
  geom_boxplot(fill = "grey",
               width = 0.5, lwd = 0.1,
               outlier.shape = 3,
               outlier.size = 0.5) +
  facet_wrap(. ~ Drug_target_type) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.text = element_text(size = 8),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 4),
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  labs(x = "Cancers",
       y = "Number of selected genes")



dev.off()


print(warnings())