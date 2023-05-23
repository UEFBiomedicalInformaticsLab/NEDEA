# Validation using the CDCDB drug combinations (RX/OTC type)

# Notes:
# (1) The two drug combinations from FDA Orange Book in C-DCDC  
# (2) Remove drug combinations that contain at least one drug already in 
#     the training set across all cancers.
# (3) Removed drug combinations that do not contain any reported drug targets
# (4) Run RWR-FGSEA and then use the selected model to predict





# Load libraries
library(foreach)
library(doParallel)
library(caret)
library(openxlsx)
library(tidyverse)
source("Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("Scripts/Functions/Functions_parallelprocesses.R")





# Read drug combinations used in training set
drugCombs_training <- list()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", disease, ".rds"))
  drugCombs <- do.call(rbind, drugCombs)
  row.names(drugCombs) <- NULL
  drugCombs_training[[disease]] <- drugCombs
}
drugCombs_training <- do.call(rbind, drugCombs_training)
row.names(drugCombs_training) <- NULL
drugCombs_training <- unique(drugCombs_training)
drugs_training <- unique(c(drugCombs_training$Drug1_DrugBank_drug_id, drugCombs_training$Drug2_DrugBank_drug_id))



# Read the Orange Book drug combinations from CDCDB
drugCombs <- readRDS("InputFiles/ReferenceList/CDCDB_drugCombinations_OrangeBook_rxotc.rds")



# Remove combinations that contain atleast one drug participating in training
remove_index <- c()
for(i in 1:nrow(drugCombs)){
  drug1 <- drugCombs[i, "Drug1_DrugBank_drug_id"]
  drug2 <- drugCombs[i, "Drug2_DrugBank_drug_id"]
  tmp <- drugCombs_training[drugCombs_training$Drug1_DrugBank_drug_id == drug1 & drugCombs_training$Drug2_DrugBank_drug_id == drug2, ]
  if(nrow(tmp) > 0){remove_index <- c(remove_index, i)} 
}
if(length(remove_index) > 0){drugCombs <- drugCombs[-remove_index,]}
rm(tmp)



# Read drug target interactions from drug bank
drug_target_ixn  <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")



# Check if the drugs forming the pairs have reported targets
drugCombs <- drugCombs[(drugCombs$Drug1_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id 
                                                & drugCombs$Drug2_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id), ]
drugCombs <- unique(drugCombs)
row.names(drugCombs) <- NULL



# Read the network on which to run RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")



# Execute RWR
nproc <- detectCores()/2 #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 

drugCombs_rwr <- foreach (i=1:nrow(drugCombs),
                           .packages = c("igraph", "tidyr", "dnet")) %dopar% {
                             tmp <- func_dNetRWR_on_drugCombination(rwr_input_network = rwr_input_network,
                                                                    drug1 = drugCombs[i,c("Drug1_DrugBank_drug_id")], 
                                                                    drug2 = drugCombs[i,c("Drug2_DrugBank_drug_id")],
                                                                    drug_target_ixn = drug_target_ixn,
                                                                    rwr_norm_input = "laplacian",
                                                                    rwr_norm_output = "none", 
                                                                    rwr_restart = 0.5)
                             tmp
                           }


for(i in 1:length(drugCombs_rwr)){
  names(drugCombs_rwr)[i] <- paste("Unk", drugCombs[i,c("Drug1_DrugBank_drug_id")], 
                                          drugCombs[i,c("Drug2_DrugBank_drug_id")], sep = "__")
}

stopCluster(cl)
unregister_dopar()








# Create predictions using the six cancer models

validation_predictions <- list()
validation_predictions[["Summary"]] <- data.frame()
for(disease in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  
  cat(paste0("\nExecuting for: ", disease, "\n"))
  
  # Run FGSEA on combined Disease2Gene and drug withdrawal based Adr2Gene library (efficacy + safety)
  enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
  names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
  enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
  names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))

  enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

  enrichment_res <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                          rwr_data = drugCombs_rwr, 
                                                          quantile_prob = 0.9,
                                                          verbose = FALSE)

  drugComb_NES <- func_extract_fgsea_result(enrichment_result = enrichment_res,
                                                 result_type = "NES",
                                                 enrichment_library = enrichment_lib)
  
 
  # Read the selected model
  model <- readRDS(paste0("OutputFiles/Final_model/FinalModel_", disease, "_rf_none_CombinedDisAdr2Gene.rds"))          

  # Prepare the input data
  data <- as.data.frame(t(drugComb_NES))
  # preProcess <- preProcess(data, method = c("zv", "center", "scale"))
  # data <- predict(object = preProcess, newdata = data)

  # Predict
  predictions <- predict(object = model, newdata = data, type = "prob") 
  predictions$Class <- predict(object = model, newdata = data, type = "raw")
  predictions <- rownames_to_column(predictions, "drugCombs")
  predictions <- predictions[order(predictions$Eff, decreasing = TRUE),]
  validation_predictions[[disease]] <- predictions
  validation_predictions[["Summary"]] <- rbind(validation_predictions[["Summary"]], 
                                                data.frame(Disease = disease,
                                                            predicted_Eff = nrow(predictions[predictions$Class == "Eff", ]),
                                                            predicted_Adv = nrow(predictions[predictions$Class == "Adv", ])))
}


# Save as file
if(!dir.exists("OutputFiles/Validation/CDCDB/")){
  dir.create("OutputFiles/Validation/CDCDB/", recursive = TRUE)
}                          
write.xlsx(validation_predictions, "OutputFiles/Validation/CDCDB/validation_CDCDB_OrangeBook_rxotc.xlsx", rowNames = FALSE, overwrite = TRUE)






# Create box plots for efficacy scores
validation_predictions <- validation_predictions[names(validation_predictions) != c("Summary")]
validation_predictions <- do.call(rbind, validation_predictions)
validation_predictions <- rownames_to_column(validation_predictions, "Disease")
validation_predictions$Disease <- gsub("\\.[0-9]+$", "", validation_predictions$Disease)




tiff("OutputFiles/Validation/CDCDB/validation_CDCDB_OrangeBook_rxotc.tiff", 
     width = 6, height = 5, 
     units = "cm", compression = "lzw", res = 1200)

ggplot(validation_predictions, aes(x = Disease, y = Eff, fill = Class)) + #
  geom_violin(lwd = 0.1) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.1, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.1,),
        strip.text = element_text(margin = margin(1,1,1,1)),
        text = element_text(linewidth = 2.5), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(linewidth = 1.5),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1)) +
  # scale_fill_manual(values = c("#008080", "#ffa500", "#00ff7f", "#00bfff", "#deb887")) + 
  ggtitle(paste0("Predictions for C-DCDB drug combinations (Orange Book, RX/OTC)")) +
  xlab("Cancers") + ylab("Efficacy scores") + #
  labs(colour = "Class : ")

dev.off()



print(warnings())