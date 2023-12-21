set.seed(5081)


# Script to check the overlap score between drug targets (known and extended) with the extended efficacy and safety library


# Load libraries
library(openxlsx)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")


scale_values <- function(x){ (x - min(x)) / (max(x) - min(x)) }

cl <- makeCluster(20)
registerDoParallel(cl)


result_final <- foreach(drug_target_type = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR", "all"),
                        .final =  function(x){setNames(x, c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR", "all"))},
                        .packages = c("tidyverse", "org.Hs.eg.db")) %:%
  foreach(disease = c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer"), 
          .final =  function(x){setNames(x, c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer"))},
          .packages = c("tidyverse", "org.Hs.eg.db")) %dopar% {
            
            
            # Select the column containing the drug target based on th user input
            switch(drug_target_type,
                   "known" = { drug_target_col <- c("drugTarget_geneSymbol") },
                   "PS" = { drug_target_col <- c("ext_PS_targets") },
                   "SIGNOR" = { drug_target_col <- c("ext_SIGNOR_targets") },
                   "NPA" = { drug_target_col <- c("ext_NPA_targets") },
                   "RI" = { drug_target_col <- c("ext_RI_targets") },
                   "KEGG" = { drug_target_col <- c("ext_KEGG_targets") },
                   "all" = { drug_target_col <- c("drugTarget_geneSymbol", "ext_PS_targets", 
                                                  "ext_SIGNOR_targets", "ext_NPA_targets", 
                                                  "ext_RI_targets", "ext_KEGG_targets") })
            
            # Read the drug targets
            drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
            drugCombs_targets <- as.data.frame(drugCombs_targets)
            row.names(drugCombs_targets) <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
            drugCombs_targets <- drugCombs_targets[, drug_target_col, drop = FALSE]
            drugCombs_targets <- as.data.frame(apply(drugCombs_targets, 1, function(x){paste(x, collapse = ",")}))
            colnames(drugCombs_targets) <- "Targets"
            drugCombs_targets <- apply(drugCombs_targets, 1, function(x){
              target_set <- unique(unlist(strsplit(x, ",")))
              
              suppressMessages(mapping <- select(org.Hs.eg.db, 
                                                 keys = target_set, 
                                                 columns = "ENSEMBL", 
                                                 keytype = "SYMBOL"))
              unique(mapping$ENSEMBL)
            })
            
            
            # Read the enrichment library
            enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries_extended/Disease2Gene_", disease, "_extendedLib.rds"))
            names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
            enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries_extended/curatedAdr2Gene_extendedLib.rds")
            names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
            enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)
            
            # Get scaled values for scoring 
            target_scaling_value <- scale_values(lengths(drugCombs_targets))
            enrichLib_scaling_value <- scale_values(lengths(enrichment_lib))
            
            
            # Calculate the overlap
            result <- data.frame()
            for(drugComb in names(drugCombs_targets)){
              
              target_set <- drugCombs_targets[[drugComb]]
              for(lib_name in names(enrichment_lib)){
                
                tmp1 <- (length(intersect(enrichment_lib[[lib_name]], target_set)) / length(enrichment_lib[[lib_name]])) + target_scaling_value[[drugComb]] + enrichLib_scaling_value[[lib_name]]
                
                result <- rbind(result, data.frame("comb_name" = drugComb,
                                                   "Enrichment_lib_name" = lib_name,
                                                   "Number_of_targets" = length(target_set),
                                                   "Score" = tmp1,
                                                   check.names = FALSE,
                                                   row.names = NULL))
                
              }
            }
            
            result <- pivot_wider(data = result, names_from = Enrichment_lib_name, values_from = Score)
            result
            
          }


if(!dir.exists("OutputFiles/Geneset_overlap_check/")){dir.create("OutputFiles/Geneset_overlap_check/", recursive = TRUE)}

for(drug_target_type in names(result_final)){
  write.xlsx(x = result_final[[drug_target_type]], 
             file = paste0("OutputFiles/Geneset_overlap_check/OverlapScore_", drug_target_type, "_target_extendedEnrichLib.xlsx"), 
             overwrite = TRUE)
}
saveRDS(result_final, "OutputFiles/Geneset_overlap_check/OverlapScore_target_extendedEnrichLib.rds")


stopCluster(cl)
unregister_dopar()


print(warnings())