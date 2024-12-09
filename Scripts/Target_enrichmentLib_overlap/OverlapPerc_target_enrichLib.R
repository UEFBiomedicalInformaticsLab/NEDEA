set.seed(5081)


# Script to check the overlap between drug targets (known and extended) with the compiled efficacy and safety library


# Load libraries
library(openxlsx)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")


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
              
              suppressMessages(mapping <- AnnotationDbi::select(org.Hs.eg.db, 
                                                 keys = target_set, 
                                                 columns = "ENSEMBL", 
                                                 keytype = "SYMBOL"))
              unique(mapping$ENSEMBL)
            })
            
            
            # Read the enrichment library
            enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
            names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
            enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
            names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
            enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)
            
            # Calculate the overlap
            result <- data.frame()
            for(drugComb in names(drugCombs_targets)){
              target_set <- drugCombs_targets[[drugComb]]
              tmp1 <- lapply(enrichment_lib, function(x){
                round((length(intersect(x, target_set))/length(x)) * 100, 3)
              })
              
              result <- rbind(result, data.frame("comb_name" = drugComb,
                                                 "Number_of_targets" = length(target_set),
                                                 rbind(unlist(tmp1)),
                                                 check.names = FALSE,
                                                 row.names = NULL))
              
            }

            result
            
          }


if(!dir.exists("OutputFiles/Geneset_overlap_check/")){dir.create("OutputFiles/Geneset_overlap_check/", recursive = TRUE)}

for(drug_target_type in names(result_final)){
  write.xlsx(x = result_final[[drug_target_type]], 
             file = paste0("OutputFiles/Geneset_overlap_check/OverlapPerc_", drug_target_type, "_target_enrichLib.xlsx"), 
             overwrite = TRUE)
}
saveRDS(result_final, "OutputFiles/Geneset_overlap_check/OverlapPerc_target_enrichLib.rds")


stopCluster(cl)
unregister_dopar()


print(warnings())