

library(org.Hs.eg.db)
source("Scripts/Functions/Functions_RWR.R")


disease <- "BreastCancer"
drug_target_type <- "known"
drug_comb_name <- "DB08881_DB00642"




# Extract the drug targets 
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, 
                                     drugCombs_targets$Drug2_DrugBank_id, 
                                     sep = "_")
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

target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drug_comb_name, drug_target_col]
target_set <- unlist(strsplit(target_set, ","))
suppressMessages(target_set <- select(org.Hs.eg.db, 
                                      keys = target_set, 
                                      columns = "ENSEMBL", 
                                      keytype = "SYMBOL"))
target_set <- target_set$ENSEMBL




# Extract the genes selected for FGSEA
rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
if(!file.exists(rwr_result_file_path)){
  stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
}
rwr_result <- readRDS(file = rwr_result_file_path)
rwr_threshold <- sapply(apply(rwr_result, 2, func_RWR_threshold), function(x){x$ELB})
rwr_data_select <- rwr_result[, drug_comb_name]





probs_1 <- sort(rwr_data_select, decreasing = TRUE)
rwr_threshold <- func_RWR_threshold(probabilities = probs_1)
rwr_threshold <- rwr_threshold$ELB
rankedGeneList_1 <- sort(probs_1[probs_1 > rwr_threshold], decreasing = TRUE)





probs_2 <- sort(rwr_data_select, decreasing = TRUE)
probs_2 <- probs_2[!names(probs_2) %in% target_set] 
rwr_threshold <- func_RWR_threshold(probabilities = probs_2)
rwr_threshold <- rwr_threshold$ELB
rankedGeneList_2 <- sort(probs_2[probs_2 > rwr_threshold], decreasing = TRUE)



p1 <- plot(probs_1)
p2 <- plot(probs_2)


length(probs_1)
length(probs_2)


length(rankedGeneList_1)
length(rankedGeneList_2)
