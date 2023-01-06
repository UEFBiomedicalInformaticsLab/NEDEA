source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_dNet_RWR_analysis.R")

disease <- "LungCancer"
drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugCombs_v3/DrugComb_", disease, "_v3.rds"))
drug1 = drugCombs$effectiveCombinations[23, "Drug1_DrugBank_drug_id"]
drug2 = drugCombs$effectiveCombinations[23, "Drug2_DrugBank_drug_id"]
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")

drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")


rwr_restart = 0.5
rwr_norm_input = "none"
rwr_norm_output = "none"
verbose = TRUE


# func_dNetRWR_on_drugCombination <- function(rwr_input_network, drug1, drug2, 
#                                             drug_target_ixn,
#                                             rwr_restart = 0.5, 
#                                             rwr_norm_input = "none", 
#                                             rwr_norm_output = "none", 
#                                             verbose = TRUE){}

rwr_input_seed <- func_RWR_seed_from_DTI(rwr_input_network = rwr_input_network, 
                                         drugs = c(drug1, drug2), 
                                         drug_target_ixn = drug_target_ixn)
  
  
  
# Execute dNet RWR
rwr_result <- suppressMessages(dRWR(g = rwr_input_network, 
                                    setSeeds = rwr_input_seed, 
                                    normalise = rwr_norm_input, 
                                    restart = rwr_restart, 
                                    normalise.affinity.matrix = rwr_norm_output))


# Prepare result for export
rwr_result <- as.data.frame(as.matrix(rwr_result))
colnames(rwr_result) <- colnames(rwr_input_seed)
rwr_result$node_name <- V(rwr_input_network)$name
rwr_result$seed <- as.logical(as.data.frame(rwr_input_seed)$union_seed[match(rwr_result$node_name, row.names(rwr_input_seed))])
rwr_result$addEffect <- rowSums(rwr_result[,query_nodes])
# rwr_result <- rwr_result[order(rwr_result$union_seed, decreasing = TRUE), ]
rwr_result <- rwr_result[, c("node_name", "node_type", "seed", query_nodes, "addEffect", "union_seed", "intersect_seed")]
row.names(rwr_result) <- NULL