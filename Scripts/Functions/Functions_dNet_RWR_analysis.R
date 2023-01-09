# Functions related to dNet Randomwalk With Restart



# Load libraries
library(igraph)
library(dnet)






# Function to create seed matrix for RWR from drug target interactions ---------
# 
# 
# Test:
# disease <- "LungCancer"
# drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugCombs_v3/DrugComb_", disease, "_v3.rds"))
# drug1 = drugCombs$effectiveCombinations[23, "Drug1_DrugBank_drug_id"]
# drug2 = drugCombs$effectiveCombinations[23, "Drug2_DrugBank_drug_id"]
# rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")
# Drug_Target_Net <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")
# func_RWR_seed_from_DTI(rwr_input_network = rwr_input_network, drugs = c(drug1, drug2), drug_target_ixn = Drug_Target_Net)



func_RWR_seed_from_DTI <- function(rwr_input_network, drugs, drug_target_ixn){
  
  if(!class(rwr_input_network) == "igraph"){stop("rwr_input_network should be an igraph object")}
  if(!class(drugs) == "character"){stop("drugs should be a character vector containing DrugBank drug IDs")}
  if(!class(drug_target_ixn) == "data.frame"){stop("drug_target_ixn should be a dataframe")}
  if(!all(c("Node1_drugbank_drug_id", "Node2_ensembl_gene_id") %in% colnames(drug_target_ixn))){
    stop("The dataframe drug_target_ixn should contain the columns 
       Node1_drugbank_drug_id and Node2_ensembl_gene_id")}
  
  
  # Extract targets for the selected drugs from the DTI 
  drug_target_ixn <- drug_target_ixn[drug_target_ixn$Node1_drugbank_drug_id %in% drugs, ]
  drug_target_ixn <- as.data.frame(sapply(drug_target_ixn, as.vector))

  # Filter to keep targets in the network
  drug_target_ixn <- drug_target_ixn[drug_target_ixn$Node2_ensembl_gene_id %in% V(rwr_input_network)$name, ]
  
  
  # Create seed matrix
  rwr_seed_matrix <- matrix(data = 0,
                            nrow = vcount(rwr_input_network),
                            ncol = length(unique(drugs)),
                            dimnames = list(V(rwr_input_network)$name, unique(drugs)))
  
  for(drug in unique(drug_target_ixn$Node1_drugbank_drug_id)){
    tmp1 <- drug_target_ixn[drug_target_ixn$Node1_drugbank_drug_id == drug, ]
    rwr_seed_matrix[tmp1$Node2_ensembl_gene_id, drug] <- 1
  }
  
  return(rwr_seed_matrix)
  
}





# Function for running dNet RWR on drug pairs ----------------------------------
# 
# 
# Test: 
# drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugCombs_v3/DrugComb_", disease, "_v3.rds"))
# drug1 = drugCombs$effectiveCombinations[23, "Drug1_DrugBank_drug_id"]
# drug2 = drugCombs$effectiveCombinations[23, "Drug2_DrugBank_drug_id"]
# rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")
# drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")
# rwr_restart = 0.5
# rwr_norm_input = "none"
# rwr_norm_output = "none"
# verbose = TRUE



func_dNetRWR_on_drugCombination <- function(rwr_input_network, drug1, drug2,
                                            drug_target_ixn,
                                            rwr_restart = 0.5,
                                            rwr_norm_input = "none",
                                            rwr_norm_output = "none",
                                            verbose = TRUE){
  
  print(paste("RWR input network size (nodes, edges):", vcount(rwr_input_network), ecount(rwr_input_network)))
  
  # Extract the seed matrix 
  rwr_input_seed <- func_RWR_seed_from_DTI(rwr_input_network = rwr_input_network, 
                                           drugs = c(drug1, drug2), 
                                           drug_target_ixn = drug_target_ixn)
  
  union_seed <- ifelse(rowSums(rwr_input_seed) > 0, 1, 0)
  intersect_seed <- ifelse(rowSums(rwr_input_seed) == 2, 1, 0)
  rwr_input_seed <- cbind(rwr_input_seed, union_seed, intersect_seed)

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
  rwr_result$is_union_seed <- as.logical(as.data.frame(rwr_input_seed)$union_seed[match(rwr_result$node_name, row.names(rwr_input_seed))])
  
  rwr_result$addEffect <- rowSums(rwr_result[,c(drug1, drug2)])
  rwr_result <- rwr_result[order(rwr_result$union_seed, decreasing = TRUE), ]
  rwr_result <- rwr_result[, c("node_name", "is_union_seed",drug1, drug2, "addEffect", "union_seed", "intersect_seed")]
  row.names(rwr_result) <- NULL
  
  
  if(verbose){
    cat("Quantile-affinities 95%) \n")
    cat(" - with Drug1(", drug1, ") :", quantile(rwr_result[,colnames(rwr_result) == drug1], 0.95), "\n")
    cat(" - with Drug2(", drug2, ") :", quantile(rwr_result[,colnames(rwr_result) == drug2], 0.95), "\n")
    cat(" - with Drug1 + Drug2:", quantile(rwr_result$addEffect, 0.95), "\n")
    cat(" - with Drug1 or Drug2:", quantile(rwr_result$union_seed, 0.95), "\n")
    cat(" - with Drug1 & Drug2:", quantile(rwr_result$intersect_seed, 0.95), "\n")
  }
  
  return(rwr_result)
}


