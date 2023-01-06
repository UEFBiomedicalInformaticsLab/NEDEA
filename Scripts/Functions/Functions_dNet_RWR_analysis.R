# Functions related to dNet Randomwalk With Restart



# Load libraries
library(igraph)
library(dnet)






# Function to create seed matrix for RWR from drug target interactions
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




















