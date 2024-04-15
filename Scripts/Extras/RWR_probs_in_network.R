
library(igraph)
library(org.Hs.eg.db)
library(RCy3)
source("Scripts/Functions/Functions_FGSEA.R")



# Define global options for this script 
disease <- "BreastCancer"
drug_target_type <- "known"
drugComb <- "DB00188_DB00398"

cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))




# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))





# Read the extended drug target data
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_targets/Drug_targets_extended_", disease, ".rds"))


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










# Extract the known targets of the drug combination

if(drug_target_type != "known"){
  drugCombs_targets_known <- drugCombs_targets[(drugCombs_targets$Drug1_DrugBank_id == strsplit(drugComb, "_")[[1]][1] & 
                                                  drugCombs_targets$Drug2_DrugBank_id == strsplit(drugComb, "_")[[1]][2]), "drugTarget_geneSymbol"]
  
  drugCombs_targets_known <- apply(drugCombs_targets_known, 1, function(x){paste(x, collapse = ",")})
  
  
  drugCombs_targets_known <- unique(unlist(strsplit(drugCombs_targets_known, ",")))
  
  drugCombs_targets_known <- select(org.Hs.eg.db, 
                                    keys = drugCombs_targets_known, 
                                    columns = "ENSEMBL", 
                                    keytype = "SYMBOL")
  drugCombs_targets_known <- drugCombs_targets_known[drugCombs_targets_known$ENSEMBL %in% V(rwr_input_network)$name, "ENSEMBL"]
}



# Extract the targets of the drug combination
drugCombs_targets <- drugCombs_targets[(drugCombs_targets$Drug1_DrugBank_id == strsplit(drugComb, "_")[[1]][1] & 
                                          drugCombs_targets$Drug2_DrugBank_id == strsplit(drugComb, "_")[[1]][2]), drug_target_col]

drugCombs_targets <- apply(drugCombs_targets, 1, function(x){paste(x, collapse = ",")})


drugCombs_targets <- unique(unlist(strsplit(drugCombs_targets, ",")))

drugCombs_targets <- select(org.Hs.eg.db, 
                            keys = drugCombs_targets, 
                            columns = "ENSEMBL", 
                            keytype = "SYMBOL")


drugCombs_targets <- drugCombs_targets[drugCombs_targets$ENSEMBL %in% V(rwr_input_network)$name, "ENSEMBL"]













# Read the RWR results
rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
if(!file.exists(rwr_result_file_path)){
  stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
}
rwr_result <- readRDS(file = rwr_result_file_path)


# Extract the gene list influenced by the drug target
rwr_data <- rwr_result[, drugComb, drop = FALSE]
rwr_threshold <- sapply(apply(rwr_data, 2, func_RWR_threshold), function(x) x$ELB)
rankedGeneList <- names(sort(rwr_data[which(rwr_data[,drugComb] > rwr_threshold[drugComb]), drugComb], 
                             decreasing = TRUE))


# Extract the probabilities for the targets and the genes selected by RWR
rwr_data <- as.data.frame(as.matrix(rwr_data[row.names(rwr_data) %in% c(drugCombs_targets, rankedGeneList), , drop = FALSE]))

rwr_data$target <- "not_target"

if(drug_target_type != "known"){
  rwr_data[(row.names(rwr_data) %in% drugCombs_targets_known), "target"] <- "known_target"
  rwr_data[(row.names(rwr_data) %in% drugCombs_targets) & !(row.names(rwr_data) %in% drugCombs_targets_known), "target"] <- "extended_target"
}

if(drug_target_type == "known"){
  rwr_data[(row.names(rwr_data) %in% drugCombs_targets), "target"] <- "known_target"
}

rwr_data[[drugComb]] <- as.numeric(rwr_data[[drugComb]])


# Extarct subnetwork
subnet <- induced.subgraph(graph = rwr_input_network, vids = c(drugCombs_targets, rankedGeneList))
subnet <- set_vertex_attr(graph = subnet, 
                          name = "RWR_prob", 
                          index = V(subnet)[V(subnet)$name %in% row.names(rwr_data)], 
                          value = rwr_data[[drugComb]])
subnet <- set_vertex_attr(graph = subnet, 
                          name = "target", 
                          index = V(subnet)[V(subnet)$name %in% row.names(rwr_data)], 
                          value = rwr_data$target)

vertex_attr_names(subnet)


createNetworkFromIgraph(subnet, title = "Plot Network", collection = "Network")




library(visNetwork)


visIgraph(subnet, ) %>%
  visConfigure(enabled = TRUE)