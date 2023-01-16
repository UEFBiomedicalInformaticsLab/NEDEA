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





# Function to perform FGSEA using the ranked genes from RWR --------------------
# 
# The number of genes to be used can be defined using quantile_prob
# 
# 
# Test:
# 
# 
# enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
# rwr_data <- effectiveComb_rwr
# enrichment_library <- enrichment_lib
# minSize = 1
# maxSize = 200
# scoreType = "std"
# verbose = TRUE
# quantile_prob = 0.9


func_fgsea_from_rwr_probCut <- function(enrichment_library, 
                                        rwr_data, 
                                        minSize = 1, 
                                        maxSize = Inf, 
                                        scoreType = "std", 
                                        quantile_prob = 0.75, 
                                        verbose = TRUE){
  
  require(fgsea)
  require(BiocParallel)
  
  if(!(quantile_prob > 0 & quantile_prob <= 1)){
    stop("quantile_prob must be between 0 and 1. Defines the percentage of ranked genes to be used for FGSEA")
  }
  
  count = 0 # For run status
  total_run <- length(rwr_data) 
  enrichment_result <- list()
  
  p <- bpstart(MulticoreParam(60))
  for(i in names(rwr_data)){
    
    # Rank genes for FGSEA. 
    ranked_genes <- rwr_data[[i]]$union_seed
    names(ranked_genes) <- rwr_data[[i]]$node_name
    ranked_genes <- ranked_genes[order(ranked_genes, decreasing = TRUE)]
    threshold <- as.numeric(quantile(unique(ranked_genes), prob = quantile_prob))
    ranked_genes <- ranked_genes[ranked_genes > threshold]
    
    
    if(verbose){
      # Print running status
      count <- count+1
      cat(paste0("\n\n- Running FGSEA for: (", count , "/", total_run, ") ", i))
      cat(paste0("\n-- Quantile score = ", threshold))
      cat(paste0("\n-- Number of genes used for enrichment = ", length(ranked_genes), "\n"))
      print(head(ranked_genes))
      print(tail(ranked_genes))
    }
    
    # Perform FGSEA
    enrichment_result[[i]] <- fgsea(pathways = enrichment_library,
                                    stats = ranked_genes,
                                    minSize = minSize, maxSize = maxSize, scoreType = scoreType, BPPARAM = p)
  }
  bpstop(p)
  return(enrichment_result)
}





# Function to extract results from FGSEA ---------------------------------------
# 
# After getting the complete FGSEA reult using the function func_fgsea_from_rwr_probCut(),
# use the below function to extract the results as a single dataframe for all combinations
# NA values replaced by 0


func_extract_fgsea_result <- function(enrichment_result, result_type = "NES", enrichment_library, rwr_data){
  switch (result_type,
          "NES" = {
            tmp <- matrix(data = 0, nrow = length(names(enrichment_library)), ncol = length(names(enrichment_result)),
                          dimnames = list(names(enrichment_library), names(enrichment_result)))
            for(i in names(enrichment_result)){
              tmp[enrichment_result[[i]]$pathway, i] <- enrichment_result[[i]]$NES
            }
            tmp[is.na(tmp)] <- 0 # Replace NA values by 0. NA due to overestimation.
            tmp <- as.data.frame(tmp)
            # tmp$mean <- rowMeans(tmp, na.rm = TRUE)
            return(tmp)
          },
          "pval" = {
            tmp <- matrix(nrow = length(names(enrichment_library)), ncol = length(names(enrichment_result)),
                          dimnames = list(names(enrichment_library), names(enrichment_result)))
            for(i in names(enrichment_result)){
              tmp[enrichment_result[[i]]$pathway, i] <- enrichment_result[[i]]$pval
            }
            tmp <- as.data.frame(tmp)
            # tmp$mean <- rowMeans(tmp)
            return(tmp)
          },
          "padj" = {
            tmp <- matrix(nrow = length(names(enrichment_library)), ncol = length(names(enrichment_result)),
                          dimnames = list(names(enrichment_library), names(enrichment_result)))
            for(i in names(enrichment_result)){
              tmp[enrichment_result[[i]]$pathway, i] <- enrichment_result[[i]]$padj
            }
            # tmp <- as.data.frame(tmp)
            # tmp$mean <- rowMeans(tmp)
            return(tmp)
          },
          "size" = {
            tmp <- matrix(nrow = length(names(enrichment_library)), ncol = length(names(enrichment_result)),
                          dimnames = list(names(enrichment_library), names(enrichment_result)))
            for(i in names(enrichment_result)){
              tmp[enrichment_result[[i]]$pathway, i] <- enrichment_result[[i]]$size
            }
            tmp <- as.data.frame(tmp)
            # tmp$mean <- rowMeans(tmp)
            return(tmp)
          },
          "ES" = {
            tmp <- matrix(nrow = length(names(enrichment_library)), ncol = length(names(enrichment_result)),
                          dimnames = list(names(enrichment_library), names(enrichment_result)))
            for(i in names(enrichment_result)){
              tmp[enrichment_result[[i]]$pathway, i] <- enrichment_result[[i]]$ES
            }
            tmp <- as.data.frame(tmp)
            # tmp$mean <- rowMeans(tmp)
            return(tmp)
          },
          cat("ERROR: Result type should be either of the following: NES, pval, padj, size, ES")
  )
}