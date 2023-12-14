# Function to calculate the proximity between the gene sets

func_calculate_proximity <- function(input_network = input_network,
                                     enrichment_library = enrichment_lib,
                                     disease = disease,
                                     drug_target_type = drug_target_type,
                                     nproc = 5){

  require(org.Hs.eg.db, quietly = TRUE)
  source("Scripts/Functions/Functions_Barabasi_metrics.R")
  
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
  
  
  result <- foreach(drugComb=drugCombs_targets$comb_name, 
                    .combine = "rbind", 
                    .packages = c("org.Hs.eg.db", "igraph"), 
                    .verbose = FALSE) %:%
    foreach(lib_name=names(enrichment_library), 
            .combine = "rbind", 
            .packages = c("org.Hs.eg.db", "igraph"), 
            .verbose = FALSE) %dopar% {
      
      source("Scripts/Functions/Functions_Barabasi_metrics.R")
      
      
      # Extract the targets of the drug combination
      if(drug_target_type == "all"){
        target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drugComb, drug_target_col, drop = FALSE]
        target_set <- apply(target_set, 1, function(x){paste(x, collapse = ",")})
      }else{
        target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drugComb, drug_target_col, drop = TRUE]
      }
      target_set <- unique(unlist(strsplit(target_set, ",")))
      suppressMessages(target_set <- select(org.Hs.eg.db, 
                                            keys = target_set, 
                                            columns = "ENSEMBL", 
                                            keytype = "SYMBOL"))
      target_set <- target_set$ENSEMBL
      
      
      # Extract the enrichment library genes
      enrichment_library_genes <- enrichment_library[[lib_name]]
      
      
      # Calculate the proximities
      prox_separation <- Barabasi_proximity_separation(gene_network = input_network, geneSet1 = target_set, geneSet2 = enrichment_library_genes)
      prox_closest <- Barabasi_proximity_closest(gene_network = input_network, geneSet1 = target_set, geneSet2 = enrichment_library_genes)
      prox_shortest <- Barabasi_proximity_shortest(gene_network = input_network, geneSet1 = target_set, geneSet2 = enrichment_library_genes)
      prox_centre <- Barabasi_proximity_centre(gene_network = input_network, geneSet1 = target_set, geneSet2 = enrichment_library_genes)
      prox_kernel <- Barabasi_proximity_kernel(gene_network = input_network, geneSet1 = target_set, geneSet2 = enrichment_library_genes)
      
      
      # Prepare as dataframe
      tmp1 <- data.frame("drugComb" = drugComb,
                         "lib_name" = lib_name,
                         "prox_Separation" = prox_separation,
                         "prox_Closest" = prox_closest,
                         "prox_Shortest" = prox_shortest,
                         "prox_Centre" = prox_centre,
                         "prox_Kernel" = prox_kernel)
      
      tmp1
    }
  

  

  # Restructure the result as matrix
  result_list <- list()
  for(proximity in c("Centre", "Closest", "Kernel", "Separation", "Shortest")){
    tmp1 <- result[, c("drugComb", "lib_name", paste0("prox_", proximity))]
    tmp1 <- pivot_wider(tmp1, names_from = "drugComb", values_from = paste0("prox_", proximity))
    tmp1 <- column_to_rownames(tmp1, "lib_name")
    result_list[[proximity]] <- tmp1
  }
  
  return(result_list)

}


