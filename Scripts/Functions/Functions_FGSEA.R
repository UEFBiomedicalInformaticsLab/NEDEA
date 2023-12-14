# Function to execute FGSEA from RWR probabilities

require(BiocParallel)
require(fgsea)
require(org.Hs.eg.db)
source("Scripts/Functions/Functions_RWR.R")


func_run_FGSEA_on_RWR <- function(rwr_data, enrichment_library, disease, drug_target_type, nproc){
  
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
  
  
  # Create an empty matrix to store the enrichment result
  enrichment_result_mat <- matrix(0, 
                                  nrow = length(enrichment_library), 
                                  ncol = ncol(rwr_data), 
                                  dimnames = list(names(enrichment_library), 
                                                  colnames(rwr_data)))
  
  # For each drug combination run FGSEA and extract the significant NES
  for(drugComb in colnames(rwr_data)){
    
    # Extract the targets of the drug combination
    target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drugComb, drug_target_col, drop = TRUE]
    if(drug_target_type == "all"){
      target_set <- apply(target_set, 1, function(x){paste(x, collapse = ",")})
    }
    target_set <- unique(unlist(strsplit(target_set, ",")))
    suppressMessages(target_set <- select(org.Hs.eg.db, 
                                          keys = target_set, 
                                          columns = "ENSEMBL", 
                                          keytype = "SYMBOL"))
    target_set <- target_set$ENSEMBL
    
    
    # Select the genes for FGSEA
    rwr_data_select <- rwr_data[, drugComb]
    rwr_threshold <- func_RWR_threshold(rwr_data_select[!names(rwr_data_select) %in% target_set])
    rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold$ELB], decreasing = TRUE)
    
    
    # Run FGSEA
    skip_iteration <- FALSE  # Initialize flag variable
    
    tryCatch({
      
      enrichment_result <- fgseaMultilevel(pathways = enrichment_library,
                                           stats = rankedGeneList,
                                           minSize = 5, 
                                           maxSize = 1000, 
                                           scoreType = "pos", 
                                           nproc = nproc,
                                           BPPARAM = MulticoreParam(progressbar = FALSE))
      
      enrichment_result$NES[which(enrichment_result$padj > 0.05)] <- 0
      
      enrichment_result_mat[enrichment_result$pathway, drugComb] <- enrichment_result$NES
      
    }, 
    error = function(e) {
      if (grepl("GSEA statistic is not defined when all genes are selected", conditionMessage(e))) {
        cat(paste("\n\t- Skipping iteration for", drugComb, "due to GSEA statistic error."))
        skip_iteration <- TRUE  # Skip to next iteration of the loop
      } else {
        # Handle other errors here if needed, or re-throw the error
        stop(e)
      }
    })
    if (skip_iteration) {
      next  # Skip to the next iteration of the loop
    }
  }
  return(enrichment_result_mat)
}