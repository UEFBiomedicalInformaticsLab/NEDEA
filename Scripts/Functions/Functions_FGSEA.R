# Function to execute FGSEA from RWR probabilities

require(BiocParallel)
require(fgsea)
source("Scripts/Functions/Functions_RWR.R")


func_run_FGSEA_on_RWR <- function(rwr_data, enrichment_library){
  
  # Extract the thresholds for RWR probabilities
  rwr_threshold <- sapply(apply(rwr_data, 2, func_RWR_threshold), function(x) x$ELB)
  
  # Create an empty matrix to store the enrichment result
  enrichment_result_mat <- matrix(0, 
                                  nrow = length(enrichment_library), 
                                  ncol = ncol(rwr_data), 
                                  dimnames = list(names(enrichment_library), 
                                                  colnames(rwr_data)))
  
  # For each drug combination run FGSEA and extract the significant NES
  for(drugComb in colnames(rwr_data)){
    
    # rankedGeneList <- sort(rwr_data[which(rwr_data[,drugComb] > rwr_threshold[drugComb]), drugComb], 
    #                        decreasing = TRUE)
    
    rwr_data_select <- rwr_data[, drugComb]
    rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold[drugComb]], decreasing = TRUE)
    

    skip_iteration <- FALSE  # Initialize flag variable
    
    tryCatch({
      
      enrichment_result <- fgseaMultilevel(pathways = enrichment_library,
                                           stats = rankedGeneList,
                                           minSize = 5, 
                                           maxSize = 500, 
                                           scoreType = "pos", 
                                           BPPARAM = SerialParam())
      
      # enrichment_result$NES[which(enrichment_result$pval > 0.01)] <- 0
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