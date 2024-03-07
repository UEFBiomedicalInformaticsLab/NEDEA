

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  # Read the drug target data
  drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
  drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, 
                                       drugCombs_targets$Drug2_DrugBank_id, 
                                       sep = "_")
  drugCombs_targets <- drugCombs_targets[,c("comb_name", "drugTarget_ensembl_id"), drop = FALSE]
  drugCombs_targets <- split(drugCombs_targets$drugTarget_ensembl_id, drugCombs_targets$comb_name)
  drugCombs_targets <- lapply(drugCombs_targets, function(x){unlist(strsplit(x, ","))})
  
  
  
  
  # Read the enrichment libraries
  enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
  enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
  names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
  enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)
  
  
  
  # Create dataframe to store results
  count_df <- as.data.frame(matrix(0, 
                                   nrow = length(drugCombs_targets), 
                                   ncol = length(enrichment_lib), 
                                   dimnames = list(names(drugCombs_targets), 
                                                   names(enrichment_lib))
  ))
  
  
  # Check the number of overlap between drug targets and the genesets
  for(drugComb in row.names(count_df)){
    for(lib_name in colnames(count_df)){
      count_df[drugComb, lib_name] <- length(intersect(drugCombs_targets[[drugComb]], 
                                                       enrichment_lib[[lib_name]]))
    }
  }
  
  # Calculate for how many gene sets there is atleast one overlapping target
  count_df$all_overlap_count <- apply(count_df, 1, function(x)sum(x > 0))
  count_df$disease_overlap_count <- apply(count_df[, grep("^\\[DISEASE\\]", colnames(count_df))], 1, function(x)sum(x > 0))
  count_df$adr_overlap_count  <- apply(count_df[, grep("^\\[ADR\\]", colnames(count_df))], 1, function(x)sum(x > 0))
  
  
  # Add drug target count
  count_df$drug_target_count <- 0
  for(drugComb in row.names(count_df)){
    count_df[drugComb, ]$drug_target_count <- length(drugCombs_targets[[drugComb]])
  }
  
  
  
  # Add category info
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, 
                                   drugCombs_cat$Drug2_DrugBank_id, 
                                   sep = "_")
  count_df$class <- drugCombs_cat$class_EffAdv[match(row.names(count_df), drugCombs_cat$comb_name)]
  
  
  
  
  tmp1 <- count_df[count_df$all_overlap_count == 0, ]
  tmp2 <- count_df[count_df$all_overlap_count == 0 & !is.na(count_df$class), ]
  
  
  # print result
  print("---------")
  print(paste0("Disease: ", disease))
  print(paste0("Total number of drug combs: ", nrow(count_df)))
  print(paste0("Without any overlap: ", nrow(tmp1)))
  print(paste0("Without any overlap but with class: ", nrow(tmp2)))
  
  print( (nrow(tmp1) / nrow(count_df)) * 100 )
  
  # print(count_df["DB00290_DB04845", c("all_overlap_count", "disease_overlap_count", "adr_overlap_count", "drug_target_count", "class")])
}

