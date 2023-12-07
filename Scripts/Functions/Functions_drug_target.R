# Functions used in extracting and extending drug targets



# Function to convert Ensembl IDs to gene symbol
convert_ensembl_to_symbol <- function(ids_str) {
  ids <- unique(unlist(strsplit(ids_str, ",")))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys = ids, 
                                     columns = "SYMBOL", 
                                     keytype = "ENSEMBL"))
  return(paste(mapping$SYMBOL, collapse = ","))
}



# Function to convert gene symbol to Ensembl IDs
convert_symbol_to_ensembl <- function(ids_str) {
  ids <- unique(unlist(strsplit(ids_str, ",")))
  suppressMessages(mapping <- select(org.Hs.eg.db, keys = ids, 
                                     columns = "ENSEMBL", 
                                     keytype = "SYMBOL"))
  return(paste(mapping_df$ensembl_gene_id, collapse = ","))
}



# Function to extend known drug targets with any 'interaction' database
extend_targets_ixns_database <- function(ids_str, ixns_data) {
  known_drug_targets <- unique(unlist(strsplit(ids_str, ",")))
  # Extend drug targets using the provided interaction database
  extended_targets <- ixns_data[ixns_data$source_genesymbol %in% known_drug_targets,]
  extended_targets <- extended_targets[which(extended_targets$curation_effort > 1 & extended_targets$is_directed),]
  all_targets <- unique(c(extended_targets$source_genesymbol, extended_targets$target_genesymbol, known_drug_targets))
  return(list("ext_targets" = paste(all_targets, collapse = ","), 
              "ext_tar_cnt" = length(all_targets)))
}



# Function to extend drug targets by using a kegg pathway
extend_targets_kegg <- function(initial_targets, kegg_data) {
  
  extended_targets <- unique(unlist(strsplit(initial_targets, ",")))
  
  # Continue until no new targets are found
  new_targets_found <- TRUE
  while (new_targets_found) {
    new_targets_found <- FALSE
    
    # Find interactions where the source gene is in our list of targets
    new_targets <- kegg_data %>%
      filter(genesymbol_source %in% extended_targets) %>%
      pull(genesymbol_target)
    
    # Add unique new targets to our list
    original_length <- length(extended_targets)
    extended_targets <- unique(c(extended_targets, new_targets))
    
    if (length(extended_targets) > original_length) {
      new_targets_found <- TRUE
    }
  }
  return(extended_targets)
}



extend_targets_kegg_list <- function(initial_targets, list_kegg) {
  # Loop to extend targets for each KEGG data frame
  results <- list()
  for (i in seq_along(list_kegg)) {
    kegg_data <- list_kegg[[i]]
    extended_targets <- extend_targets_kegg(initial_targets, kegg_data)
    results[[names(list_kegg)[i]]] <- extended_targets
  }
  return(list("ext_kegg_targets" = paste(Reduce('union',results), collapse = ","), 
              "ext_kegg_tar_cnt" = length(Reduce('union',results))))
  return()
}




# Download the protein-protein interactions based on KEGG

# download_KEGG_ixns <- function(drug_targets, pathway_ids) {
#   # overlap <- rep(0,length(pathway_ids))
#   list_kegg_pw <- list()
#   nn <- c()
#   for(i in 1:length(pathway_ids)) {
#     print(pathway_ids[i])
#     hsap <- kegg_pathway_download(pathway_ids[i], process = FALSE)
#     kegg_nets <- kegg_process(hsap$entries, hsap$relations, simplify = FALSE)
#     list_kegg_pw[[length(list_kegg_pw)+1]] <- kegg_nets
#     print(length(intersect(drug_targets, union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target))))
#     print(length(union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target)))
#     if(length(union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target)) > 0)
#       nn <- c(nn, pathway_ids[i])
#   }
#   names(list_kegg_pw) <- nn
#   return(list_kegg_pw)
# }

download_KEGG_ixns <- function(pathway_ids){
  list_kegg_pw <- list()
  # summary_df <- data.frame()
  
  for(path in pathway_ids){
    # print(path)
    kegg_net <- kegg_pathway_download(pathway_id = path, process = FALSE)
    kegg_net <- kegg_process(kegg_net$entries, kegg_net$relations, simplify = FALSE)
    
    # summary_df <- rbind(summary_df, 
    #                     data.frame("Pathway" = path,
    #                                "Num_genes_in_pathway" = length(union(kegg_net$genesymbol_source, kegg_net$genesymbol_target)),
    #                                "Num_drug_targets_in_pathway" = length(intersect(drug_targets, union(kegg_net$genesymbol_source, kegg_net$genesymbol_target)))
    #                     ))
    if(length(union(kegg_net$genesymbol_source, kegg_net$genesymbol_target)) > 0){
      list_kegg_pw[[path]] <- kegg_net
    }
  }
  
  # cat("\n\nSummary of the downloaded KEGG interactions:\n")
  # print(summary_df)
  return(list_kegg_pw)
}