# Function to extend the enrichment library



# Function to extend enrichment library using Steiner Tree
func_extendEnrichmentLib_bySteinerTree <- function(enrichment_library, input_network, ST_type, nproc = 5){
  
  require(foreach)
  require(doParallel)
  source("Scripts/Functions/Functions_parallelprocesses.R")
  
  
  if(!is.igraph(input_network)){ stop("input_network should be igraph network", .call = FALSE) }
  
  if(!ST_type %in% c("EXA", "SP", "KB", "RSP", "SPM", "ASP")){ stop("Wrong ST_type", .call = FALSE) }
  
  
  cl <- makeCluster(nproc)
  registerDoParallel(cl) 
  
  
  
  res <- foreach(lib_name=names(enrichment_library), 
                 .packages = c("igraph", "SteinerNet"),
                 .final =  function(x){setNames(x, names(enrichment_library))},
                 .verbose = FALSE) %dopar% {
                   
                   enrichment_library_genes <- enrichment_library[[lib_name]]
                   
                   tmp1 <- V(input_network)$name[V(input_network)$name %in% enrichment_library_genes]
                   
                   if(length(tmp1) > 0){
                     steiner_tree <- steinertree(type = ST_type, 
                                                 terminals = tmp1, 
                                                 graph = input_network)
                     V(steiner_tree[[2]])$name
                   }else{NA}
                   
                 }
  
  
  stopCluster(cl)
  unregister_dopar()
  
  return(res)
}




# Function to extend enrichment library using RWR
func_extendEnrichmentLib_byRWR <- function(enrichment_library, input_network, nproc = 5){
  
  require(foreach)
  require(doParallel)
  require(dnet)
  source("Scripts/Functions/Functions_parallelprocesses.R")
  source("Scripts/Functions/Functions_RWR.R")
  
  
  if(!is.igraph(input_network)){ stop("input_network should be igraph network", .call = FALSE) }
  
  # Create seed matrix with the enrichment genes
  rwr_seed_matrix  <- do.call(cbind, lapply(enrichment_library,
                                            function(x){
                                              sapply(V(input_network)$name, function(y) ifelse(y %in% x, 1, 0))
                                            }))
  
  # Run RWR using dRWR
  if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
  cl <- makeCluster(nproc)
  registerDoParallel(cl) 
  
  rwr_result <- dRWR(g = input_network, 
                     setSeeds = rwr_seed_matrix,
                     normalise = "row",
                     restart = 0.5, 
                     normalise.affinity.matrix = "none",
                     multicores = nproc)
  
  colnames(rwr_result) <- colnames(rwr_seed_matrix)
  rownames(rwr_result) <- rownames(rwr_seed_matrix)
  
  
  stopCluster(cl)
  unregister_dopar()
  
  res <- list()
  # For each drug combination run FGSEA and extract the significant NES
  for(lib_name in colnames(rwr_result)){
    
    # Extract the seed genes used
    enrichment_library_genes <- enrichment_library[[lib_name]]
    
    # Select the genes for FGSEA
    rwr_result_select <- rwr_result[, lib_name]
    rwr_threshold <- func_RWR_threshold(rwr_result_select[!names(rwr_result_select) %in% enrichment_library_genes])
    rankedGeneList <- sort(rwr_result_select[rwr_result_select > rwr_threshold$ELB], decreasing = TRUE)
    res[[lib_name]] <- names(rankedGeneList)
  }
  
  return(res)
  
}