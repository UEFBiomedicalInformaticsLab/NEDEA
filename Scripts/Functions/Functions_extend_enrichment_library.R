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

