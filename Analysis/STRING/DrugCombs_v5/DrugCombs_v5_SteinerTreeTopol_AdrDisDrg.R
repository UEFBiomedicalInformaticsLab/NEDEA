rm(list = ls())



# Extract Steiner tree from adr, disease and drug combination related genes and calculate its topological properties (version 5)


# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(igraph)
library(RCy3)
library(SteinerNet)
library(tidyverse)
library(ggfree)
library(openxlsx)
source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_parallelprocesses.R")





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}


if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))




# Read the network on which to calculate the proxoimities
gene_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")




# Get the disease gene list
# Removing lst from sources which have too high number of genes
disease_genes <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
disease_genes <- disease_genes[grep(pattern = "literature|RNA|Enrichr", x = names(disease_genes), invert = TRUE)]
disease_genes <- unique(unlist(disease_genes, use.names = FALSE))
# disease_net <- subgraph(graph = gene_network, v = V(gene_network)[V(gene_network)$name %in% disease_genes])


# Get the list of genes associated to drug withdrawal related ADRs
adr_genes <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
adr_genes <- unique(unlist(adr_genes, use.names = FALSE))
# adr_net <- subgraph(graph = gene_network, v = V(gene_network)[V(gene_network)$name %in% adr_genes])


















# Load RWR results
load(file = paste0("Analysis/STRING/DrugCombs_v5/", disease, "/dNetRWR050_", disease, ".rda"))


effectiveComb_rwr <- drugCombs_rwr_res_final$effectiveCombinations
adverseComb_rwr <- drugCombs_rwr_res_final$adverseCombinations

# Check if all results have valid number of columns
effectiveComb_rwr <- effectiveComb_rwr[lapply(effectiveComb_rwr, ncol) == 7]
adverseComb_rwr <- adverseComb_rwr[lapply(adverseComb_rwr, ncol) == 7]
rm(drugCombs_rwr_res_final)


cat("\n\nNumbe rof drug combinations: ")
cat(paste0("\n - Effective combinations = ", length(effectiveComb_rwr)))
cat(paste0("\n - Adverse combinations = ", length(adverseComb_rwr), "\n\n"))


drugCombs <- c(effectiveComb_rwr, adverseComb_rwr)

if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 


results <- foreach(i=names(drugCombs), .combine = "rbind", .packages = c("igraph", "SteinerNet", "ggfree")) %dopar% {
  
  # Generate network starting the known drug targets of the combination
  # As only a few targets are available for the drug combinations,
  # we use RWR to extend the gene list. 
  drugTarget_genes <- drugCombs[[i]]$union_seed
  names(drugTarget_genes) <- drugCombs[[i]]$node_name
  drugTarget_genes <- drugTarget_genes[order((drugTarget_genes), decreasing = TRUE)]
  drugTarget_genes <- names(drugTarget_genes[drugTarget_genes > quantile(unique(drugTarget_genes), prob = c(.95))])
  # drugTarget_net <- subgraph(graph = gene_network, v = V(gene_network)[V(gene_network)$name %in% drugTarget_genes])
  
  
  selected_genes <- unique(c(adr_genes, disease_genes, drugTarget_genes))
  
  
  # Extract the Steiner tree connecting teh ADR, disease and drug related genes
  steinertree_res <- steinertree(graph = gene_network,
                                 terminals = V(gene_network)$name[V(gene_network)$name %in% selected_genes],
                                 type = "SP")
  
  subNet <- steinertree_res[[2]]
  # subNet <- subgraph(graph = gene_network, v = V(subNet)$name)
  
  
  # Create a color mapping for nodes
  node_colors <- as.data.frame(matrix(0, nrow = vcount(subNet), ncol = 4,
                                      dimnames = list(V(subNet)$name,
                                                      c("adr", "disease", "drug", "Steiner_node"))))
  node_colors[row.names(node_colors) %in% adr_genes, "adr"] <- "red"
  node_colors[row.names(node_colors) %in% disease_genes, "disease"] <- "blue"
  node_colors[row.names(node_colors) %in% drugTarget_genes, "drug"] <- "green"
  node_colors[!(row.names(node_colors) %in% selected_genes), "Steiner_node"] <- "yellow"
  node_colors[node_colors == 0] <- NA
  node_colors$color <- apply(node_colors, 1, function(x)blend.colors(na.omit(x)))
  
  subNet <- set_vertex_attr(graph = subNet, 
                            name = "color", 
                            index = V(subNet), 
                            value = node_colors$color)
  
  
  
  # # Open in cytoscape
  # createNetworkFromIgraph(igraph = subNet, title = i)
  
  # Calculate the Steiner tree properties
  network_parameters <- data.frame(drugComb = i,
                                   number_steiner_trees = length(steinertree_res) - 1 ,
                                   number_terminal_nodes = length(selected_genes),
                                   number_steiner_nodes = length(V(subNet)$name[!V(subNet)$name %in% selected_genes]),
                                   subNet_nodeNumber = vcount(subNet),
                                   subNet_edgeNumber = ecount(subNet),
                                   subNet_diameter = diameter(graph = subNet, directed = FALSE),
                                   subNet_radius = radius(graph = subNet),
                                   subNet_averagePathLength = mean_distance(graph = subNet),
                                   subNet_clusteringCoefficient = transitivity(graph = subNet, type = "global", isolates = "zero"),
                                   subNet_density = edge_density(graph = subNet),
                                   subNet_connectedComponents = count_components(graph = subNet),
                                   subNet_largestComponentSize = max(components(subNet)$csize),
                                   subNet_powerFit = if(vcount(subNet)>1){fit_power_law(x = igraph::degree(subNet)+1)$alpha}else{"NA"},
                                   subNet_degree = centr_degree(graph = subNet)$centralization,
                                   subNet_closeness = centr_clo(graph = subNet)$centralization,
                                   subNet_betweenness = centr_betw(graph = subNet)$centralization,
                                   subNet_eigen = centr_eigen(graph = subNet)$centralization)
  
  network_parameters
  
}

print(results)
print(class(results))
results <- as.data.frame(t(column_to_rownames(results, "drugComb")))
results <- rownames_to_column(results, "features")

write.xlsx(results, file = paste0("Analysis/STRING/DrugCombs_v5/", disease, "/SteinerTreeTopol_AdrDisDrg_", disease, ".xlsx"), overwrite = TRUE)
saveRDS(results, paste0("Analysis/STRING/DrugCombs_v5/", disease, "/SteinerTreeTopol_AdrDisDrg_", disease, ".rds"))


stopCluster(cl)
unregister_dopar()


print(warnings())