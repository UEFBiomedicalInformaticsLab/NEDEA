rm(list = ls())



# Calculate different Barabasi metrics between targets of the two drugs and the ADR related genes



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(tidyverse)
library(openxlsx)
library(igraph)
source("Scripts/Functions/Functions_Barabasi_metrics.R")
source("Scripts/Functions/Functions_parallelprocesses.R")





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
  stop("--nproc should be be an integer.", call.=FALSE)
}
}


# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Read the drug combinations
drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugCombs_v5/DrugComb_", disease, "_v5.rds"))
cat("\n\nNumber of drug combinations:\n")
lapply(drugCombs, nrow)


drugCombs <- do.call(rbind, drugCombs)
drugCombs$Class <- substr(x = row.names(drugCombs), start = 1, stop = 3) 
drugCombs$Class <- gsub(pattern = "^eff", replacement = "Eff", x = drugCombs$Class)
drugCombs$Class <- gsub(pattern = "^adv", replacement = "Adv", x = drugCombs$Class)
row.names(drugCombs) <- NULL



# Read drug target interactions
drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")


# Get the list of genes associated to drug withdrawal related ADRs
adr_genes <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
adr_genes <- unique(unlist(adr_genes, use.names = FALSE))
# adr_net <- subgraph(graph = gene_network, v = V(gene_network)[V(gene_network)$name %in% adr_genes])



# Read the network on which to calculate the proximities
gene_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")
print(paste("Input network size (before clean):", vcount(gene_network), ecount(gene_network)))
isolated <- names(which(degree(gene_network)==0)) # remove isolated nodes
gene_network <- delete.vertices(gene_network, isolated)
cc <- components(gene_network)
largest_components <- V(gene_network)[cc$membership == which.max(cc$csize)] # remove connected components
gene_network <- induced_subgraph(gene_network, largest_components)
print(paste("Input network size (after clean):", vcount(gene_network), ecount(gene_network)))


# Using external cluster instead of clustering through each function call
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 


proximity_matrix  <- foreach(i=1:nrow(drugCombs),
                             .combine = rbind,
                             .packages = c("igraph", "tidyr")) %dopar%  {
                               
                               drug1 <- drugCombs[i, "Drug1_DrugBank_drug_id"]
                               drug2 <- drugCombs[i, "Drug2_DrugBank_drug_id"]
                               
                               cat(paste0("\nCalculating for: ", drugCombs[i,"Class"], "__", drug1, "__", drug2, "\n"))
                               
                               
                               # Retrieve the targets for the two drugs
                               drug1_targets <- drug_target_ixn[drug_target_ixn$Node1_drugbank_drug_id == drug1, "Node2_ensembl_gene_id"]
                               drug2_targets <- drug_target_ixn[drug_target_ixn$Node1_drugbank_drug_id == drug2, "Node2_ensembl_gene_id"]
                               drug_targets <- unique(c(drug1_targets, drug2_targets))
  
                               # Calculate the proximities
                               proximity_closest <- Barabasi_proximity_closest(gene_network, drug_targets, adr_genes)
                               proximity_shortest <- Barabasi_proximity_shortest(gene_network, drug_targets, adr_genes)
                               proximity_centre <- Barabasi_proximity_centre(gene_network, drug_targets, adr_genes)
                               proximity_kernel <- Barabasi_proximity_kernel(gene_network, drug_targets, adr_genes)
                               proximity_separation <- Barabasi_proximity_separation(gene_network, drug_targets, adr_genes)
                               
                               # Prepare for export as dataframe
                               tmp <- data.frame(drugComb = paste(drugCombs[i,"Class"], drug1, drug2, sep = "__"),
                                                 drug1_target_number = length(drug1_targets),
                                                 drug2_target_number = length(drug2_targets),
                                                 adr_gene_number = length(adr_genes),
                                                 proximity_closest = proximity_closest,
                                                 proximity_shortest = proximity_shortest,
                                                 proximity_centre = proximity_centre,
                                                 proximity_kernel = proximity_kernel,
                                                 proximity_separation = proximity_separation)
                               tmp
                             }

stopCluster(cl)
unregister_dopar()



rownames(proximity_matrix) <- proximity_matrix[, "drugComb"]
proximity_matrix <- proximity_matrix[,!colnames(proximity_matrix) %in% "drugComb"]
proximity_matrix <- as.data.frame(t(proximity_matrix))
proximity_matrix <- rownames_to_column(proximity_matrix, "features")


saveRDS(proximity_matrix, file = paste0("Analysis/STRING/DrugCombs_v5/", disease, "/BarabasiProx_DrugAdr_", disease, ".rds"))
write.xlsx(proximity_matrix, file = paste0("Analysis/STRING/DrugCombs_v5/", disease, "/BarabasiProx_DrugAdr_", disease, ".xlsx"), overwrite = TRUE)


print(warnings())