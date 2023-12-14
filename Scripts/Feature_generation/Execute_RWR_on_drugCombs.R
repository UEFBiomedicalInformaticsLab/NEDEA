set.seed(5081)


# Script to execute RWR on the drug combinations


# Load libraries
library(unixtools)
library(optparse)
library(igraph)
library(org.Hs.eg.db)
library(dnet)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known", 
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer.", call.=FALSE)
  }
}


# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type
nproc <- opt$nproc


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))


# Read the extended drug target data
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))


# Select the column containing the drug target based on th user input
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


# Extract the targets for the drug combinations
seed_matrix_targets <- drugCombs_targets[,drug_target_col, drop = FALSE]
seed_matrix_targets <- apply(seed_matrix_targets, 1, function(x){paste(x, collapse = ",")})


# Convert the target sets to a matrix, each column corresponds to a drug pair
rwr_seed_matrix <- do.call(cbind, lapply(seed_matrix_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys = target_set, 
                                     columns = "ENSEMBL", 
                                     keytype = "SYMBOL"))
  sapply(V(rwr_input_network)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

colnames(rwr_seed_matrix) <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")


# Run RWR using dRWR

if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 

rwr_result <- dRWR(g = rwr_input_network, 
                   setSeeds = rwr_seed_matrix,
                   normalise = "row",
                   restart = 0.5, 
                   normalise.affinity.matrix = "none",
                   multicores = nproc)

colnames(rwr_result) <- colnames(rwr_seed_matrix)
rownames(rwr_result) <- rownames(rwr_seed_matrix)


stopCluster(cl)
unregister_dopar()



if(!dir.exists("OutputFiles/RWR_results/")){
  dir.create("OutputFiles/RWR_results/", recursive = TRUE)
} 
saveRDS(rwr_result, paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds"))



print(warnings())