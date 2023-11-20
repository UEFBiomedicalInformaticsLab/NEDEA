set.seed(5081)


# Script to extend the efficacy and safety gene sets by RWR


# Load libraries
library(unixtools)
library(optparse)
library(igraph)
library(dnet)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")
source("Scripts/Functions/Functions_RWR.R")
# source("Scripts/Functions/Functions_FGSEA.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


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


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))


# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))



# Read the libraries and extract gene list
enrichment_lib_efficacy <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
enrichment_lib_efficacy <- unique(unlist(enrichment_lib_efficacy, use.names = FALSE))

enrichment_lib_safety <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
enrichment_lib_safety <- unique(unlist(enrichment_lib_safety, use.names = FALSE))


# prepare the seed matrix for RWR
rwr_seed_matrix  <- do.call(cbind, lapply(list("Efficacy" = enrichment_lib_efficacy,
                                               "Safety" = enrichment_lib_safety), 
                                          function(x){
                                            sapply(V(rwr_input_network)$name, function(y) ifelse(y %in% x, 1, 0))
                                          }))


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


# Calculate the RWR probability threshold
rwr_threshold <- sapply(apply(rwr_result, 2, func_RWR_threshold), function(x) x$ELB)


# Extract the genes above the selected threshold
extended_enrichment_lib <- list()
for(col_select in colnames(rwr_result)){
  rwr_data_select <- rwr_result[, col_select]
  rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold[col_select]], decreasing = TRUE)
  extended_enrichment_lib[[col_select]] <- names(rankedGeneList)
}




# Save the libraries
if(!dir.exists("Extend_Eff_Saf/Enrichment_analysis_libraries/")){dir.create("Extend_Eff_Saf/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(extended_enrichment_lib, paste0("Extend_Eff_Saf/Enrichment_analysis_libraries/", disease, "_extended_EfficacySafety_lib.rds"))


print(warnings())