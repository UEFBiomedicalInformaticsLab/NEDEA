set.seed(5081)


# Script to calculate network proximity to the efficacy and safety library


# Load libraries
library(unixtools)
library(optparse)
library(igraph)
library(foreach)
library(doParallel)
library(tidyverse)
source("Scripts/Functions/Functions_proximity.R")
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
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call. = FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer.", call.=FALSE)
  }
}

if(is.null(opt$nproc)){
  opt$nproc <- detectCores()/2
  }


# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type
nproc <- opt$nproc

cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


# Read the network on which to calculate proximity
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))



cl <- makeCluster(nproc)
registerDoParallel(cl)


proximity_result_final <- list()

# Proximity with efficacy library
cat("\n--- Calculating proximity with efficacy library")
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib) <- paste0("[DISEASE] ", names(enrichment_lib))

proximity_result <- func_calculate_proximity(input_network = input_network, 
                                             enrichment_library = enrichment_lib,
                                             disease = disease,
                                             drug_target_type = drug_target_type,
                                             nproc = nproc)

proximity_result_final[["efficacy"]] <- proximity_result


# Proximity with safety library
cat("\n--- Calculating proximity with safety library")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
names(enrichment_lib) <- paste0("[ADR] ", names(enrichment_lib))

proximity_result <- func_calculate_proximity(input_network = input_network, 
                                             enrichment_library = enrichment_lib,
                                             disease = disease,
                                             drug_target_type = drug_target_type,
                                             nproc = nproc)

proximity_result_final[["safety"]] <- proximity_result


# Proximity with combined efficacy-safety library
cat("\n--- Calculating proximity with combined efficacy-safety library")
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

proximity_result <- func_calculate_proximity(input_network = input_network, 
                                             enrichment_library = enrichment_lib,
                                             disease = disease,
                                             drug_target_type = drug_target_type,
                                             nproc = nproc)

proximity_result_final[["combinedEfficacySafety"]] <- proximity_result

stopCluster(cl)
unregister_dopar()

if(!dir.exists("OutputFiles/Network_proximity_results/")){
  dir.create("OutputFiles/Network_proximity_results/", recursive = TRUE)
} 
saveRDS(proximity_result_final, paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))



print(warnings())