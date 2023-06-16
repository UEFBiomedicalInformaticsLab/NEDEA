set.seed(5081)



# Train ML models using Barabasi metrics between targets of the two drugs and the ADR related genes



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(openxlsx)
library(tidyverse)
source("Scripts/Functions/Functions_modelTraning_drugComb.R")
source("Scripts/Functions/Functions_parallelprocesses.R")



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: none. Default: none.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
  
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$data_balance_method %in% c("none")){
  stop("--data_balance_method: No data balancing methods currently included. Default: none")
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}


# Define global options for this script 
disease <- opt$disease
data_balance_method <- opt$data_balance_method
nproc <- opt$nproc

cat(paste0("\n\nTraining model for: ", disease, "\n"))
cat(paste0("\nData balance method: ", data_balance_method, "\n"))



# Read the train-test split
train_test_split <- readRDS(paste0("OutputFiles/Model_train/", disease, "/ML_dataSplit_", disease, ".rds"))


# Using external cluster instead of clustering through each function call
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 





# Read the Barabasi metrics for drug-ADR genes
proximity_matrix <- readRDS(paste0("OutputFiles/Model_train/", disease, "/BarabasiProx_DrugAdr_", disease, ".rds"))



cat("\n- Feature: all Barabasi Proximities (Drug-ADR)\n")
proximity_all <- bind_rows(proximity_matrix, .id = "proximity")
proximity_all$features <- paste0("[", proximity_all$proximity, "] ", proximity_all$features)
proximity_all <- proximity_all[, -1]
proximity_all[sapply(proximity_all, is.infinite)] <- NA
proximity_all <- proximity_all[!rowSums(is.na(proximity_all)) == ncol(proximity_all)-1 ,]
proximity_all <- proximity_all[, !colSums(is.na(proximity_all)) > 0]
BarabasiProx_all_model <- func_train_model(feature_matrix = proximity_all, 
                                           train_test_split = train_test_split, 
                                           data_balance_method = data_balance_method, 
                                           allow_parallel = FALSE, 
                                           nproc = NULL)


cat("\n- Feature: Barabasi Proximity (Closest) (Drug-ADR)\n")
proximity_closest <- proximity_matrix[["proximity_closest"]]
proximity_closest[sapply(proximity_closest, is.infinite)] <- NA
proximity_closest <- proximity_closest[!rowSums(is.na(proximity_closest)) == ncol(proximity_closest)-1 ,]
proximity_closest <- proximity_closest[, !colSums(is.na(proximity_closest)) > 0]
BarabasiProx_closest_model <- func_train_model(feature_matrix = proximity_closest, 
                                               train_test_split = train_test_split, 
                                               data_balance_method = data_balance_method, 
                                               allow_parallel = FALSE, 
                                               nproc = NULL)


cat("\n- Feature: Barabasi Proximity (Shortest) (Drug-ADR)\n")
proximity_shortest <- proximity_matrix[["proximity_shortest"]]
proximity_shortest[sapply(proximity_shortest, is.infinite)] <- NA
proximity_shortest <- proximity_shortest[!rowSums(is.na(proximity_shortest)) == ncol(proximity_shortest)-1 ,]
proximity_shortest <- proximity_shortest[, !colSums(is.na(proximity_shortest)) > 0]
BarabasiProx_shortest_model <- func_train_model(feature_matrix = proximity_shortest, 
                                                train_test_split = train_test_split, 
                                                data_balance_method = data_balance_method, 
                                                allow_parallel = FALSE, 
                                                nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_centre) (Drug-ADR)\n")
proximity_centre <- proximity_matrix[["proximity_centre"]]
proximity_centre[sapply(proximity_centre, is.infinite)] <- NA
proximity_centre <- proximity_centre[!rowSums(is.na(proximity_centre)) == ncol(proximity_centre)-1 ,]
proximity_centre <- proximity_centre[, !colSums(is.na(proximity_centre)) > 0]
BarabasiProx_centre_model <- func_train_model(feature_matrix = proximity_centre, 
                                              train_test_split = train_test_split, 
                                              data_balance_method = data_balance_method, 
                                              allow_parallel = FALSE, 
                                              nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_kernel) (Drug-ADR)\n")
proximity_kernel <- proximity_matrix[["proximity_kernel"]]
proximity_kernel[sapply(proximity_kernel, is.infinite)] <- NA
proximity_kernel <- proximity_kernel[!rowSums(is.na(proximity_kernel)) == ncol(proximity_kernel)-1 ,]
proximity_kernel <- proximity_kernel[, !colSums(is.na(proximity_kernel)) > 0]
BarabasiProx_kernel_model <- func_train_model(feature_matrix = proximity_kernel, 
                                              train_test_split = train_test_split, 
                                              data_balance_method = data_balance_method, 
                                              allow_parallel = FALSE, 
                                              nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_separation) (Drug-ADR)\n")
proximity_separation <- proximity_matrix[["proximity_separation"]]
proximity_separation[sapply(proximity_separation, is.infinite)] <- NA
proximity_separation <- proximity_separation[!rowSums(is.na(proximity_separation)) == ncol(proximity_separation)-1 ,]
proximity_separation <- proximity_separation[, !colSums(is.na(proximity_separation)) > 0]
BarabasiProx_separation_model <- func_train_model(feature_matrix = proximity_separation, 
                                                  train_test_split = train_test_split, 
                                                  data_balance_method = data_balance_method, 
                                                  allow_parallel = FALSE, 
                                                  nproc = NULL)




# Save the results
final_results <- list()
final_results$BbsiProx_all <- BarabasiProx_all_model$modelling_results
final_results$BbsiProx_closest <- BarabasiProx_closest_model$modelling_results
final_results$BbsiProx_shortest <- BarabasiProx_shortest_model$modelling_results
final_results$BbsiProx_centre <- BarabasiProx_centre_model$modelling_results
final_results$BbsiProx_kernel <- BarabasiProx_kernel_model$modelling_results
final_results$BbsiProx_separation <- BarabasiProx_separation_model$modelling_results

saveRDS(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_BarabasiProx_DrugAdr_", disease, ".rds"))


final_results <- list()
final_results$BbsiProx <- BarabasiProx_all_model$result_summary_tables
final_results$BbsiProx_closest <- BarabasiProx_closest_model$result_summary_tables
final_results$BbsiProx_shortest <- BarabasiProx_shortest_model$result_summary_tables
final_results$BbsiProx_centre <- BarabasiProx_centre_model$result_summary_tables
final_results$BbsiProx_kernel <- BarabasiProx_kernel_model$result_summary_tables
final_results$BbsiProx_separation <- BarabasiProx_separation_model$result_summary_tables
final_results <- unlist(final_results, recursive = FALSE)
write.xlsx(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_BarabasiProx_DrugAdr_", disease, ".xlsx"), overwrite = TRUE)



stopCluster(cl)
unregister_dopar()
Sys.sleep(60)



print(warnings())