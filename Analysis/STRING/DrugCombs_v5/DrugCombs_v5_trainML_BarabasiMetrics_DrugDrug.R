rm(list = ls())



# Train ML models using Barabasi metrics between drug-pairs as features (version 5)



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(openxlsx)
source("Scripts/Functions/Functions_modelTraning_drugComb.R")
source("Scripts/Functions/Functions_parallelprocesses.R")




# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: SMOTE, downSample, upSample, or none. Default: none.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
  
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$data_balance_method %in% c("SMOTE", "downSample", "upSample", "none")){
  stop("--data_balance_method must be SMOTE, downSample, upSample or none")
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
train_test_split <- readRDS(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/ML_dataSplit_", disease, ".rds"))


# Using external cluster instead of clustering through each function call
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 





# Read the Barabasi metrics for drug pairs
proximity_matrix <- readRDS(paste0("Analysis/STRING/DrugCombs_v5/", disease, "/BarabasiProx_DrugDrug_", disease, ".rds"))
proximity_matrix[proximity_matrix == "Inf"] <- NaN
proximity_matrix <- proximity_matrix[, !colSums(proximity_matrix == "NaN")>0]




cat("\n- Feature: all Barabasi Proximities (Drug-Drug)\n")
BarabasiProx_all_model <- func_repeated_train(feature_matrix = proximity_matrix, 
                                             train_test_split = train_test_split, 
                                             data_balance_method = data_balance_method, 
                                             allow_parallel = FALSE, 
                                             nproc = NULL)


cat("\n- Feature: Barabasi Proximity (Closest) (Drug-Drug)\n")
proximity_closest <- proximity_matrix[proximity_matrix$features == "proximity_closest", ]
BarabasiProx_closest_model <- func_repeated_train(feature_matrix = proximity_closest, 
                                              train_test_split = train_test_split, 
                                              data_balance_method = data_balance_method, 
                                              allow_parallel = FALSE, 
                                              nproc = NULL)


cat("\n- Feature: Barabasi Proximity (Shortest) (Drug-Drug)\n")
proximity_shortest <- proximity_matrix[proximity_matrix$features == "proximity_shortest", ]
BarabasiProx_shortest_model <- func_repeated_train(feature_matrix = proximity_shortest, 
                                              train_test_split = train_test_split, 
                                              data_balance_method = data_balance_method, 
                                              allow_parallel = FALSE, 
                                              nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_centre) (Drug-Drug)\n")
proximity_centre <- proximity_matrix[proximity_matrix$features == "proximity_centre", ]
BarabasiProx_centre_model <- func_repeated_train(feature_matrix = proximity_centre, 
                                                   train_test_split = train_test_split, 
                                                   data_balance_method = data_balance_method, 
                                                   allow_parallel = FALSE, 
                                                   nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_kernel) (Drug-Drug)\n")
proximity_kernel <- proximity_matrix[proximity_matrix$features == "proximity_kernel", ]
BarabasiProx_kernel_model <- func_repeated_train(feature_matrix = proximity_kernel, 
                                                 train_test_split = train_test_split, 
                                                 data_balance_method = data_balance_method, 
                                                 allow_parallel = FALSE, 
                                                 nproc = NULL)


cat("\n- Feature: Barabasi Proximity (proximity_separation) (Drug-Drug)\n")
proximity_separation <- proximity_matrix[proximity_matrix$features == "proximity_separation", ]
BarabasiProx_separation_model <- func_repeated_train(feature_matrix = proximity_separation, 
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

saveRDS(final_results, paste0("Analysis/STRING/DrugCombs_v5/", disease, "/models_", data_balance_method, "_BarabasiProx_DrugDrug_", disease, ".rds"))


final_results <- list()
final_results$BbsiProx <- BarabasiProx_all_model$result_summary_tables
final_results$BbsiProx_closest <- BarabasiProx_closest_model$result_summary_tables
final_results$BbsiProx_shortest <- BarabasiProx_shortest_model$result_summary_tables
final_results$BbsiProx_centre <- BarabasiProx_centre_model$result_summary_tables
final_results$BbsiProx_kernel <- BarabasiProx_kernel_model$result_summary_tables
final_results$BbsiProx_separation <- BarabasiProx_separation_model$result_summary_tables
final_results <- unlist(final_results, recursive = FALSE)
write.xlsx(final_results, paste0("Analysis/STRING/DrugCombs_v5/", disease, "/models_", data_balance_method, "_BarabasiProx_DrugDrug_", disease, ".xlsx"), overwrite = TRUE)



stopCluster(cl)
unregister_dopar()
Sys.sleep(60)

print(warnings())