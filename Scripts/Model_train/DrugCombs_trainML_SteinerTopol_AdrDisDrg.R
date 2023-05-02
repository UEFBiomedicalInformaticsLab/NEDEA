# Train ML models using topology of Steiner tree of ADR, disease and drug target related genes



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
train_test_split <- readRDS(paste0("OutputFiles/Model_train/", disease, "/ML_dataSplit_", disease, ".rds"))


# Using external cluster instead of clustering through each function call
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 





# Read the FGSEA results (Efficacy-Safety library)
net_topol <- readRDS(paste0("OutputFiles/Model_train/", disease, "/SteinerTreeTopol_AdrDisDrg_", disease, ".rds"))

features_to_use <- c("number_terminal_nodes", 
                     "number_steiner_nodes",
                     "subNet_diameter",
                     "subNet_averagePathLength",
                     "subNet_clusteringCoefficient",
                     "subNet_density")

net_topol <- net_topol[net_topol$features %in% features_to_use, ]
row.names(net_topol) <- NULL

cat("\n- Feature: Steiner Tree topology (ADR-Disease-DrugTarget)\n")
SteinerTopol_model <- func_repeated_train(feature_matrix = net_topol, 
                                          train_test_split = train_test_split, 
                                          data_balance_method = data_balance_method, 
                                          allow_parallel = FALSE, 
                                          nproc = NULL)





# Save the results
final_results <- list()
final_results$SteinerTopol <- SteinerTopol_model$modelling_results
saveRDS(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_SteinerTopol_AdrDisDrg_", disease, ".rds"))


final_results <- list()
final_results$SteinerTopol <- SteinerTopol_model$result_summary_tables
final_results <- unlist(final_results, recursive = FALSE)
write.xlsx(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_SteinerTopol_AdrDisDrg_", disease, ".xlsx"), overwrite = TRUE)



stopCluster(cl)
unregister_dopar()
Sys.sleep(60)



print(warnings())