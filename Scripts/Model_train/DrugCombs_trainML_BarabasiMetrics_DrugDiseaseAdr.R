# Train ML models using Barabasi proximity between targets of the two drugs, the disease related genes
# and the ADR related genes combined





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
              help = "The method to be used to balance imbalanced data. Possible values: SMOTE, downSample, upSample, or none. Default: none.", metavar = "character"),
  make_option(c("--proximity"), type = "character", default = "separation", 
              help = "The proximity type to use. Possible values: closest, shortest, centre, kernel, separation. Default: separation", metavar = "character"),
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

if(!opt$proximity %in% c("closest", "shortest", "centre", "kernel", "separation")){
  stop("--proximity must be closest, shortest, centre, kernel, or separation")
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
proximity <- opt$proximity
nproc <- opt$nproc

cat(paste0("\n\nTraining model for: ", disease, "\n"))
cat(paste0("\nData balance method: ", data_balance_method, "\n"))
cat(paste0("\nProximity: ", proximity, "\n"))
proximity <- paste0("proximity_", proximity)



# Read the train-test split
train_test_split <- readRDS(paste0("OutputFiles/Model_train/", disease, "/ML_dataSplit_", disease, ".rds"))


# Using external cluster instead of clustering through each function call
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 



# Read the Barabasi metrics 
proximity_matrix_DrgDis <- readRDS(paste0("OutputFiles/Model_train/", disease, "/BarabasiProx_DrugDisease_", disease, ".rds"))
proximity_matrix_DrgDis <- proximity_matrix_DrgDis[[proximity]]
proximity_matrix_DrgDis$features <- paste0("[DrgDis] ", proximity_matrix_DrgDis$features)

proximity_matrix_DrgAdr <- readRDS(paste0("OutputFiles/Model_train/", disease, "/BarabasiProx_DrugAdr_", disease, ".rds"))
proximity_matrix_DrgAdr <- proximity_matrix_DrgAdr[[proximity]]
proximity_matrix_DrgAdr$features <- paste0("[DrgAdr] ", proximity_matrix_DrgAdr$features)

# proximity_matrix_DrgDrg <- readRDS(paste0("OutputFiles/Model_train/", disease, "/BarabasiProx_DrugDrug_", disease, ".rds"))
# proximity_matrix_DrgDrg <- proximity_matrix_DrgDrg[proximity_matrix_DrgDrg$features == proximity, ]
# proximity_matrix_DrgDrg$features <- gsub(pattern = proximity, replacement = "[DrgDrg] Proximity between targets of drugs" , x = proximity_matrix_DrgDrg$features)


# proximity_matrix <- do.call(rbind, list(proximity_matrix_DrgDis, proximity_matrix_DrgAdr, proximity_matrix_DrgDrg))
proximity_matrix <- do.call(rbind, list(proximity_matrix_DrgDis, proximity_matrix_DrgAdr))
# proximity_matrix[proximity_matrix == "Inf"] <- NA
proximity_matrix[sapply(proximity_matrix, is.infinite)] <- NA
proximity_matrix <- proximity_matrix[!rowSums(is.na(proximity_matrix)) == ncol(proximity_matrix)-1 ,]
proximity_matrix <- proximity_matrix[, !colSums(is.na(proximity_matrix)) > 0]


# Train model
BarabasiProx_model <- func_repeated_train(feature_matrix = proximity_matrix, 
                                              train_test_split = train_test_split, 
                                              data_balance_method = data_balance_method, 
                                              allow_parallel = FALSE, 
                                              nproc = NULL)

proximity <- gsub("^proximity_", "", proximity)


# Save the results
final_results <- list()
final_results$BbsiProx <- BarabasiProx_model$modelling_results
saveRDS(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_BarabasiProx_DrgDisAdr_", disease, "_", proximity, ".rds"))


final_results <- list()
final_results$BbsiProx <- BarabasiProx_model$result_summary_tables
final_results <- unlist(final_results, recursive = FALSE)
write.xlsx(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_BarabasiProx_DrgDisAdr_", disease, "_", proximity, ".xlsx"), overwrite = TRUE)



stopCluster(cl)
unregister_dopar()
Sys.sleep(60)



print(warnings())