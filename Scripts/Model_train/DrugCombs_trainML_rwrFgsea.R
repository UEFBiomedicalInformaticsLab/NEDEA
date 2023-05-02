# Train ML models using RWR-FGSEA results as features 




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
fgsea_result <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))

cat("\n- Feature: Disease2Gene library\n")
Disease2Gene_model <- func_repeated_train(feature_matrix = fgsea_result$NES_Disease2Gene, 
                                          train_test_split = train_test_split, 
                                          data_balance_method = data_balance_method, 
                                          allow_parallel = FALSE, 
                                          nproc = NULL)

cat("\n- Feature: WithdrawalAdr2Gene library\n")
WithdrawalAdr2Gene_model <- func_repeated_train(feature_matrix = fgsea_result$NES_WithdrawalAdr2Gene, 
                                          train_test_split = train_test_split, 
                                          data_balance_method = data_balance_method, 
                                          allow_parallel = FALSE, 
                                          nproc = NULL)

cat("\n- Feature: CombinedDisAdr2Gene library\n")
CombinedDisAdr2Gene_model <- func_repeated_train(feature_matrix = fgsea_result$NES_CombinedDisAdr2Gene, 
                                          train_test_split = train_test_split, 
                                          data_balance_method = data_balance_method, 
                                          allow_parallel = FALSE, 
                                          nproc = NULL)





# Read the FGSEA results (Common library)
fgsea_result <- readRDS(paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))

cat("\n- Feature: keggPath library\n")
keggPath_model <- func_repeated_train(feature_matrix = fgsea_result$NES_keggPath, 
                                      train_test_split = train_test_split, 
                                      data_balance_method = data_balance_method, 
                                      allow_parallel = FALSE, 
                                      nproc = NULL)

cat("\n- Feature: SMPDbPath_DrugMet library\n")
SMPDbPath_DrugMet_model <- func_repeated_train(feature_matrix = fgsea_result$NES_SMPDbPath_DrugMet, 
                                               train_test_split = train_test_split, 
                                               data_balance_method = data_balance_method, 
                                               allow_parallel = FALSE, 
                                               nproc = NULL)

cat("\n- Feature: NES_SMPDbPath_DrugAction library\n")
SMPDbPath_DrugAction_model <- func_repeated_train(feature_matrix = fgsea_result$NES_SMPDbPath_DrugAction, 
                                                      train_test_split = train_test_split, 
                                                      data_balance_method = data_balance_method, 
                                                      allow_parallel = FALSE, 
                                                      nproc = NULL)

cat("\n- Feature: NES_miscGeneSet library\n")
miscGeneSet_model <- func_repeated_train(feature_matrix = fgsea_result$NES_miscGeneSet, 
                                             train_test_split = train_test_split, 
                                             data_balance_method = data_balance_method, 
                                             allow_parallel = FALSE, 
                                             nproc = NULL)





# Save the results
final_results <- list()
final_results$Dis2Gene <- Disease2Gene_model$modelling_results
final_results$WdrlAdr2Gene <- WithdrawalAdr2Gene_model$modelling_results
final_results$CombDisAdr2Gene <- CombinedDisAdr2Gene_model$modelling_results
final_results$keggPath <- keggPath_model$modelling_results
final_results$SMPDbPath_DrugMet <- SMPDbPath_DrugMet_model$modelling_results
final_results$SMPDbPath_DrugAction <- SMPDbPath_DrugAction_model$modelling_results
final_results$miscGeneSet <- miscGeneSet_model$modelling_results
saveRDS(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_rwrFgsea_", disease, ".rds"))


final_results <- list()
final_results$Dis2Gene <- Disease2Gene_model$result_summary_tables
final_results$WdrlAdr2Gene <- WithdrawalAdr2Gene_model$result_summary_tables
final_results$CombDisAdr2Gene <- CombinedDisAdr2Gene_model$result_summary_tables
final_results$keggPath <- keggPath_model$result_summary_tables
final_results$SMPDbPath_DrugMet <- SMPDbPath_DrugMet_model$result_summary_tables
final_results$SMPDbPath_DrugAction <- SMPDbPath_DrugAction_model$result_summary_tables
final_results$miscGeneSet <- miscGeneSet_model$result_summary_tables
final_results <- unlist(final_results, recursive = FALSE)
write.xlsx(final_results, paste0("OutputFiles/Model_train/", disease, "/models_", data_balance_method, "_rwrFgsea_", disease, ".xlsx"), overwrite = TRUE)



stopCluster(cl)
unregister_dopar()
Sys.sleep(60)



print(warnings())