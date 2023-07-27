set.seed(5081)



# External predictions using the CDCDB drug combinations (AACT)

# Notes:
# (1) The two drug combinations from AACT in C-DCDC being used
# (2) Remove drug combinations that contained already in 
#     the training set of the cancer.
# (3) Removed drug combinations that do not contain any reported drug targets
# (4) Run RWR-FGSEA and then use the selected model to predict





# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(openxlsx)
library(tidyverse)
library(caret)
source("Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("Scripts/Functions/Functions_parallelprocesses.R")



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--data_balance_method"), type = "character", default = "none", 
              help = "The method to be used to balance imbalanced data. Possible values: none. Default: none.", metavar = "character"),
  make_option(c("--model"), type = "character", default = NULL, 
              help = "The modelling technique to use. Possible values: glmnet, nb, rf, svmRadial. Default: NULL", metavar = "character"),
  make_option(c("--feature_type"), type = "character", default = NULL, 
              help = "The feature type to use for modelling. Possible values: Disease2Gene, WithdrawalAdr2Gene, CombinedDisAdr2Gene, keggPath, SMPDbPath_DrugMet, SMPDbPath_DrugAction, miscGeneSet. Default: NULL", metavar = "character"),
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

if(is.null(opt$model)){
  print_help(opt_parser)
  stop("--model argument needed", call.=FALSE)
}

if(!opt$model %in% c("glmnet", "knn", "nb", "rf", "svmRadial")){
  stop("--model must be glmnet, knn, nb, rf, svmRadial")
}

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
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
model <- opt$model
feature_type <- opt$feature_type
nproc <- opt$nproc


# Read drug combinations used in training set
drugCombs_training <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", disease, ".rds"))
drugCombs_training <- do.call(rbind, drugCombs_training)
row.names(drugCombs_training) <- NULL
drugs_training <- unique(c(drugCombs_training$Drug1_DrugBank_drug_id, drugCombs_training$Drug2_DrugBank_drug_id))



# Read the AACT drug combinations from CDCDB
drugCombs <- readRDS("InputFiles/ReferenceList/CDCDB_drugCombinations_AACT.rds")
drugCombs <- drugCombs[grep("-1", drugCombs$Drug1_DrugBank_drug_id, invert = TRUE),]
drugCombs <- drugCombs[grep("-1", drugCombs$Drug2_DrugBank_drug_id, invert = TRUE),]
drugCombs <- drugCombs[grep("DBSALT", drugCombs$Drug1_DrugBank_drug_id, invert = TRUE),]
drugCombs <- drugCombs[grep("DBSALT", drugCombs$Drug2_DrugBank_drug_id, invert = TRUE),]


# Remove combinations that contains selected terms
DrugBank_drug_category <- read.csv("Databases/DrugBank/drug_categories.csv")
DrugBank_drug_category <- DrugBank_drug_category[DrugBank_drug_category$category %in% 
                                                   c("Amino Acids", "Coenzymes", "Enzymes and Coenzymes", 
                                                     "Hormones"), ]
drugCombs <- drugCombs[!(drugCombs$Drug1_DrugBank_drug_id %in% DrugBank_drug_category$parent_key | 
                           drugCombs$Drug2_DrugBank_drug_id %in% DrugBank_drug_category$parent_key), ]

# Remove combinations that contain drug pair participating in training
remove_index <- c()
for(i in 1:nrow(drugCombs)){
  drug1 <- drugCombs[i, "Drug1_DrugBank_drug_id"]
  drug2 <- drugCombs[i, "Drug2_DrugBank_drug_id"]
  tmp <- drugCombs_training[(drugCombs_training$Drug1_DrugBank_drug_id == drug1 & drugCombs_training$Drug2_DrugBank_drug_id == drug2) | 
                              (drugCombs_training$Drug1_DrugBank_drug_id == drug2 & drugCombs_training$Drug2_DrugBank_drug_id == drug1), ]
  if(nrow(tmp) > 0){print(i); remove_index <- c(remove_index, i)} 
}
if(length(remove_index) > 0){drugCombs <- drugCombs[-remove_index,]}
rm(tmp)



# Read DrugBank drug annotations
# Keep combinations formed by only small molecular drugs
DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", ] 


drugCombs <- drugCombs[(drugCombs$Drug1_DrugBank_drug_id %in% DrugBank_Drugs$primary_key & 
                          drugCombs$Drug2_DrugBank_drug_id %in% DrugBank_Drugs$primary_key), ]


# Read drug target interactions from drug bank
drug_target_ixn  <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")

# Check if the drugs forming the pairs have reported targets
drugCombs <- drugCombs[(drugCombs$Drug1_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id 
                        & drugCombs$Drug2_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id), ]
drugCombs <- unique(drugCombs)
row.names(drugCombs) <- NULL
print(dim(drugCombs))



# Read drug-drug interactions from DrugBank involving risk or severity for side effects
DrugBank_drugInteractions <- readRDS("InputFiles/ReferenceList/DrugBank_drugInteractions_withRiskSeverity.rds")


# Remove drug pairs that are already reported as adverse pairs in DrugBank
remove_index <- c()
for(i in 1:nrow(drugCombs)){
  drug1 <- drugCombs[i, "Drug1_DrugBank_drug_id"] 
  drug2 <- drugCombs[i, "Drug2_DrugBank_drug_id"] 
  
  tmp <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id == drug1 & DrugBank_drugInteractions$Drug2_DrugBank_drug_id == drug2)|
                                     (DrugBank_drugInteractions$Drug1_DrugBank_drug_id == drug2 & DrugBank_drugInteractions$Drug2_DrugBank_drug_id == drug1),] 
  if(nrow(tmp) > 0){remove_index <- c(remove_index, i)} 
} 
if(length(remove_index) > 0){drugCombs <- drugCombs[-remove_index,]}
rm(tmp)
cat(paste0("\n\nNumber of drug combinations to be tested: ", nrow(drugCombs), "\n\n"))



# Read the network on which to run RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")

# Execute RWR
# nproc <- detectCores()/2 #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 

drugCombs_rwr <- foreach (i=1:nrow(drugCombs),
                          .packages = c("igraph", "tidyr", "dnet")) %dopar% {
                            tmp <- func_dNetRWR_on_drugCombination(rwr_input_network = rwr_input_network,
                                                                   drug1 = drugCombs[i,c("Drug1_DrugBank_drug_id")], 
                                                                   drug2 = drugCombs[i,c("Drug2_DrugBank_drug_id")],
                                                                   drug_target_ixn = drug_target_ixn,
                                                                   rwr_norm_input = "laplacian",
                                                                   rwr_norm_output = "none", 
                                                                   rwr_restart = 0.5)
                            tmp
                          }


for(i in 1:length(drugCombs_rwr)){
  names(drugCombs_rwr)[i] <- paste(drugCombs[i,c("nct_id")], 
                                   drugCombs[i,c("Drug1_DrugBank_drug_id")], 
                                   drugCombs[i,c("Drug2_DrugBank_drug_id")], sep = "__")
}

stopCluster(cl)
unregister_dopar()





# Read the FGSEA enrichment library
switch(feature_type,
       "Disease2Gene" = {
         enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
       }, 
       "WithdrawalAdr2Gene" = {
         enrichment_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
       }, 
       "CombinedDisAdr2Gene" = {
         enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
         names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
         enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
         names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
         enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)
       }, 
       "keggPath" = {
         enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/CHG_keggPath2Gene_lib.rds"))
       }, 
       "SMPDbPath_DrugMet" = {
         enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
         enrichment_lib <- enrichment_lib$`Drug Metabolism`
       }, 
       "SMPDbPath_DrugAction" = {
         enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
         enrichment_lib <- enrichment_lib$`Drug Action`
       }, 
       "miscGeneSet" = {
         enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/miscellaneous_gene_lib.rds"))
         cat("\n\n\n- Executing FGSEA for miscellaneous_gene_lib\n")
       }
)



# Execute FGSEA
enrichment_res <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                              rwr_data = drugCombs_rwr, 
                                              quantile_prob = 0.9,
                                              verbose = FALSE)

drugComb_NES <- func_extract_fgsea_result(enrichment_result = enrichment_res,
                                          result_type = "NES",
                                          enrichment_library = enrichment_lib)

if(!dir.exists("OutputFiles/External_predictions/CDCDB_AACT/features/")){
  dir.create("OutputFiles/External_predictions/CDCDB_AACT/features/", recursive = TRUE)
}  
saveRDS(drugComb_NES, paste0("OutputFiles/External_predictions/CDCDB_AACT/features/CDCDB_AACT__", disease, "_", feature_type, ".rds"))





# Read the selected model for predictions
final_model <- readRDS(paste0("OutputFiles/Final_model/FinalModel_", disease, "_", model, "_", data_balance_method, "_", feature_type, ".rds"))          





# Prepare the input data
data <- as.data.frame(t(drugComb_NES))
data <- data[, colnames(data) %in% colnames(final_model$trainingData)]
preProcess <- preProcess(data, method = c("center", "scale"))
data <- predict(object = preProcess, newdata = data)

# Predicted
predictions <- predict(object = final_model, newdata = data, type = "prob") 
print(predictions)
row.names(predictions) <- row.names(data)
predictions$predicted_class <- predict(object = final_model, newdata = data, type = "raw")
predictions <- rownames_to_column(predictions, "drugCombs")
predictions <- separate(data = predictions, 
                        col = "drugCombs", 
                        into = c("nct_id", "Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id"), 
                        remove = FALSE)
predictions$Drug1_name <- DrugBank_Drugs$name[match(predictions$Drug1_DrugBank_drug_id, DrugBank_Drugs$primary_key)]
predictions$Drug2_name <- DrugBank_Drugs$name[match(predictions$Drug2_DrugBank_drug_id, DrugBank_Drugs$primary_key)]
predictions <- predictions[order(predictions$Eff, decreasing = TRUE),]
predictions <- merge(predictions, drugCombs, by = c("nct_id", "Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id"), all.x = TRUE)

predictions <- predictions %>% select(c("drugCombs", "nct_id", 
                               "Drug1_DrugBank_drug_id", "Drug1_name", 
                               "Drug2_DrugBank_drug_id", "Drug2_name", 
                               "Adv", "Eff", "predicted_class"), everything())

predictions <- predictions[, !colnames(predictions) %in% c("drugbank_identifier", "n_drugs")]

# Save as file
if(!dir.exists("OutputFiles/External_predictions/CDCDB_AACT/")){
  dir.create("OutputFiles/External_predictions/CDCDB_AACT/", recursive = TRUE)
}                          
write.xlsx(predictions, paste0("OutputFiles/External_predictions/CDCDB_AACT/ExtPred_CDCDB_AACT__", disease, "_", model, "_", data_balance_method, "_", feature_type, ".xlsx"), 
           rowNames = FALSE, overwrite = TRUE)



print(warnings())