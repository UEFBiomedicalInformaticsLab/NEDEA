set.seed(5081)



# External predictions using the CDD drug combinations 



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



# Download cancer drug list from Cancer Drugs Database
if(!dir.exists("Databases/CancerDrugsDatabase/"))dir.create("Databases/CancerDrugsDatabase/", recursive = TRUE)
if(!file.exists("Databases/CancerDrugsDatabase/cancerdrugsdb.txt")){
  download.file(url = "https://acfdata.coworks.be/cancerdrugsdb.txt",
                destfile = "Databases/CancerDrugsDatabase/cancerdrugsdb.txt", method = "wget")
}


# Read the  list of cancer drugs from cancer drugs database
CDD_drugs <- read.table("https://acfdata.coworks.be/cancerdrugsdb.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)
CDD_drugs$DrugBank.ID <- gsub("^<a href=\"https://go.drugbank.com/drugs/[A-D0-9]+\">", "", CDD_drugs$DrugBank.ID)
CDD_drugs$DrugBank.ID <- gsub("</a>", "", CDD_drugs$DrugBank.ID)
if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/", recursive = TRUE)
saveRDS(CDD_drugs, "InputFiles/ReferenceList/CDD_drugs.rds")                                       


# Remove drugs used in training 
CDD_drugs <- CDD_drugs[!(CDD_drugs$DrugBank.ID %in% drugs_training),]



# Read drug target interactions from drug bank
drug_target_ixn  <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")

# Check if the drugs forming the pairs have reported targets
CDD_drugs <- CDD_drugs[CDD_drugs$DrugBank.ID %in% drug_target_ixn$Node1_drugbank_drug_id, ]



# Generate all possible combinations
drugCombs <- as.data.frame(t(combn(CDD_drugs$DrugBank.ID, 2, simplify = TRUE)))
colnames(drugCombs) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")



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
  names(drugCombs_rwr)[i] <- paste("Unk", drugCombs[i,c("Drug1_DrugBank_drug_id")], 
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

if(!dir.exists("OutputFiles/External_predictions/CDD/features/")){
  dir.create("OutputFiles/External_predictions/CDD/features/", recursive = TRUE)
}  
saveRDS(drugComb_NES, paste0("OutputFiles/External_predictions/CDD/features/CDD_drugCombs__", disease, "_", feature_type, ".rds"))





# Read the selected model for predictions
final_model <- readRDS(paste0("OutputFiles/Final_model/FinalModel_", disease, "_", model, "_", data_balance_method, "_", feature_type, ".rds"))          





# Prepare the input data
data <- as.data.frame(t(drugComb_NES))
data <- data[, colnames(data) %in% colnames(final_model$trainingData)]
preProcess <- preProcess(data, method = c("center", "scale"))
data <- predict(object = preProcess, newdata = data)

# Predicted
predictions <- predict(object = final_model, newdata = data, type = "prob") 
# print(predictions)
row.names(predictions) <- row.names(data)
predictions$predicted_class <- predict(object = final_model, newdata = data, type = "raw")
predictions <- rownames_to_column(predictions, "drugCombs")
predictions <- separate(data = predictions, 
                        col = "drugCombs", 
                        into = c("actual_class", "Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id"), 
                        remove = FALSE)
predictions$Drug1_name <- CDD_drugs$Product[match(predictions$Drug1_DrugBank_drug_id, CDD_drugs$DrugBank.ID)]
predictions$Drug2_name <- CDD_drugs$Product[match(predictions$Drug2_DrugBank_drug_id, CDD_drugs$DrugBank.ID)]
predictions <- predictions[, c("drugCombs", "actual_class", 
                               "Drug1_DrugBank_drug_id", "Drug1_name", 
                               "Drug2_DrugBank_drug_id", "Drug2_name", 
                               "Adv", "Eff", "predicted_class")]
predictions <- predictions[order(predictions$Eff, decreasing = TRUE),]

predictions$Drug1_indications <- CDD_drugs$Indications[match(predictions$Drug1_DrugBank_drug_id, CDD_drugs$DrugBank.ID)]
predictions$Drug2_indications <- CDD_drugs$Indications[match(predictions$Drug2_DrugBank_drug_id, CDD_drugs$DrugBank.ID)]


# Save as file
if(!dir.exists("OutputFiles/External_predictions/CDD/")){
  dir.create("OutputFiles/External_predictions/CDD/", recursive = TRUE)
}                          
write.xlsx(predictions, paste0("OutputFiles/External_predictions/CDD/ExtPred_CDD_drugCombs__", disease, "_", model, "_", data_balance_method, "_", feature_type, ".xlsx"), 
           rowNames = FALSE, overwrite = TRUE)



print(warnings())