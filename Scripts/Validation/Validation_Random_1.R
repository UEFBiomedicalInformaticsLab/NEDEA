# Validation using randomly generated drug combinations
#
# Notes:
# (1) The effective and adverse drug combinations were generated such that
#      the targets of the drugs in each category were specific to either disease 
#      related genes or ADR related genes
# (2) To ensure that the genes/targets in the two category do not overlap,
#     only the genes with near highest distances considered.
# (3) Filtering done to remove combinations which inlcude atleast one drug in the 
#     training set and has known targets
# (4) Since the number of possible combinations were too many, we randomly selected n combinations





# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(caret)
library(openxlsx)
library(tidyverse)
library(MLmetrics)
library(igraph)
source("Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("Scripts/Functions/Functions_parallelprocesses.R")








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
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc





# Read drug combinations used in training set
drugCombs_training <- list()
for(tmp_dis in c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")){
  tmp1 <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", tmp_dis, ".rds"))
  tmp1 <- do.call(rbind, tmp1)
  row.names(tmp1) <- NULL
  drugCombs_training[[tmp_dis]] <- tmp1
}
drugCombs_training <- do.call(rbind, drugCombs_training)
row.names(drugCombs_training) <- NULL
drugCombs_training <- unique(drugCombs_training)
drugs_training <- unique(c(drugCombs_training$Drug1_DrugBank_drug_id, drugCombs_training$Drug2_DrugBank_drug_id))
rm(list = c("tmp_dis", "tmp1"))





if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 



cat(paste0("\nExecuting for: ", disease, "\n"))


# Read the efficacy and safety related genes
enrichment_lib_dis <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_dis) <- paste0("[DISEASE] ", names(enrichment_lib_dis))
enrichment_lib_dis <- unique(unlist(enrichment_lib_dis, use.names = FALSE))

enrichment_lib_adr <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
names(enrichment_lib_adr) <- paste0("[ADR] ", names(enrichment_lib_adr))
enrichment_lib_adr <- unique(unlist(enrichment_lib_adr, use.names = FALSE))


# Read the network on which RWR will be run
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")


# Calculate the distances between the disease related genes and the ADR related genes
# Select gene pairs with the among highest distance
distance_table <- data.frame()
for(gene in enrichment_lib_dis){
  sp <- as.data.frame(distances(graph = rwr_input_network,
                        v = V(rwr_input_network)[V(rwr_input_network)$name %in% gene],
                        to = V(rwr_input_network)[V(rwr_input_network)$name %in% enrichment_lib_adr],
                        mode = "all",
                        algorithm = "dijkstra"))
  distance_table <- rbind(distance_table, sp)
}
distance_table <- rownames_to_column(distance_table, "dis")
distance_table <- reshape(distance_table,
                            direction = "long",
                            varying = colnames(distance_table)[colnames(distance_table) != "dis"],
                            v.names = "sp",
                            timevar = "adr",
                            times = colnames(distance_table)[colnames(distance_table) != "dis"])

row.names(distance_table) <- NULL

distance_table <- distance_table[distance_table$sp >= max(distance_table$sp)-2,]


cat(paste0("\nDisease genes: ", length(unique(distance_table$dis))))
cat(paste0("\nADR genes: ", length(unique(distance_table$adr))))


# Read drug target interactions from drug bank
drug_target_ixn  <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")


# Generate combinations based on drugs targeting disease genes
dis_drugs <- unique(drug_target_ixn[drug_target_ixn$Node2_ensembl_gene_id %in% unique(distance_table$dis), "Node1_drugbank_drug_id"])
cat(paste0("\nDisease drugs: ", length(dis_drugs)))
dis_drugCombs <- expand.grid(dis_drugs, dis_drugs)
dis_drugCombs <- dis_drugCombs[dis_drugCombs$Var1 != dis_drugCombs$Var2, ]
colnames(dis_drugCombs) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")
dis_drugCombs$Class <- "Eff"

# Generate combinations based on drugs targeting ADR genes
adr_drugs <- unique(drug_target_ixn[drug_target_ixn$Node2_ensembl_gene_id %in% unique(distance_table$adr), "Node1_drugbank_drug_id"])
cat(paste0("\ADR drugs: ", length(adr_drugs)))
adr_drugCombs <- expand.grid(adr_drugs, adr_drugs)
adr_drugCombs <- adr_drugCombs[adr_drugCombs$Var1 != adr_drugCombs$Var2, ]
colnames(adr_drugCombs) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")
adr_drugCombs$Class <- "Adv"


# Randomly select n combinations from each category for analysis
dis_drugCombs <- dis_drugCombs[sample(1:nrow(dis_drugCombs), 100, replace = FALSE),]
adr_drugCombs <- adr_drugCombs[sample(1:nrow(adr_drugCombs), 100, replace = FALSE),]
drugCombs <- rbind(dis_drugCombs, adr_drugCombs)


# Remove combinations that contain atleast one drug participating in training
remove_index <- c()
for(i in 1:nrow(drugCombs)){
  drug1 <- drugCombs[i, "Drug1_DrugBank_drug_id"]
  drug2 <- drugCombs[i, "Drug2_DrugBank_drug_id"]
  tmp <- drugCombs_training[drugCombs_training$Drug1_DrugBank_drug_id == drug1 & drugCombs_training$Drug2_DrugBank_drug_id == drug2, ]
  if(nrow(tmp) > 0){remove_index <- c(remove_index, i)} 
}
if(length(remove_index) > 0){drugCombs <- drugCombs[-remove_index,]}
rm(tmp)


# Read drug target interactions from drug bank
drug_target_ixn  <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")


# Check if the drugs forming the pairs have reported targets
drugCombs <- drugCombs[(drugCombs$Drug1_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id 
                                                & drugCombs$Drug2_DrugBank_drug_id %in% drug_target_ixn$Node1_drugbank_drug_id), ]
drugCombs <- unique(drugCombs)
drugCombs <- data.frame(lapply(drugCombs, as.character), stringsAsFactors = FALSE)
row.names(drugCombs) <- NULL


# Execute RWR
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
  names(drugCombs_rwr)[i] <- paste(drugCombs[i,c("Class")], 
                                   drugCombs[i,c("Drug1_DrugBank_drug_id")], 
                                   drugCombs[i,c("Drug2_DrugBank_drug_id")], sep = "__")
}


# Run FGSEA on combined Disease2Gene and drug withdrawal based Adr2Gene library (efficacy + safety)
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))

enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

enrichment_res <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                        rwr_data = drugCombs_rwr, 
                                                        quantile_prob = 0.9,
                                                        verbose = FALSE)

drugComb_NES <- func_extract_fgsea_result(enrichment_result = enrichment_res,
                                               result_type = "NES",
                                               enrichment_library = enrichment_lib)

if(!dir.exists("OutputFiles/Validation/Random_1/")){
  dir.create("OutputFiles/Validation/Random_1/", recursive = TRUE)
} 
write.csv(drugComb_NES, paste0("OutputFiles/Validation/Random_1/drugComb_NES_", disease, ".csv"))


# Read the selected model
model <- readRDS(paste0("OutputFiles/Final_model/FinalModel_", disease, "_rf_none_CombinedDisAdr2Gene.rds"))          

# Prepare the input data
data <- as.data.frame(t(drugComb_NES))
# preProcess <- preProcess(data, method = c("zv", "center", "scale"))
# data <- predict(object = preProcess, newdata = data)

# Predict
predictions <- predict(object = model, newdata = data, type = "prob") 
predictions$Class <- predict(object = model, newdata = data, type = "raw")
predictions <- rownames_to_column(predictions, "drugCombs")
predictions <- predictions[order(predictions$Eff, decreasing = TRUE),]
predictions$actual_class <- substr(predictions$drugCombs, 1, 3)
validation_predictions <- list()
validation_predictions[["predictions"]] <- predictions

confusionMatrix <- confusionMatrix(table(predictions$Class, predictions$actual_class), positive = "Eff")
rocauc <- AUC(y_pred = predictions$Eff, y_true = ifelse(predictions$actual_class == "Eff", 1, 0))
prauc <- PRAUC(y_pred = predictions$Eff, y_true = ifelse(predictions$actual_class == "Eff", 1, 0))

validation_predictions[["stats"]] <- data.frame(PositiveClass = confusionMatrix$positive,
                                                ROCAUC = as.numeric(rocauc),
                                                PRAUC = as.numeric(prauc),                           
                                                Sensitivity =  unname (confusionMatrix$byClass["Sensitivity"]),
                                                Specificity =  unname (confusionMatrix$byClass["Specificity"]),
                                                Precision =  unname (confusionMatrix$byClass["Precision"]),
                                                Recall =  unname (confusionMatrix$byClass["Recall"]),
                                                F1 =  unname (confusionMatrix$byClass["F1"]),
                                                BalancedAccuracy =  unname (unname (confusionMatrix$byClass["Balanced Accuracy"])))

write.xlsx(validation_predictions, paste0("OutputFiles/Validation/Random_1/validation_", disease, ".xlsx"))





stopCluster(cl)
unregister_dopar()


print(warnings())