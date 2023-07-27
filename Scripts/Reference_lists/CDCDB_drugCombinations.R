set.seed(5081)



# Drug combinations from C-DCDB (https://icc.ise.bgu.ac.il/medical_ai/CDCDB/)



# Load libraries
library(tidyverse)



# Download drug combinations from C-DCDB
if(!dir.exists("Databases/CDCDB/"))dir.create("Databases/CDCDB/")
if(!file.exists("Databases/CDCDB/11.07.2023.zip")){
  warning(paste0("ERROR: C-DCDB file not found !!! \n", 
                 "Download 11.07.2023.zip file manually from https://icc.ise.bgu.ac.il/medical_ai/CDCDB/"))
}





# Extract the drug combinations reported in the FDA Orange Book (RX/OTC type)
if(!file.exists("Databases/CDCDB/orangebook_combs_df.csv")){
  unzip("Databases/CDCDB/11.07.2023.zip", files = "orangebook_combs_df.csv", exdir = "Databases/CDCDB/")
}
CDCDB_drugCombs <- read.csv("Databases/CDCDB/orangebook_combs_df.csv")
CDCDB_drugCombs <- CDCDB_drugCombs[CDCDB_drugCombs$TYPE != "DISCN", ] # Remove discontinued drug combinations
CDCDB_drugCombs <- CDCDB_drugCombs[, "drugbank_ids", drop = FALSE]
CDCDB_drugCombs$drugbank_ids <- gsub('\\"', "", CDCDB_drugCombs$drugbank_ids)
CDCDB_drugCombs$drugbank_ids <- gsub("\\[|\\]", "", CDCDB_drugCombs$drugbank_ids)
CDCDB_drugCombs <- data.frame(str_split(CDCDB_drugCombs$drugbank_ids, ",|, ", simplify = TRUE))
CDCDB_drugCombs <- CDCDB_drugCombs[rowSums(CDCDB_drugCombs != "") == 2, ]
CDCDB_drugCombs <- CDCDB_drugCombs[, 1:2]
colnames(CDCDB_drugCombs) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")
CDCDB_drugCombs <- as.data.frame(apply(CDCDB_drugCombs, 2, function(x) gsub("\\s+", "", x)))
CDCDB_drugCombs <- CDCDB_drugCombs[CDCDB_drugCombs$Drug1_DrugBank_drug_id != CDCDB_drugCombs$Drug2_DrugBank_drug_id,] # remove pairs in which both drugs are same
CDCDB_drugCombs <- CDCDB_drugCombs[(CDCDB_drugCombs$Drug1_DrugBank_drug_id != "-1" | CDCDB_drugCombs$Drug2_DrugBank_drug_id != "-1"), ]                                       
row.names(CDCDB_drugCombs) <- NULL
if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(CDCDB_drugCombs, "InputFiles/ReferenceList/CDCDB_drugCombinations_OrangeBook_rxotc.rds")
                                       
                                       
    

                                       
# Extract the drug combinations reported in the FDA Orange Book (DISCN type)
if(!file.exists("Databases/CDCDB/orangebook_combs_df.csv")){
  unzip("Databases/CDCDB/11.07.2023.zip", files = "orangebook_combs_df.csv", exdir = "Databases/CDCDB/")
}
CDCDB_drugCombs <- read.csv("Databases/CDCDB/orangebook_combs_df.csv")
CDCDB_drugCombs <- CDCDB_drugCombs[CDCDB_drugCombs$TYPE == "DISCN", ] # Remove discontinued drug combinations
CDCDB_drugCombs <- CDCDB_drugCombs[, "drugbank_ids", drop = FALSE]
CDCDB_drugCombs$drugbank_ids <- gsub('\\"', "", CDCDB_drugCombs$drugbank_ids)
CDCDB_drugCombs$drugbank_ids <- gsub("\\[|\\]", "", CDCDB_drugCombs$drugbank_ids)
CDCDB_drugCombs <- data.frame(str_split(CDCDB_drugCombs$drugbank_ids, ",|, ", simplify = TRUE))
CDCDB_drugCombs <- CDCDB_drugCombs[rowSums(CDCDB_drugCombs != "") == 2, ]
CDCDB_drugCombs <- CDCDB_drugCombs[, 1:2]
colnames(CDCDB_drugCombs) <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")
CDCDB_drugCombs <- as.data.frame(apply(CDCDB_drugCombs, 2, function(x) gsub("\\s+", "", x)))
CDCDB_drugCombs <- CDCDB_drugCombs[CDCDB_drugCombs$Drug1_DrugBank_drug_id != CDCDB_drugCombs$Drug2_DrugBank_drug_id,] # remove pairs in which both drugs are same
CDCDB_drugCombs <- CDCDB_drugCombs[(CDCDB_drugCombs$Drug1_DrugBank_drug_id != "-1" | CDCDB_drugCombs$Drug2_DrugBank_drug_id != "-1"), ]                                       
row.names(CDCDB_drugCombs) <- NULL

if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(CDCDB_drugCombs, "InputFiles/ReferenceList/CDCDB_drugCombinations_OrangeBook_discn.rds")




                                       
# Extract conditions from AACT
if(!file.exists("Databases/CDCDB/conditions_df.csv")){
  unzip("Databases/CDCDB/11.07.2023.zip", files = "conditions_df.csv", exdir = "Databases/CDCDB/")
}
conditions <- read.csv("Databases/CDCDB/conditions_df.csv")
conditions <- conditions[grep("cancer|carcinoma|sarcoma", conditions$condition, ignore.case = TRUE), ]

# Read the drug combinations for the selected NCT IDs
if(!file.exists("Databases/CDCDB/design_group_df.csv")){
  unzip("Databases/CDCDB/11.07.2023.zip", files = "design_group_df.csv", exdir = "Databases/CDCDB/")
}
design <- read.csv("Databases/CDCDB/design_group_df.csv")
design <- design[, c("nct_id", "drugbank_identifier")]
design$drugbank_identifier <- gsub('\\"', "", design$drugbank_identifier)
design$drugbank_identifier <- gsub("\\[|\\]", "", design$drugbank_identifier)
design <- design[grep("-1", design$drugbank_identifier, invert = TRUE), ]
tmp1 <- data.frame(str_split(design$drugbank_identifier, ",|, ", simplify = TRUE))
tmp1 <- as.data.frame(apply(tmp1, 2, function(x) gsub("\\s+", "", x)))
tmp1$n_drugs <- rowSums(tmp1 != "") 
design <- cbind(design, tmp1)
design <- design[design$n_drugs == 2, ]
design <- design[, !(colSums(design == "") == nrow(design))]
drugCombs <- design[design$nct_id %in% conditions$nct_id,]
rm(tmp1)

# Add info regarding the conditions to the selected drug combinations
drugCombs$condition <- conditions$condition_downcase[match(drugCombs$nct_id , conditions$nct_id)]
names(drugCombs)[grep("^X", colnames(drugCombs))] <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")                                   
drugCombs <- unique(drugCombs)
row.names(drugCombs) <- NULL

# Extract the trials informations
if(!file.exists("Databases/CDCDB/trials_df.csv")){
  unzip("Databases/CDCDB/11.07.2023.zip", files = "trials_df.csv", exdir = "Databases/CDCDB/")
}
trials <- read.csv("Databases/CDCDB/trials_df.csv")
drugCombs <- merge(drugCombs, trials, by = "nct_id", all.x = TRUE)

if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(drugCombs, "InputFiles/ReferenceList/CDCDB_drugCombinations_AACT.rds")

              


                         
                                                                       
# # Extract conditions from AACT (toxicity ended trials)
# if(!file.exists("Databases/CDCDB/conditions_df.csv")){
#   unzip("Databases/CDCDB/11.07.2023.zip", files = "conditions_df.csv", exdir = "Databases/CDCDB/")
# }
# conditions <- read.csv("Databases/CDCDB/conditions_df.csv")
# conditions <- conditions[grep("cancer", conditions$condition, ignore.case = TRUE), ]
# 
# # Read the drug combinations for the selected NCT IDs
# if(!file.exists("Databases/CDCDB/design_group_df.csv")){
#   unzip("Databases/CDCDB/11.07.2023.zip", files = "design_group_df.csv", exdir = "Databases/CDCDB/")
# }
# design <- read.csv("Databases/CDCDB/design_group_df.csv")
# design <- design[, c("nct_id", "drugbank_identifier")]
# design$drugbank_identifier <- gsub('\\"', "", design$drugbank_identifier)
# design$drugbank_identifier <- gsub("\\[|\\]", "", design$drugbank_identifier)
# design <- design[grep("-1", design$drugbank_identifier, invert = TRUE), ]
# tmp1 <- data.frame(str_split(design$drugbank_identifier, ",|, ", simplify = TRUE))
# tmp1 <- as.data.frame(apply(tmp1, 2, function(x) gsub("\\s+", "", x)))
# tmp1$n_drugs <- rowSums(tmp1 != "") 
# design <- cbind(design, tmp1)
# design <- design[design$n_drugs == 2, ]
# design <- design[, !(colSums(design == "") == nrow(design))]
# drugCombs <- design[design$nct_id %in% conditions$nct_id,]
# rm(tmp1)
# 
# # Extract the trials that were withdrawn/suspended/terminated due to toxicity
# if(!file.exists("Databases/CDCDB/trials_df.csv")){
#   unzip("Databases/CDCDB/11.07.2023.zip", files = "trials_df.csv", exdir = "Databases/CDCDB/")
# }
# trials <- read.csv("Databases/CDCDB/trials_df.csv")
# trials <- trials[trials$overall_status %in% c("Suspended", "Terminated", "Withdrawn"),]
# trials <- trials[grep("toxicity|toxicities", trials$why_stopped, ignore.case = TRUE), ]
# drugCombs <- unique(drugCombs[drugCombs$nct_id %in% trials$nct_id,])
# 
# # Add info regarding the conditions to the selected drug combinations
# drugCombs$condition <- conditions$condition_downcase[match(drugCombs$nct_id , conditions$nct_id)]
# colnames(drugCombs)[grep("^X", colnames(drugCombs))] <- c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")                                       
# row.names(drugCombs) <- NULL
# 
# if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
# saveRDS(drugCombs, "InputFiles/ReferenceList/CDCDB_drugCombinations_AACT_toxic.rds")                                       


                                       
print(warnings())