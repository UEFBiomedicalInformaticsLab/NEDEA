set.seed(5081)



# Script to compile results from external predictions using the CDCDB drug combinations (RX/OTC type)
# Notes:
# (a) In the final output, checks if the drugs are already approved for certain type of cancers.


# Load libraries
library(openxlsx)
library(tidyverse)
library(sparklyr)
library(sparklyr.nested)
library(unixtools)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Compile the results into one file
files <- list.files(path = "OutputFiles/External_predictions/CDCDB",
                    pattern = "^ExtPred_CDCDB_OrangeBook_rxotc__", 
                    full.names = TRUE)


extPred_drugCombs <- data.frame()
for(file in files){
  disease <- strsplit(x = file, split = "__")[[1]][[2]]
  disease <- strsplit(x = disease, split = "_")[[1]][[1]]
  
  model <- strsplit(x = file, split = "__")[[1]][[2]]
  model <- strsplit(x = model, split = "_")[[1]][[2]]
  
  balance <- strsplit(x = file, split = "__")[[1]][[2]]
  balance <- strsplit(x = balance, split = "_")[[1]][[3]]
  
  featureType <- strsplit(x = file, split = "__")[[1]][[2]]
  featureType <- strsplit(x = featureType, split = "_")[[1]][[4]]
  featureType <- strsplit(x = featureType, split = "\\.")[[1]][[1]]
  
  
  tmp1 <- read.xlsx(file)
  tmp1$disease <- disease
  tmp1$model <- model
  tmp1$balance <- balance
  tmp1$featureType <- featureType
  
  extPred_drugCombs <- rbind(extPred_drugCombs, tmp1)
}
rm(tmp1)



# Retrieve all drugs information from OpenTargets
if(!dir.exists("Databases/OpenTargets/")){dir.create("Databases/OpenTargets/", recursive = TRUE)}
if(!dir.exists("Databases/OpenTargets/indication")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/parquet/indication")
}
if(!dir.exists("Databases/OpenTargets/diseases")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/parquet/diseases")
}

spark_install()
sc <- spark_connect(master = "local")


# Simplify to a data frame
OpenTargets_Drug_indications <- spark_read_parquet(sc, "OpenTargets_Drug_indications", "Databases/OpenTargets/indication/")

OpenTargets_Drug_approvedIndications <- as.data.frame(sdf_unnest_longer(OpenTargets_Drug_indications, "approvedIndications"))
OpenTargets_Drug_approvedIndications <- OpenTargets_Drug_approvedIndications[, c("id", "approvedIndications")]

OpenTargets_Drug_indications <- OpenTargets_Drug_indications %>% sdf_unnest(indications) 
OpenTargets_Drug_indications <- as.data.frame(OpenTargets_Drug_indications)
OpenTargets_Drug_indications <- OpenTargets_Drug_indications[, c("id", "disease", "efoName", "maxPhaseForIndication")]
OpenTargets_Drug_indications <- OpenTargets_Drug_indications[grep("cancer|carcinoma|sarcoma|leukemia", 
                                                                  OpenTargets_Drug_indications$efoName,
                                                                  ignore.case = TRUE), ]

# Map drugs to DrugBankID using ChEMBL IDs
DrugBank_Drugs_idMap <- read.csv("Databases/DrugBank/drug_external_identifiers.csv", header = TRUE)
DrugBank_Drugs_idMap <- DrugBank_Drugs_idMap[DrugBank_Drugs_idMap$resource == "ChEMBL", ]
OpenTargets_Drug_indications$DrugBank_drug_id <- DrugBank_Drugs_idMap$parent_key[match(OpenTargets_Drug_indications$id, DrugBank_Drugs_idMap$identifier)]
OpenTargets_Drug_approvedIndications$DrugBank_drug_id <- DrugBank_Drugs_idMap$parent_key[match(OpenTargets_Drug_approvedIndications$id, DrugBank_Drugs_idMap$identifier)]


# Map disease names
OpenTargets_Diseases <- spark_read_parquet(sc, "OpenTargets_Diseases", "Databases/OpenTargets/diseases/")
OpenTargets_Diseases <- as.data.frame(OpenTargets_Diseases)
OpenTargets_Diseases <- OpenTargets_Diseases[, c("id", "name")]
OpenTargets_Drug_approvedIndications$efoName <- OpenTargets_Diseases$name[match(OpenTargets_Drug_approvedIndications$approvedIndications, OpenTargets_Diseases$id)]
OpenTargets_Drug_approvedIndications <- OpenTargets_Drug_approvedIndications[grep("cancer|carcinoma|sarcoma|leukemia", 
                                                                                  OpenTargets_Drug_approvedIndications$efoName, 
                                                                                  ignore.case = TRUE), ]



# Add annotations
extPred_drugCombs$Drug1_indication <- NA
extPred_drugCombs$Drug2_indication <- NA
extPred_drugCombs$Drug1_maxPhaseForIndication <- NA
extPred_drugCombs$Drug2_maxPhaseForIndication <- NA
extPred_drugCombs$Drug1_approvedIndication <- NA
extPred_drugCombs$Drug2_approvedIndication <- NA


for(i in 1:nrow(extPred_drugCombs)){
  Drug1 <- extPred_drugCombs[i, "Drug1_DrugBank_drug_id"]
  Drug2 <- extPred_drugCombs[i, "Drug2_DrugBank_drug_id"]
  disease <- extPred_drugCombs[i, "disease"]
  
  
  tmp1 <- OpenTargets_Drug_indications[(OpenTargets_Drug_indications$DrugBank_drug_id %in% Drug1), ]
  
  tmp2 <- OpenTargets_Drug_indications[(OpenTargets_Drug_indications$DrugBank_drug_id %in% Drug2), ]
  
  
  tmp3 <- OpenTargets_Drug_approvedIndications[(OpenTargets_Drug_approvedIndications$DrugBank_drug_id %in% Drug1), ]
  
  tmp4 <- OpenTargets_Drug_approvedIndications[(OpenTargets_Drug_approvedIndications$DrugBank_drug_id %in% Drug2), ]
  
  extPred_drugCombs[i, "Drug1_indication"] <- ifelse(nrow(tmp1) > 0, paste(unique(tmp1$efoName), collapse = ", "), NA)
  extPred_drugCombs[i, "Drug2_indication"] <- ifelse(nrow(tmp2) > 0, paste(unique(tmp2$efoName), collapse = ", "), NA)
  
  extPred_drugCombs[i, "Drug1_maxPhaseForIndication"] <- ifelse(nrow(tmp1) > 0, max(tmp1$maxPhaseForIndication), NA)
  extPred_drugCombs[i, "Drug2_maxPhaseForIndication"] <- ifelse(nrow(tmp2) > 0, max(tmp2$maxPhaseForIndication), NA)
  
  extPred_drugCombs[i, "Drug1_approvedIndication"] <- ifelse(nrow(tmp1) > 0, paste(unique(tmp3$efoName), collapse = ", "), NA)
  extPred_drugCombs[i, "Drug2_approvedIndication"] <- ifelse(nrow(tmp2) > 0, paste(unique(tmp4$efoName), collapse = ", "), NA)
}



# Save as file
if(!dir.exists("OutputFiles/External_predictions/CDCDB/")){
  dir.create("OutputFiles/External_predictions/CDCDB/", recursive = TRUE)
}                          
write.xlsx(extPred_drugCombs, "OutputFiles/External_predictions/CDCDB/Compiled_ExtPred_CDCDB_OrangeBook_rxotc.xlsx", 
           rowNames = FALSE, overwrite = TRUE)

