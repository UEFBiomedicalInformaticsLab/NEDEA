set.seed(5081)



# Drug-Disease association from OpenTargets
# Notes: 
# (1) Only small molecular drug combinations



# Load libraries
library(tidyverse)
library(sparklyr)
library(sparklyr.nested)



# spark_install()
# Retrieve all drugs information from OpenTargets
if(!dir.exists("Databases/OpenTargets/")){dir.create("Databases/OpenTargets/", recursive = TRUE)}
if(!dir.exists("Databases/OpenTargets/molecule")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/parquet/molecule
")
}

sc <- spark_connect(master = "local")
OpenTargets_Drugs <- spark_read_parquet(sc, "OpenTargets_Drugs", "Databases/OpenTargets/molecule/")

# List columns in the tibble
columns <- OpenTargets_Drugs %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()


# Extract the selected columns
OpenTargets_Drug_Disease_Net <- OpenTargets_Drugs %>%
  dplyr::select(c("id", "inchiKey", "drugType", "name", "maximumClinicalTrialPhase", "hasBeenWithdrawn", "isApproved", "tradeNames", "linkedDiseases")) %>%
  collect()


# Simplify to a data frame
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net %>% unnest(linkedDiseases) %>% 
  group_by(id) %>% 
  mutate(col=seq_along(id)) %>% 
  spread(key=col, value=linkedDiseases) %>%
  unnest("1") %>% unnest(tradeNames)

OpenTargets_Drug_Disease_Net <- as.data.frame(OpenTargets_Drug_Disease_Net)
names(OpenTargets_Drug_Disease_Net)[names(OpenTargets_Drug_Disease_Net) == "1"] <- "linkedDiseases"
names(OpenTargets_Drug_Disease_Net)[names(OpenTargets_Drug_Disease_Net) == "2"] <- "linkedDiseasesCount"
names(OpenTargets_Drug_Disease_Net)[names(OpenTargets_Drug_Disease_Net) == "drugType"] <- "Drug_type"

spark_disconnect(sc)

# Retain only small molecule drugs
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net[OpenTargets_Drug_Disease_Net$Drug_type %in% "Small molecule", ]

# Add column for clinical status
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net %>%
  mutate(Drug_Disease_clinical_status = case_when(
    (isApproved == TRUE & hasBeenWithdrawn == TRUE) ~ "Withdrawn",
    (isApproved == TRUE & hasBeenWithdrawn == FALSE) ~ "Approved",
    (isApproved == FALSE & hasBeenWithdrawn == TRUE) ~ "Withdrawn",
    (isApproved == FALSE & hasBeenWithdrawn == FALSE) ~ paste0("Phase ", maximumClinicalTrialPhase)
  ))


# Map drugs to DrugBankID using ChEMBL IDs
DrugBank_Drugs_idMap <- read.csv("Databases/DrugBank/drug_external_identifiers.csv", header = TRUE)
DrugBank_Drugs_idMap <- DrugBank_Drugs_idMap[DrugBank_Drugs_idMap$resource == "ChEMBL", ]
OpenTargets_Drug_Disease_Net$DrugBank_drug_id <- DrugBank_Drugs_idMap$parent_key[match(OpenTargets_Drug_Disease_Net$id, DrugBank_Drugs_idMap$identifier)]


OpenTargets_Drug_Disease_Net <- unique(OpenTargets_Drug_Disease_Net[, c("DrugBank_drug_id", "linkedDiseases", "Drug_Disease_clinical_status", "Drug_type")])
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net[!is.na(OpenTargets_Drug_Disease_Net$DrugBank_drug_id), ]
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net[!is.na(OpenTargets_Drug_Disease_Net$linkedDiseases), ]

colnames(OpenTargets_Drug_Disease_Net)[1:2] <- c("Node1_drugbank_drug_id", "Node2_disease_id")
OpenTargets_Drug_Disease_Net$Node1_type <- "drug"
OpenTargets_Drug_Disease_Net$Node2_type <- "disease"
OpenTargets_Drug_Disease_Net$Edge_type <- "undirected"
OpenTargets_Drug_Disease_Net <- OpenTargets_Drug_Disease_Net %>%
  select(c("Node1_drugbank_drug_id", "Node1_type", "Node2_disease_id", "Node2_type", "Edge_type"), everything())
rownames(OpenTargets_Drug_Disease_Net) <- NULL

if(!dir.exists("InputFiles/Associations/")){dir.create("InputFiles/Associations/", recursive = TRUE)} 
saveRDS(OpenTargets_Drug_Disease_Net, "InputFiles/Associations/OpenTargets_Drug_Disease_Net.rds")



print(warnings())