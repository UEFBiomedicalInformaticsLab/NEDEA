set.seed(5081)
rm(list = ls())


# Drug-Target association from DrugBank


# Load libraries
library(openxlsx)
library(tidyr)
library(stringr)
library(biomaRt)

# Download the complete DrugBank database as XML and parse it using dbparser
# Following that the data is saved as CSV files

if(!dir.exists("Databases/DrugBank")){dir.create("Databases/DrugBank/", recursive = TRUE)} 
if(!file.exists("Databases/DrugBank/drugbank_all_full_database.xml.zip")){
  warning(paste0("ERROR: DrugBank database file not found !!! \n", 
                 "Download file in terminal using:\n", 
                 "\t curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://go.drugbank.com/releases/5-1-8/downloads/all-full-database"))
}

if(!file.exists("Databases/DrugBank/drug.csv")){
  library(dbparser)
  read_drugbank_xml_db("Databases/DrugBank/drugbank_all_full_database.xml.zip")
  targets <- targets(save_csv = TRUE, csv_path = "Databases/DrugBank/",  override_csv = TRUE)
  drugs <- drugs(save_csv = TRUE, csv_path = "Databases/DrugBank/",  override_csv = TRUE)
}


# Extract the biological entities (BE) id from DrugBank for all targets and map to Ensembl Gene ID
# The targets_polypep_ex_ident() function extracts mapping the BE IDs to several external IDs
if(!file.exists("Databases/DrugBank/targets_polypeptides_ext_id.csv")){
  DrugBank_Targets_idMap <- targets_polypep_ex_ident(save_csv = TRUE, csv_path = "Databases/DrugBank/",  override_csv = TRUE)
}
DrugBank_Targets_idMap <- read.csv("Databases/DrugBank/targets_polypeptides_ext_id.csv", header = TRUE)
DrugBank_Targets_idMap <- reshape(as.data.frame(DrugBank_Targets_idMap), idvar = "parent_key", timevar = "resource", direction = "wide")
colnames(DrugBank_Targets_idMap) <- gsub("identifier.", "", colnames(DrugBank_Targets_idMap))
DrugBank_Targets_idMap <- DrugBank_Targets_idMap[grep("_HUMAN$", x = DrugBank_Targets_idMap$`UniProt Accession`),]


# Create mapping from UniProtKB ID to Ensembl GEne ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
uniprotID_2_ensemblID <- getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"), mart = ensembl, filters = "uniprotswissprot", values = DrugBank_Targets_idMap$UniProtKB)


DrugBank_Targets_idMap$ensembl_gene_id <- uniprotID_2_ensemblID$ensembl_gene_id[match(DrugBank_Targets_idMap$UniProtKB, uniprotID_2_ensemblID$uniprotswissprot)]
# All Uniprot Accessions do not map to Ensembl gene ID and leaves NA 


# # Extract the DrugBank primary/parent key to cas number mappting
# DrugBank_Drugs_idMap <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)


# Extract all drug target interactions in human
DrugBank_Drug_Target <- read.csv("Databases/DrugBank/targets.csv", header = TRUE)
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$organism == "Humans", ]
colnames(DrugBank_Drug_Target)[colnames(DrugBank_Drug_Target) == "parent_key"] <- "Node1_drugbank_drug_id"
DrugBank_Drug_Target$Node2_ensembl_gene_id <- DrugBank_Targets_idMap$ensembl_gene_id[match(DrugBank_Drug_Target$id, DrugBank_Targets_idMap$parent_key)]
DrugBank_Drug_Target[DrugBank_Drug_Target == ""] <- NA

DrugBank_Drug_Target_Net <- na.exclude(DrugBank_Drug_Target[, c("Node1_drugbank_drug_id", "Node2_ensembl_gene_id")])
DrugBank_Drug_Target_Net$Node1_type <- "drug"
DrugBank_Drug_Target_Net$Node2_type <- "gene"
DrugBank_Drug_Target_Net$Edge_type <- "undirected"
DrugBank_Drug_Target_Net <- DrugBank_Drug_Target_Net[, c("Node1_drugbank_drug_id", "Node1_type", "Node2_ensembl_gene_id", "Node2_type", "Edge_type")]


if(!dir.exists("InputFiles/Associations/")){dir.create("InputFiles/Associations/", recursive = TRUE)} 
saveRDS(DrugBank_Drug_Target_Net, "InputFiles/Associations/DrugBank_Drug_Target_Net.rds")


print(warnings())