set.seed(5081)


# Extract the drug-target associations from DrugBank


# Load libraries
library(openxlsx)
library(tidyverse)
library(biomaRt)
library(dbparser)


#####


# Download the complete DrugBank database as XML and parse it using dbparser

if(!dir.exists("Databases/DrugBank")){dir.create("Databases/DrugBank/", recursive = TRUE)} 
if(!file.exists("Databases/DrugBank/drugbank_all_full_database.xml.zip")){ 
  warning(paste0("ERROR: DrugBank database file not found !!! \n", 
                 "Download file in terminal using:\n", 
                 "\t curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://go.drugbank.com/releases/5-1-10/downloads/all-full-database"))
}

if(!file.exists("Databases/DrugBank/parsed_DrugBank_data.rds")){
  require(dbparser)
  dvobj <- parseDrugBank(db_path = "Databases/DrugBank/drugbank_all_full_database.xml.zip",
                         drug_options = drug_node_options(),
                         parse_salts = TRUE,
                         parse_products = TRUE,
                         references_options = references_node_options(),
                         cett_options = cett_nodes_options())
  saveRDS(dvobj, "Databases/DrugBank/parsed_DrugBank_data.rds")
}


#####


# Extract the biological entities (BE) id from DrugBank for all targets and map to Ensembl Gene ID
DrugBank_Targets_idMap <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_Targets_idMap <- DrugBank_Targets_idMap$cett$targets$polypeptides$external_identy
DrugBank_Targets_idMap <- reshape(as.data.frame(DrugBank_Targets_idMap), idvar = "parent_key", timevar = "resource", direction = "wide")
colnames(DrugBank_Targets_idMap) <- gsub("identifier.", "", colnames(DrugBank_Targets_idMap))
DrugBank_Targets_idMap <- DrugBank_Targets_idMap[grep("_HUMAN$", x = DrugBank_Targets_idMap$`UniProt Accession`),]


# Create mapping from UniProtKB ID to Ensembl GEne ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
uniprotID_2_ensemblID <- getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"), mart = ensembl, filters = "uniprotswissprot", values = DrugBank_Targets_idMap$UniProtKB)


DrugBank_Targets_idMap$ensembl_gene_id <- uniprotID_2_ensemblID$ensembl_gene_id[match(DrugBank_Targets_idMap$UniProtKB, uniprotID_2_ensemblID$uniprotswissprot)]
# All Uniprot Accessions do not map to Ensembl gene ID and leaves NA 


# Extract all drug target interactions in human
DrugBank_Drug_Target <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_Drug_Target <- DrugBank_Drug_Target$cett$targets$general_information
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$organism == "Humans", ]
colnames(DrugBank_Drug_Target)[colnames(DrugBank_Drug_Target) == "parent_key"] <- "drugbank_drug_id"
DrugBank_Drug_Target$ensembl_gene_id <- DrugBank_Targets_idMap$ensembl_gene_id[match(DrugBank_Drug_Target$id, DrugBank_Targets_idMap$parent_key)]
DrugBank_Drug_Target[DrugBank_Drug_Target == ""] <- NA

DrugBank_Drug_Target_Net <- na.exclude(DrugBank_Drug_Target[, c("drugbank_drug_id", "ensembl_gene_id")])


if(!dir.exists("InputFiles/Associations/")){dir.create("InputFiles/Associations/", recursive = TRUE)} 
saveRDS(DrugBank_Drug_Target_Net, "InputFiles/Associations/DrugBank_Drug_Target_associations.rds")


#####


print(warnings())