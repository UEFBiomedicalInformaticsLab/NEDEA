set.seed(5081)
rm(list = ls())


# Drug-Target association from Therapeutic Target Database (TTD)


# Load libraries
library(openxlsx)
library(tidyverse)
# source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_ID_Conversion.R")


# Read drug - target mapping file
if(!dir.exists("Databases/TTD/")){dir.create("Databases/TTD/", recursive = TRUE)} 
if(!file.exists("Databases/TTD/P1-07-Drug-TargetMapping.xlsx")){
  download.file(url = "http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-07-Drug-TargetMapping.xlsx",
                destfile = "Databases/TTD/P1-07-Drug-TargetMapping.xlsx", method = "wget")
}
TTD_Drug_Gene <- read.xlsx("Databases/TTD/P1-07-Drug-TargetMapping.xlsx")




# Map target ID to Ensembl gene ID
if(!dir.exists("Databases/TTD/")){dir.create("Databases/TTD/", recursive = TRUE)} 
if(!file.exists("Databases/TTD/P2-01-TTD_uniprot_all.txt")){
  download.file(url = "http://db.idrblab.net/ttd/sites/default/files/ttd_database/P2-01-TTD_uniprot_all.txt",
                destfile = "Databases/TTD/P2-01-TTD_uniprot_all.txt", method = "wget")
}
TTD_Targets <- read.table("Databases/TTD/P2-01-TTD_uniprot_all.txt", sep = "\t", skip = 22, fill = TRUE, quote = "")
ttdTargetId_2_uniprotEntryName <- reshape(as.data.frame(TTD_Targets), idvar = "V1", timevar = "V2", direction = "wide")
ttdTargetId_2_uniprotEntryName <- as.data.frame(separate_rows(ttdTargetId_2_uniprotEntryName, V3.UNIPROID, sep = "; "))[,-1]
colnames(ttdTargetId_2_uniprotEntryName) <- c("ttd_target_id", "uniprot_entry_name", "target_name", "target_type")
ttdTargetId_2_uniprotEntryName <- ttdTargetId_2_uniprotEntryName[(str_detect(ttdTargetId_2_uniprotEntryName$uniprot_entry_name, "HUMAN$")),] # Retain only targets in human


tmp <- func_uniprot_mapping(ttdTargetId_2_uniprotEntryName$uniprot_entry_name, "ID", "ENSEMBL_ID")
ttdTargetId_2_ensemblGeneId <- merge(ttdTargetId_2_uniprotEntryName, tmp, by.x = "uniprot_entry_name", by.y = "From", all = TRUE)
colnames(ttdTargetId_2_ensemblGeneId)[5] <- "Target_ensembl_gene_id"
rm(tmp)
TTD_Drug_Gene <- TTD_Drug_Gene[(TTD_Drug_Gene$TargetID %in% ttdTargetId_2_ensemblGeneId$ttd_target_id),] # Retain drug-target association only for human targets
TTD_Drug_Gene <- merge(TTD_Drug_Gene, ttdTargetId_2_ensemblGeneId[, c("ttd_target_id", "Target_ensembl_gene_id")], by.x = "TargetID", by.y = "ttd_target_id")