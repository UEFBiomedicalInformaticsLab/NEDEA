set.seed(5081)
rm(list = ls())



# Drug combinations from DrugCombDb (http://drugcombdb.denglab.org/main)



# Drug combinations from DrugCombDb [Two drug combinations from assay]
if(!dir.exists("Databases/DrugCombDb/"))dir.create("Databases/DrugCombDb/")
if(!file.exists("Databases/DrugCombDb/Syner&Antag_voting.csv")){
  download.file(url = "http://drugcombdb.denglab.org/download/Syner&Antag_voting.csv",
                destfile = "Databases/DrugCombDb/Syner&Antag_voting.csv", method = "wget")
}

DrugCombDb_drugCombCat <- read.csv("Databases/DrugCombDb/Syner&Antag_voting.csv", header = TRUE)



# Map drugs to DrugBankID using PubChem CIDs
if(!dir.exists("Databases/DrugCombDb/"))dir.create("Databases/DrugCombDb/")
if(!file.exists("Databases/DrugCombDb/drug_chemical_info.csv")){
  download.file(url = "http://drugcombdb.denglab.org/download/drug_chemical_info.csv",
                destfile = "Databases/DrugCombDb/drug_chemical_info.csv", method = "wget")
}

DrugCombDb_drugInfo <- read.csv("Databases/DrugCombDb/drug_chemical_info.csv", header = TRUE)
DrugCombDb_drugInfo$PubChem_cid <- gsub("^CIDs", "", DrugCombDb_drugInfo$cIds)
DrugCombDb_drugInfo$PubChem_cid <- gsub("^0+", "", DrugCombDb_drugInfo$PubChem_cid)

DrugBank_Drugs_idMap <- read.csv("Databases/DrugBank/drug_external_identifiers.csv", header = TRUE)
DrugBank_Drugs_idMap <- DrugBank_Drugs_idMap[DrugBank_Drugs_idMap$resource == "PubChem Compound", ]
DrugCombDb_drugInfo$DrugBank_drug_id <- DrugBank_Drugs_idMap$parent_key[match(DrugCombDb_drugInfo$PubChem_cid, DrugBank_Drugs_idMap$identifier)]

DrugCombDb_drugCombCat$Drug1_DrugBank_drug_id <- DrugCombDb_drugInfo$DrugBank_drug_id[match(DrugCombDb_drugCombCat$Drug1, DrugCombDb_drugInfo$drugName)]
DrugCombDb_drugCombCat$Drug2_DrugBank_drug_id <- DrugCombDb_drugInfo$DrugBank_drug_id[match(DrugCombDb_drugCombCat$Drug2, DrugCombDb_drugInfo$drugName)]


DrugCombDb_drugCombCat_all <- unique(na.exclude(DrugCombDb_drugCombCat[, c("Drug1", "Drug2_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))

DrugCombDb_drugComb <- list()
DrugCombDb_drugComb$synergy <- DrugCombDb_drugCombCat_all[DrugCombDb_drugCombCat_all$classification == "synergy", ]
DrugCombDb_drugComb$antagonism <- DrugCombDb_drugCombCat_all[DrugCombDb_drugCombCat_all$classification == "antagonism", ]
row.names(DrugCombDb_drugComb$synergy) <- row.names(DrugCombDb_drugComb$antagonism) <- NULL
if(!dir.exists("InputFiles/ReferenceList/")){dir.create("InputFiles/ReferenceList/", recursive = TRUE)} 
saveRDS(DrugCombDb_drugComb, "InputFiles/ReferenceList/DrugCombDb_drugCombinations.rds")
print("All drug combinations")
print(lapply(DrugCombDb_drugComb, nrow))




# Lung cancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
LungCancer_cellLines = c("EKVX", "HOP-62", "HOP-92", "NCI-H226", "NCI-H322M", 
                         "NCI-H460", "NCI-H522", "A549", "A427", "NCIH1650",
                         "NCIH2122", "NCIH520", "SKMES1")

DrugCombDb_drugCombCat_LungCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% LungCancer_cellLines), ] 
DrugCombDb_drugCombCat_LungCancer <- unique(na.exclude(DrugCombDb_drugCombCat_LungCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_LungCancer_drugComb <- list()
DrugCombDb_LungCancer_drugComb$synergy <- DrugCombDb_drugCombCat_LungCancer[DrugCombDb_drugCombCat_LungCancer$classification == "synergy", ]
DrugCombDb_LungCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_LungCancer[DrugCombDb_drugCombCat_LungCancer$classification == "antagonism", ]
row.names(DrugCombDb_LungCancer_drugComb$synergy) <- row.names(DrugCombDb_LungCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_LungCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_LungCancer_drugCombinations.rds")
print("LungCancer")
print(lapply(DrugCombDb_LungCancer_drugComb, nrow))





# Leukemia specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
Leukemia_cellLines = c("KBM-7")

DrugCombDb_drugCombCat_Leukemia <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% Leukemia_cellLines), ] 
DrugCombDb_drugCombCat_Leukemia <- unique(na.exclude(DrugCombDb_drugCombCat_Leukemia[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_Leukemia_drugComb <- list()
DrugCombDb_Leukemia_drugComb$synergy <- DrugCombDb_drugCombCat_Leukemia[DrugCombDb_drugCombCat_Leukemia$classification == "synergy", ]
DrugCombDb_Leukemia_drugComb$antagonism <- DrugCombDb_drugCombCat_Leukemia[DrugCombDb_drugCombCat_Leukemia$classification == "antagonism", ]
row.names(DrugCombDb_Leukemia_drugComb$synergy) <- row.names(DrugCombDb_Leukemia_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_Leukemia_drugComb, "InputFiles/ReferenceList/DrugCombDb_Leukemia_drugCombinations.rds")
print("Leukemia")
print(lapply(DrugCombDb_Leukemia_drugComb, nrow))





# BreastCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
BreastCancer_cellLines = c("MCF7", "UACC-257", "MDA-MB-468", "T-47D", "BT-549", 
                           "HS 578T", "KPL1", "EFM192B", "MDAMB436", "OCUBM", 
                           "MDA-MB-231",  "ZR751")

DrugCombDb_drugCombCat_BreastCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% BreastCancer_cellLines), ] 
DrugCombDb_drugCombCat_BreastCancer <- unique(na.exclude(DrugCombDb_drugCombCat_BreastCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_BreastCancer_drugComb <- list()
DrugCombDb_BreastCancer_drugComb$synergy <- DrugCombDb_drugCombCat_BreastCancer[DrugCombDb_drugCombCat_BreastCancer$classification == "synergy", ]
DrugCombDb_BreastCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_BreastCancer[DrugCombDb_drugCombCat_BreastCancer$classification == "antagonism", ]
row.names(DrugCombDb_BreastCancer_drugComb$synergy) <- row.names(DrugCombDb_BreastCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_BreastCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_BreastCancer_drugCombinations.rds")
print("BreastCancer")
print(lapply(DrugCombDb_BreastCancer_drugComb, nrow))





# ProstateCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
ProstateCancer_cellLines = c("DU-145", "PC-3", "LNCAP", "VCAP")

DrugCombDb_drugCombCat_ProstateCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% ProstateCancer_cellLines), ] 
DrugCombDb_drugCombCat_ProstateCancer <- unique(na.exclude(DrugCombDb_drugCombCat_ProstateCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_ProstateCancer_drugComb <- list()
DrugCombDb_ProstateCancer_drugComb$synergy <- DrugCombDb_drugCombCat_ProstateCancer[DrugCombDb_drugCombCat_ProstateCancer$classification == "synergy", ]
DrugCombDb_ProstateCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_ProstateCancer[DrugCombDb_drugCombCat_ProstateCancer$classification == "antagonism", ]
row.names(DrugCombDb_ProstateCancer_drugComb$synergy) <- row.names(DrugCombDb_ProstateCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_ProstateCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_ProstateCancer_drugCombinations.rds")
print("ProstateCancer")
print(lapply(DrugCombDb_ProstateCancer_drugComb, nrow))





# ColonCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
ColonCancer_cellLines = c("HCT116", "LOVO", "RKO", "SW620", "COLO 205", 
                          "HCC-2998", "HCT-15", "HCT-116", "HT29", 
                          "KM12", "SW-620", "COLO320DM", "DLD1")

DrugCombDb_drugCombCat_ColonCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% ColonCancer_cellLines), ] 
DrugCombDb_drugCombCat_ColonCancer <- unique(na.exclude(DrugCombDb_drugCombCat_ColonCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_ColonCancer_drugComb <- list()
DrugCombDb_ColonCancer_drugComb$synergy <- DrugCombDb_drugCombCat_ColonCancer[DrugCombDb_drugCombCat_ColonCancer$classification == "synergy", ]
DrugCombDb_ColonCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_ColonCancer[DrugCombDb_drugCombCat_ColonCancer$classification == "antagonism", ]
row.names(DrugCombDb_ColonCancer_drugComb$synergy) <- row.names(DrugCombDb_ColonCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_ColonCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_ColonCancer_drugCombinations.rds")
print("ColonCancer")
print(lapply(DrugCombDb_ColonCancer_drugComb, nrow))





# LiverCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
LiverCancer_cellLines = c("Huh-7", "Mak", "DD2", "3D7", "HB3")

DrugCombDb_drugCombCat_LiverCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% LiverCancer_cellLines), ] 
DrugCombDb_drugCombCat_LiverCancer <- unique(na.exclude(DrugCombDb_drugCombCat_LiverCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_LiverCancer_drugComb <- list()
DrugCombDb_LiverCancer_drugComb$synergy <- DrugCombDb_drugCombCat_LiverCancer[DrugCombDb_drugCombCat_LiverCancer$classification == "synergy", ]
DrugCombDb_LiverCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_LiverCancer[DrugCombDb_drugCombCat_LiverCancer$classification == "antagonism", ]
row.names(DrugCombDb_LiverCancer_drugComb$synergy) <- row.names(DrugCombDb_LiverCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_LiverCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_LiverCancer_drugCombinations.rds")
print("LiverCancer")
print(lapply(DrugCombDb_LiverCancer_drugComb, nrow))





# OvaryCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
OvaryCancer_cellLines = c("IGROV1", "NCI\\\\/ADR-RES", "OVCAR-4", "OVCAR-8", 
                          "SK-OV-3", "OVCAR-5", "A2780", "CAOV3", 
                          "ES2", "OV90", "OVCAR3", "PA1", "UWB1289", 
                          "UWB1289+BRCA1")

DrugCombDb_drugCombCat_OvaryCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% OvaryCancer_cellLines), ] 
DrugCombDb_drugCombCat_OvaryCancer <- unique(na.exclude(DrugCombDb_drugCombCat_OvaryCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_OvaryCancer_drugComb <- list()
DrugCombDb_OvaryCancer_drugComb$synergy <- DrugCombDb_drugCombCat_OvaryCancer[DrugCombDb_drugCombCat_OvaryCancer$classification == "synergy", ]
DrugCombDb_OvaryCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_OvaryCancer[DrugCombDb_drugCombCat_OvaryCancer$classification == "antagonism", ]
row.names(DrugCombDb_OvaryCancer_drugComb$synergy) <- row.names(DrugCombDb_OvaryCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_OvaryCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_OvaryCancer_drugCombinations.rds")
print("OvaryCancer")
print(lapply(DrugCombDb_OvaryCancer_drugComb, nrow))





# SkinCancer specific
# Cell lines retrieved from: https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm
SkinCancer_cellLines = c("LOX IMVI", "MALME-3M", "M14", "MDA-MB-435",
                         "SK-MEL-5", "SK-MEL-2", "SK-MEL-28", 
                         "A2058", "A375", "HT144", "RPMI7951", "SKMEL30", 
                         "UACC62", "MMAC-SF")

DrugCombDb_drugCombCat_SkinCancer <- DrugCombDb_drugCombCat[(DrugCombDb_drugCombCat$Cell.line %in% SkinCancer_cellLines), ] 
DrugCombDb_drugCombCat_SkinCancer <- unique(na.exclude(DrugCombDb_drugCombCat_SkinCancer[, c("Drug1", "Drug1_DrugBank_drug_id", "Drug2", "Drug2_DrugBank_drug_id", "classification")]))
DrugCombDb_SkinCancer_drugComb <- list()
DrugCombDb_SkinCancer_drugComb$synergy <- DrugCombDb_drugCombCat_SkinCancer[DrugCombDb_drugCombCat_SkinCancer$classification == "synergy", ]
DrugCombDb_SkinCancer_drugComb$antagonism <- DrugCombDb_drugCombCat_SkinCancer[DrugCombDb_drugCombCat_SkinCancer$classification == "antagonism", ]
row.names(DrugCombDb_SkinCancer_drugComb$synergy) <- row.names(DrugCombDb_SkinCancer_drugComb$antagonism) <- NULL
saveRDS(DrugCombDb_SkinCancer_drugComb, "InputFiles/ReferenceList/DrugCombDb_SkinCancer_drugCombinations.rds")
print("SkinCancer")
print(lapply(DrugCombDb_SkinCancer_drugComb, nrow))


print(warnings())