# Drug combinations from FIMM DrugComb (https://drugcomb.fimm.fi/)





# Load libraries
library(unixtools)
source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_ID_Conversion.R")





# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")





# Download drug combinations from FIMM DrugComb 
if(!dir.exists("Databases/FimmDrugComb/"))dir.create("Databases/FimmDrugComb/")
if(!file.exists("Databases/FimmDrugComb/summary_v_1_5.csv")){
  download.file(url = "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv",
                destfile = "Databases/FimmDrugComb/summary_v_1_5.csv", method = "wget")
}

# Read the DrugComb data
FimmDrugComb_drugCombCat <- read.csv("Databases/FimmDrugComb/summary_v_1_5.csv", header = TRUE)


# Extract list of participating drugs
FimmDrugComb_drugs <- sort(unique(c(FimmDrugComb_drugCombCat$drug_row, FimmDrugComb_drugCombCat$drug_col))) #4622


FimmDrugComb_drugs_byStudy <-list()
for(study in sort(unique(FimmDrugComb_drugCombCat$study_name))){
  tmp <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$study_name == study, ]
  FimmDrugComb_drugs_byStudy[[study]] <- sort(unique(c(tmp$drug_row, tmp$drug_col)))
}


# sort(table(substr(FimmDrugComb_drugs, 1, 3)), decreasing = TRUE)



# Check if drug names have been previously mapped. 
# Saves time in retrieving
if(file.exists("InputFiles/ReferenceList/FimmDrugComb_drugs_idMap.rds")){
  FimmDrugComb_drugs_idMap_old <- readRDS("InputFiles/ReferenceList/FimmDrugComb_drugs_idMap.rds")
  FimmDrugComb_drugs <- FimmDrugComb_drugs[FimmDrugComb_drugs %in% FimmDrugComb_drugs_idMap$drug_names]
}else{FimmDrugComb_drugs_idMap_old <- data.frame()}





# Map using DrugBank supplied identifiers
DrugBank_Drugs_idMap <- read.csv("Databases/DrugBank/drug_external_identifiers.csv", header = TRUE)
map1 <- unique(DrugBank_Drugs_idMap[DrugBank_Drugs_idMap$identifier %in% FimmDrugComb_drugs,c("identifier", "parent_key")])
colnames(map1) <- c("drug_names", "DrugBank_drug_id")
row.names(map1) <- NULL
FimmDrugComb_drugs <- FimmDrugComb_drugs[!FimmDrugComb_drugs %in% map1$drug_names]





# Map using metaboliteIDmapping R package
library(metaboliteIDmapping)
tmp <- metabolitesMapping[grep(pattern = "^DB[0-9]+$", x = metabolitesMapping$Drugbank), ]
map2 <- as.data.frame(unique(tmp[tmp$Name %in% FimmDrugComb_drugs, c("Name", "Drugbank")]))
colnames(map2) <- c("drug_names", "DrugBank_drug_id")
row.names(map2) <- NULL
FimmDrugComb_drugs <- FimmDrugComb_drugs[!FimmDrugComb_drugs %in% map2$drug_names]


tmp <- metabolitesMapping[grep(pattern = "^DB[0-9]+$", x = metabolitesMapping$Drugbank), ]
map3 <- as.data.frame(unique(tmp[tmp$CAS %in% FimmDrugComb_drugs, c("CAS", "Drugbank")]))
colnames(map3) <- c("drug_names", "DrugBank_drug_id")
row.names(map3) <- NULL
FimmDrugComb_drugs <- FimmDrugComb_drugs[!FimmDrugComb_drugs %in% map3$drug_names]
rm(tmp)





# Map using WebChem
library(webchem)

# names_2_cas <- cir_query(identifier = FimmDrugComb_drugs[1:20], representation = "cas", match = "all", verbose = TRUE)
# names_2_cas <- func_list_2_df(names_2_cas)
# names_2_cas <- names_2_cas[!names_2_cas$V2 == "NA", ]
# colnames(names_2_cas) <- c("drug_name", "cas")
# saveRDS(names_2_cas, "tmp_dir/names_2_cas.rds")


if(file.exists("tmp_dir/names_2_cas.rds")){
  names_2_cas <- readRDS("tmp_dir/names_2_cas.rds")
}else{names_2_cas <- data.frame()}
FimmDrugComb_drugs <- FimmDrugComb_drugs[!FimmDrugComb_drugs %in% names_2_cas$drug_name]
names_2_cas_new <- cir_query(identifier = FimmDrugComb_drugs, representation = "cas", match = "all", verbose = TRUE)
names_2_cas_new <- func_list_2_df(names_2_cas_new)
names_2_cas_new <- names_2_cas_new[!names_2_cas_new$V2 == "NA", ]
colnames(names_2_cas_new) <- c("drug_name", "cas")
names_2_cas <- rbind(names_2_cas, names_2_cas_new)
saveRDS(names_2_cas, "tmp_dir/names_2_cas.rds")


#Map the retrieved CAS to DrugBank drug ID
DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$cas_number %in% names_2_cas$cas, c("cas_number", "primary_key")]
map4 <- merge(names_2_cas, DrugBank_Drugs, by.x = "cas", by.y = "cas_number")
map4 <- unique(map4[, c("drug_name", "primary_key")])
colnames(map4) <- c("drug_names", "DrugBank_drug_id")
row.names(map4) <- NULL


FimmDrugComb_drugs_idMap <- unique(rbind(map1, map2, map3, map4))
FimmDrugComb_drugs_idMap <- rbind(FimmDrugComb_drugs_idMap_old, FimmDrugComb_drugs_idMap)
row.names(FimmDrugComb_drugs_idMap) <- NULL

if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(FimmDrugComb_drugs_idMap, "InputFiles/ReferenceList/FimmDrugComb_drugs_idMap.rds")


print(warnings())