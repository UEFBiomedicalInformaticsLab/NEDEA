# Drug combinations from FIMM DrugComb (https://drugcomb.fimm.fi/)





# Load libraries
library(httr)
library(jsonlite)
httr::set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
source("Scripts/Functions/Functions_ID_Conversion.R")





# Download drug combinations from FIMM DrugComb 
if(!dir.exists("Databases/FimmDrugComb/"))dir.create("Databases/FimmDrugComb/")
if(!file.exists("Databases/FimmDrugComb/summary_v_1_5.csv")){
  download.file(url = "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv",
                destfile = "Databases/FimmDrugComb/summary_v_1_5.csv", method = "wget")
}


# Read the DrugComb data
FimmDrugComb_drugCombCat <- read.csv("Databases/FimmDrugComb/summary_v_1_5.csv", header = TRUE)
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_row != "NULL", ]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_col != "NULL", ]


# Get cell line information
FimmDrugComb_cellLine <- GET("https://api.drugcomb.org/cell_lines")
FimmDrugComb_cellLine <- fromJSON(rawToChar(FimmDrugComb_cellLine$content))


# Download disease IDs from NCI Thesaurus (NCIt) 
if(!dir.exists("Databases/NCIt/"))dir.create("Databases/NCIt/")
if(!file.exists("Databases/NCIt/Thesaurus.txt")){
  download.file(url = "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/22.12d_Release/Thesaurus_22.12d.FLAT.zip",
                destfile = "Databases/NCIt/Thesaurus_22.12d.FLAT.zip", method = "wget")
  unzip("Databases/NCIt/Thesaurus_22.12d.FLAT.zip", exdir = "Databases/NCIt/", file = "Thesaurus.txt")
}
NCIthesaurus <- read.table("Databases/NCIt/Thesaurus.txt", sep = "\t", header = FALSE, comment.char = "", fill = TRUE, quote = "")
colnames(NCIthesaurus) <- c("code", "concept IRI", "parents", "synonyms", "definition", "display name", "concept status", "semantic type")


# Filter drug combinations tested on cancer related cell lines
NCIthesaurus <- NCIthesaurus[NCIthesaurus$code %in% FimmDrugComb_cellLine$disease_id, ]
NCIthesaurus <- NCIthesaurus[grep("cancer|carcinoma|sarcoma|lymphoma|leukemia|melanoma", NCIthesaurus$synonyms, ignore.case = TRUE), ]

FimmDrugComb_cellLine <- FimmDrugComb_cellLine[FimmDrugComb_cellLine$disease_id %in% NCIthesaurus$code, ]

FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$cell_line_name %in% FimmDrugComb_cellLine$name, ]



# Map drugs to DrugBank drug ID
FimmDrugComb_drugs <- GET("https://api.drugcomb.org/drugs")
FimmDrugComb_drugs <- fromJSON(rawToChar(FimmDrugComb_drugs$content))


FimmDrugComb_drugCombCat$Drug1_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_row, FimmDrugComb_drugs$dname)]
FimmDrugComb_drugCombCat$Drug2_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_col, FimmDrugComb_drugs$dname)]


FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug1_DrugBank_drug_id != "NA"),]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug2_DrugBank_drug_id != "NA"),]



# Assign categories based on the scores
# Positive scores infer synergism and negative scores infer antagonism
# For our study selecting only drug pairs voted synegism/antagonism in all four score types
FimmDrugComb_drugCombCat$synergy_zip_class <- ifelse(FimmDrugComb_drugCombCat$synergy_zip > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_zip < 0, -1, 0))
FimmDrugComb_drugCombCat$synergy_loewe_class <- ifelse(FimmDrugComb_drugCombCat$synergy_loewe > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_loewe < 0, -1, 0))
FimmDrugComb_drugCombCat$synergy_hsa_class <- ifelse(FimmDrugComb_drugCombCat$synergy_hsa > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_hsa < 0, -1, 0))
FimmDrugComb_drugCombCat$synergy_bliss_class <- ifelse(FimmDrugComb_drugCombCat$synergy_bliss > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_bliss < 0, -1, 0))

FimmDrugComb_drugCombCat$synergy_class_sum <- rowSums(FimmDrugComb_drugCombCat[, 
                                                      c("synergy_zip_class", "synergy_loewe_class", 
                                                        "synergy_hsa_class", "synergy_bliss_class")], na.rm = TRUE)

FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[abs(FimmDrugComb_drugCombCat$synergy_class_sum) == 4, ]
FimmDrugComb_drugCombCat$drug_class <- NULL
FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$synergy_class_sum == 4, "drug_class"] <- "synergy"
FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$synergy_class_sum == -4, "drug_class"] <- "antagonism"


FimmDrugComb_drugCombCat <- split(FimmDrugComb_drugCombCat, f = FimmDrugComb_drugCombCat$tissue_name)
FimmDrugComb_drugCombCat <- lapply(FimmDrugComb_drugCombCat, function(x){split(x, f = x$drug_class)})

if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(FimmDrugComb_drugCombCat, "InputFiles/ReferenceList/FimmDrugComb_drugCombinations.rds")


# FimmDrugComb_drugCombCat$lung$synergy$study_name




# Count number of combinations
tmp <- lapply(FimmDrugComb_drugCombCat, function(x)lapply(x, function(x){nrow(unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")]))}))
tmp <- do.call(rbind, tmp)
write.csv(tmp, "Databases/FimmDrugComb/FimmDrugComb_drug_categories_by_tissue.csv")





# LungCancer specific
tmp <- FimmDrugComb_drugCombCat$lung
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
  })
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_LungCancer_drugCombinations.rds")





# BreastCancer specific
tmp <- FimmDrugComb_drugCombCat$breast
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
})
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_BreastCancer_drugCombinations.rds")





# KidneyCancer specific
tmp <- FimmDrugComb_drugCombCat$kidney
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
})
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_KidneyCancer_drugCombinations.rds")





# OvaryCancer specific
tmp <- FimmDrugComb_drugCombCat$ovary
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
})
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_OvaryCancer_drugCombinations.rds")





# ProstateCancer specific
tmp <- FimmDrugComb_drugCombCat$prostate
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
})
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_ProstateCancer_drugCombinations.rds")





# SkinCancer specific
tmp <- FimmDrugComb_drugCombCat$skin
tmp <- lapply(tmp, function(x){
  x <- unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")])
  row.names(x) <- NULL
  x
})
saveRDS(tmp, "InputFiles/ReferenceList/FimmDrugComb_SkinCancer_drugCombinations.rds")



# Print the cell lines considered for each cancer types              
drugCombs <- readRDS("InputFiles/ReferenceList/FimmDrugComb_drugCombinations.rds")
cell_lines <- lapply(drugCombs, function(x){lapply(x, function(y){unique(y$cell_line_name)})})
cell_lines <- unlist(cell_lines, recursive = FALSE)
cell_lines <- cell_lines[grep("synergy", names(cell_lines))]

cell_lines <- cell_lines[grep("breast|kidney|lung|ovary|prostate|skin", ignore.case = TRUE, names(cell_lines))]
names(cell_lines) <- gsub(".synergy$", "", names(cell_lines))
tmp <- do.call(cbind.fill, cell_lines)
colnames(tmp) <- names(cell_lines)
write.csv(tmp, "OutputFiles/Tables/DrugComb_cellLines_training.csv")
              
              
              
print(warnings())