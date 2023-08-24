set.seed(5081)

# TO ADD
library(dplyr)
library(tidyr)
library(ggplot2)

# Drug combinations from FIMM DrugComb (https://drugcomb.fimm.fi/)



# Load libraries
library(httr)
library(jsonlite)
httr::set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
source("Scripts/Functions/Functions_data_manipulation.R")



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

# TO ADD
# FimmDrugComb_drugCombCat$synergy_loewe[which(FimmDrugComb_drugCombCat$synergy_loewe == "\\N")] <- 0
# FimmDrugComb_drugCombCat$synergy_loewe <- as.numeric(FimmDrugComb_drugCombCat$synergy_loewe) # because of the "\\N" this synergistic score is uploaded as a character

# Get cell line information
FimmDrugComb_cellLine <- GET("https://api.drugcomb.org/cell_lines")
FimmDrugComb_cellLine <- fromJSON(rawToChar(FimmDrugComb_cellLine$content))

# TO ADD
# FimmDrugComb_cellLine$tissue <- sapply(FimmDrugComb_cellLine$ccle_name, function(x) substr(x, gregexpr('\\_', x)[[1]]+1, nchar(x)))

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

# TO ADD
FimmDrugComb_drugCombCat <- merge(FimmDrugComb_drugCombCat, FimmDrugComb_cellLine[,c('name','tissue')], 
                                  by.x="cell_line_name", by.y="name", all.x=TRUE)

# Identify synergistic, antagonist and neutral drugs

# When you aggregate synergy scores across multiple cell lines of the same tumor tissue, 
# you're essentially estimating the "typical" effect of a drug pair on that tumor type. 
# This can be useful for drawing broader conclusions, especially if there is substantial variability between cell lines of the same tissue type. 
# The mean of these average scores then becomes a representative metric for that drug pair's effectiveness against a particular tumor type.

# If you compute the expected mean synergistic score for a given cancer type (by averaging all the average scores for that cancer type), 
# you're essentially setting a "baseline" or "benchmark" of synergy for that cancer type.

# Given this setup, you can apply the same reasoning to classify drug interactions at this aggregated level:

# Synergy: Observed average effect > Expected mean effect for the cancer type.
# Antagonism: Observed average effect < Expected mean effect for the cancer type.
# Additivity: Observed average effect = Expected mean effect for the cancer type.

# Note
# If your mean synergy score is positive and an individual score is just below this mean, 
# this might indicate that the interaction is less synergistic than typical for that tumor type, 
# but not necessarily antagonistic. It might still be synergistic relative to a global baseline (like 0 or 1, depending on the scale of your scores), 
# but less so than the average for that tumor type.

# TO ADD <---------------------------------------------------- big section>
# compute the average of synergistic scores: each drug combinatation evaluated across the different cancer cell lines matching a cancer type/tissue
dcc_agg_all <- FimmDrugComb_drugCombCat %>%
  group_by(drug_row, drug_col, tissue)  %>%
  summarise(AvgSynSM = mean(S_mean, na.rm=TRUE),
            AvgSynZIP = mean(synergy_zip, na.rm=TRUE),
            AvgSynLoewe = mean(synergy_loewe, na.rm=TRUE),
            AvgSynHSA = mean(synergy_hsa, na.rm=TRUE),
            AvgSynBliss = mean(synergy_bliss, na.rm=TRUE), .groups="drop")

# plot the distributions of average-scores
df_long <- dcc_agg_all %>%
  pivot_longer(
    cols = starts_with("AvgSyn"), 
    names_to = "SynergyType", 
    values_to = "Value"
  )

ggplot(df_long, aes(x=tissue, y=Value, fill=tissue)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ SynergyType, scales = "free_y") +
  labs(
    title = "Average Synergy Scores Across Cancer Types",
    x = "Cancer Type",
    y = "Average Synergy Value"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(-150, 150)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# compute general statistics for the distributions of average-scores
cancer_stats <- dcc_agg_all %>%
  group_by(tissue) %>%
  summarise(meanSynS = mean(AvgSynSM, na.rm=TRUE),
            sdSynS = sd(AvgSynSM, na.rm=TRUE),
            meanSynZIP = mean(AvgSynZIP, na.rm=TRUE),
            sdSynZIP = sd(AvgSynZIP, na.rm=TRUE),
            meanSynLoewe = mean(AvgSynLoewe, na.rm=TRUE),
            sdSynLoewe = sd(AvgSynLoewe, na.rm=TRUE),
            meanSynHSA = mean(AvgSynHSA, na.rm=TRUE),
            sdSynHSA = sd(AvgSynHSA, na.rm=TRUE),
            meanSynBliss = mean(AvgSynBliss, na.rm=TRUE),
            sdSynBliss = sd(AvgSynBliss, na.rm=TRUE))

# Identify synergistic, antagonist and neutral drugs
df_with_stats <- left_join(dcc_agg_all, cancer_stats, by="tissue")
df_with_stats$cS = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cZIP = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cLoewe = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cHSA = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cBliss = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$totSyn = rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
k <- 1
for(i in 1:nrow(df_with_stats)) {
  # S
  if(df_with_stats$AvgSynSM[i] > df_with_stats$meanSynS[i] + k * df_with_stats$sdSynS[i] & df_with_stats$AvgSynSM[i] > 0) df_with_stats$cS[i] = 1
  if(df_with_stats$AvgSynSM[i] < df_with_stats$meanSynS[i] + k * df_with_stats$sdSynS[i] & df_with_stats$AvgSynSM[i] < 0) df_with_stats$cS[i] = -1
  # ZIP
  if(df_with_stats$AvgSynZIP[i] > df_with_stats$meanSynZIP[i] + k * df_with_stats$sdSynZIP[i] & df_with_stats$AvgSynZIP[i] > 0) df_with_stats$cZIP[i] = 1
  if(df_with_stats$AvgSynZIP[i] < df_with_stats$meanSynZIP[i] + k * df_with_stats$sdSynZIP[i] & df_with_stats$AvgSynZIP[i] < 0) df_with_stats$cZIP[i] = -1
  # LOEWE
  if(df_with_stats$AvgSynLoewe[i] > df_with_stats$meanSynLoewe[i] + k * df_with_stats$sdSynLoewe[i] & df_with_stats$AvgSynLoewe[i] > 0) df_with_stats$cLoewe[i] = 1
  if(df_with_stats$AvgSynLoewe[i] < df_with_stats$meanSynLoewe[i] + k * df_with_stats$sdSynLoewe[i] & df_with_stats$AvgSynLoewe[i] < 0) df_with_stats$cLoewe[i] = -1
  # HSA
  if(df_with_stats$AvgSynHSA[i] > df_with_stats$meanSynHSA[i] + k * df_with_stats$sdSynHSA[i] & df_with_stats$AvgSynHSA[i] > 0) df_with_stats$cHSA[i] = 1
  if(df_with_stats$AvgSynHSA[i] < df_with_stats$meanSynHSA[i] + k * df_with_stats$sdSynHSA[i] & df_with_stats$AvgSynHSA[i] < 0) df_with_stats$cHSA[i] = -1
  # BLISS
  if(df_with_stats$AvgSynBliss[i] > df_with_stats$meanSynBliss[i] + k * df_with_stats$sdSynBliss[i] & df_with_stats$AvgSynBliss[i] > 0) df_with_stats$cBliss[i] = 1
  if(df_with_stats$AvgSynBliss[i] < df_with_stats$meanSynBliss[i] + k * df_with_stats$sdSynBliss[i] & df_with_stats$AvgSynBliss[i] < 0) df_with_stats$cZIP[i] = -1
  # consensus
  df_with_stats$totSyn[i] = df_with_stats$cS[i] + df_with_stats$cZIP[i] + df_with_stats$cLoewe[i] + df_with_stats$cHSA[i] + df_with_stats$cBliss[i]
}

FimmDrugComb_drugCombCat2 <- df_with_stats

# <---------------------------------------------------- end big section>


# Map drugs to DrugBank drug ID
FimmDrugComb_drugs <- GET("https://api.drugcomb.org/drugs")
FimmDrugComb_drugs <- fromJSON(rawToChar(FimmDrugComb_drugs$content))


FimmDrugComb_drugCombCat$Drug1_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_row, FimmDrugComb_drugs$dname)]
FimmDrugComb_drugCombCat$Drug2_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_col, FimmDrugComb_drugs$dname)]


FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug1_DrugBank_drug_id != "NA"),]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug2_DrugBank_drug_id != "NA"),]

# TO ADD
FimmDrugComb_drugCombCat2$drug_class <- rep('neutral', nrow(FimmDrugComb_drugCombCat2))
FimmDrugComb_drugCombCat2$drug_class[which(FimmDrugComb_drugCombCat2$totSyn > 3)] = "synergy"
FimmDrugComb_drugCombCat2$drug_class[which(FimmDrugComb_drugCombCat2$totSyn < -3)] = "antagonism"


# Assign categories based on the scores
# Positive scores infer synergism and negative scores infer antagonism
# For our study selecting only drug pairs voted synegism/antagonism in all four score types
#FimmDrugComb_drugCombCat$synergy_zip_class <- ifelse(FimmDrugComb_drugCombCat$synergy_zip > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_zip < 0, -1, 0))
#FimmDrugComb_drugCombCat$synergy_loewe_class <- ifelse(FimmDrugComb_drugCombCat$synergy_loewe > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_loewe < 0, -1, 0))
#FimmDrugComb_drugCombCat$synergy_hsa_class <- ifelse(FimmDrugComb_drugCombCat$synergy_hsa > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_hsa < 0, -1, 0))
#FimmDrugComb_drugCombCat$synergy_bliss_class <- ifelse(FimmDrugComb_drugCombCat$synergy_bliss > 0, 1, ifelse(FimmDrugComb_drugCombCat$synergy_bliss < 0, -1, 0))

#FimmDrugComb_drugCombCat$synergy_class_sum <- rowSums(FimmDrugComb_drugCombCat[, 
                                                                               c("synergy_zip_class", "synergy_loewe_class", 
                                                                                 "synergy_hsa_class", "synergy_bliss_class")], na.rm = TRUE)

#FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[abs(FimmDrugComb_drugCombCat$synergy_class_sum) == 4, ]
#FimmDrugComb_drugCombCat$drug_class <- NULL
#FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$synergy_class_sum == 4, "drug_class"] <- "synergy"
#FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$synergy_class_sum == -4, "drug_class"] <- "antagonism"


#FimmDrugComb_drugCombCat <- split(FimmDrugComb_drugCombCat, f = FimmDrugComb_drugCombCat$tissue_name)
#FimmDrugComb_drugCombCat <- lapply(FimmDrugComb_drugCombCat, function(x){split(x, f = x$drug_class)})

if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(FimmDrugComb_drugCombCat, "InputFiles/ReferenceList/FimmDrugComb_drugCombinations.rds")



# Count number of combinations
tmp <- lapply(FimmDrugComb_drugCombCat, function(x)lapply(x, function(x){nrow(unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")]))}))
tmp <- do.call(rbind, tmp)
write.csv(tmp, "Databases/FimmDrugComb/FimmDrugComb_drug_categories_by_tissue.csv")


# TO ADD


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
if(!dir.exists("OutputFiles/Tables/"))dir.create("OutputFiles/Tables/", recursive = TRUE)
write.csv(tmp, "OutputFiles/Tables/DrugComb_cellLines_training.csv")



print(warnings())
