####################################################################################################### DrugCombo FIMM
set.seed(5081)
# Drug combinations from FIMM DrugComb (https://drugcomb.fimm.fi/)

# Load libraries
library(httr)
library(jsonlite)
httr::set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
source("Scripts/Functions/Functions_data_manipulation.R")

setwd("/Users/vittfo/Documents/DRUGCOMBO/")

# Download drug combinations from FIMM DrugComb 
if(!dir.exists("Databases/FimmDrugComb/"))dir.create("Databases/FimmDrugComb/")
if(!file.exists("Databases/FimmDrugComb/summary_v_1_5.csv")){
  download.file(url = "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv",
                destfile = "Databases/FimmDrugComb/summary_v_1_5.csv", method = "wget")
}

# Read the DrugComb data
FimmDrugComb_drugCombCat <- read.csv("/Users/vittfo/Downloads/summary_v_1_5.csv", header = TRUE)
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_row != "NULL", ]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_col != "NULL", ]
FimmDrugComb_drugCombCat$synergy_loewe[which(FimmDrugComb_drugCombCat$synergy_loewe == "\\N")] = 0
FimmDrugComb_drugCombCat$synergy_loewe = as.numeric(FimmDrugComb_drugCombCat$synergy_loewe)

# Get cell line information
FimmDrugComb_cellLine <- GET("https://api.drugcomb.org/cell_lines")
FimmDrugComb_cellLine <- fromJSON(rawToChar(FimmDrugComb_cellLine$content))
FimmDrugComb_cellLine$tissue = sapply(FimmDrugComb_cellLine$ccle_name, function(x) substr(x, gregexpr('\\_', x)[[1]]+1, nchar(x)))


# Download disease IDs from NCI Thesaurus (NCIt) 
if(!dir.exists("Databases/NCIt/"))dir.create("Databases/NCIt/")
if(!file.exists("Databases/NCIt/Thesaurus.txt")){
  download.file(url = "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/22.12d_Release/Thesaurus_22.12d.FLAT.zip",
                destfile = "Databases/NCIt/Thesaurus_22.12d.FLAT.zip", method = "wget")
  unzip("Databases/NCIt/Thesaurus_22.12d.FLAT.zip", exdir = "Databases/NCIt/", file = "Thesaurus.txt")
}
NCIthesaurus <- read.table("/Users/vittfo/Downloads/Thesaurus.txt", sep = "\t", header = FALSE, comment.char = "", fill = TRUE, quote = "")
colnames(NCIthesaurus) <- c("code", "concept IRI", "parents", "synonyms", "definition", "display name", "concept status", "semantic type")

# Filter drug combinations tested on cancer related cell lines
NCIthesaurus <- NCIthesaurus[NCIthesaurus$code %in% FimmDrugComb_cellLine$disease_id, ]
NCIthesaurus <- NCIthesaurus[grep("cancer|carcinoma|sarcoma|lymphoma|leukemia|melanoma", NCIthesaurus$synonyms, ignore.case = TRUE), ]

FimmDrugComb_cellLine <- FimmDrugComb_cellLine[FimmDrugComb_cellLine$disease_id %in% NCIthesaurus$code, ]

FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$cell_line_name %in% FimmDrugComb_cellLine$name, ]

FimmDrugComb_drugCombCat <- merge(FimmDrugComb_drugCombCat, FimmDrugComb_cellLine[,c('name','tissue')], 
                                  by.x="cell_line_name", by.y="name", all.x=TRUE)

# * ZIP (Zero Interaction Potency):
# The ZIP model calculates the expected effect of two drugs if there were no interaction between them. 
# The observed combined effect is then compared to this expected effect.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * Loewe Additivity:
# The Loewe model uses a combination index (CI) derived from the expected dose of each drug 
# to achieve a particular effect if used alone versus when used in combination.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * Bliss Independence:
# Bliss independence calculates the expected effect assuming that two drugs act independently.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * HSA (Highest Single Agent):
# HSA compares the effect of a drug combination to the effect of the most effective drug in the combination used alone.
# 
# Synergy: Combined effect > Effect of the most effective single agent.
# Antagonism: Combined effect < Effect of the most effective single agent.
# Additivity: Combined effect = Effect of the most effective single agent.
#
# * FIMM Scores: The FIMM group introduced some scores to quantify drug interactions:
#   - S-score (Synergy score): This score reflects the difference between the observed and expected responses. 
#     A score near zero indicates additivity. A positive score indicates synergy, 
#     and a negative score may indicate antagonism, although the interpretation of negative scores can sometimes be ambiguous.

########################################################################################
####  1. Assess drug synergy scores across cancer cell lines.                       ####
####     It generates synergistic scores for each drug pair and cancer type/tissue. ####
########################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

dcc_agg_all <- FimmDrugComb_drugCombCat %>%
  group_by(drug_row, drug_col, tissue)  %>%
  summarise(AvgSynSM = mean(S_mean, na.rm=TRUE),
            AvgSynZIP = mean(synergy_zip, na.rm=TRUE),
            AvgSynLoewe = mean(synergy_loewe, na.rm=TRUE),
            AvgSynHSA = mean(synergy_hsa, na.rm=TRUE),
            AvgSynBliss = mean(synergy_bliss, na.rm=TRUE), .groups="drop")

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

# Identify synergistic, antagonist and neutral drugs

# When you aggregate synergy scores across multiple cell lines of the same tumor tissue or cancer type, 
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

df_with_stats <- left_join(dcc_agg_all, cancer_stats, by="tissue")
#
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


# Map drugs to DrugBank drug ID
FimmDrugComb_drugs <- GET("https://api.drugcomb.org/drugs")
FimmDrugComb_drugs <- fromJSON(rawToChar(FimmDrugComb_drugs$content))

FimmDrugComb_drugCombCat2 <- df_with_stats

FimmDrugComb_drugCombCat2$Drug1_DBank_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat2$drug_row, FimmDrugComb_drugs$dname)]
FimmDrugComb_drugCombCat2$Drug2_DBank_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat2$drug_col, FimmDrugComb_drugs$dname)]
FimmDrugComb_drugCombCat2 <- FimmDrugComb_drugCombCat2[(FimmDrugComb_drugCombCat2$Drug1_DBank_id != "NA" & 
                                                          FimmDrugComb_drugCombCat2$Drug2_DBank_id != "NA"),]

FimmDrugComb_drugCombCat2$drug_class <- rep('neutral', nrow(FimmDrugComb_drugCombCat2))
FimmDrugComb_drugCombCat2$drug_class[which(FimmDrugComb_drugCombCat2$totSyn >= 3)] = "synergy"
FimmDrugComb_drugCombCat2$drug_class[which(FimmDrugComb_drugCombCat2$totSyn <= -3)] = "antagonism"


if(!dir.exists("InputFiles/ReferenceList/"))dir.create("InputFiles/ReferenceList/")
saveRDS(FimmDrugComb_drugCombCat2, "InputFiles/ReferenceList/FimmDrugComb_drugCombinations.rds")


# Count number of combinations
#tmp <- lapply(FimmDrugComb_drugCombCat, function(x)lapply(x, function(x){nrow(unique(x[, c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")]))}))
#tmp <- do.call(rbind, tmp)
#write.csv(tmp, "Databases/FimmDrugComb/FimmDrugComb_drug_categories_by_tissue.csv")


########################################################################################
####  2. Use drugbank info on drug combinations to report annoations regarding:     ####
####     drug adverse reaction and pharmacokinetics information.                    ####
########################################################################################

ddi = read.csv("/Users/vittfo/Downloads/drug_drug_interactions.csv")
colnames(ddi)[c(1,4)] <- c("Drug1_DBank_id", "Drug2_DBank_id")

find_top_k_frequent_sequences <- function(strings_list, n = 2, k = 10) {
  require(stringi)  # For string operations
  
  # Generate n-grams
  generate_ngrams <- function(text, n) {
    words <- unlist(strsplit(text, split = " "))
    l <- length(words)
    ngrams <- sapply(1:(l-n+1), function(i) paste(words[i:(i+n-1)], collapse = " "))
    return(ngrams)
  }
  
  all_ngrams <- unlist(lapply(strings_list, function(x) generate_ngrams(x, n)))
  
  # Count and sort n-grams
  ngram_freq <- table(all_ngrams)
  sorted_ngram_freq <- sort(ngram_freq, decreasing = TRUE)
  
  # Return the top k frequent n-grams and their counts
  top_k <- min(k, length(sorted_ngram_freq))
  list(ngrams = names(sorted_ngram_freq)[1:top_k], counts = as.integer(sorted_ngram_freq[1:top_k]))
}
find_top_k_sequences_with_keyword <- function(strings_list, keyword, n = 2, k = 10) {
  require(stringi)
  
  # Generate n-grams
  generate_ngrams <- function(text, n) {
    words <- unlist(strsplit(text, split = " "))
    l <- length(words)
    ngrams <- sapply(1:(l-n+1), function(i) paste(words[i:(i+n-1)], collapse = " "))
    return(ngrams)
  }
  
  # Filter sequences with the keyword
  filter_keyword <- function(ngram_list, keyword) {
    ngram_list[sapply(ngram_list, function(ngram) grepl(keyword, ngram))]
  }
  
  all_ngrams <- unlist(lapply(strings_list, function(x) generate_ngrams(x, n)))
  
  keyword_ngrams <- filter_keyword(all_ngrams, keyword)
  
  # If no sequences with keyword, return empty list
  if (length(keyword_ngrams) == 0) {
    return(list(ngrams = character(0), counts = integer(0)))
  }
  
  # Count and sort n-grams
  ngram_freq <- table(keyword_ngrams)
  sorted_ngram_freq <- sort(ngram_freq, decreasing = TRUE)
  
  # Return the top k frequent n-grams and their counts
  top_k <- min(k, length(sorted_ngram_freq))
  list(ngrams = names(sorted_ngram_freq)[1:top_k], counts = as.integer(sorted_ngram_freq[1:top_k]))
}

top_freq = find_top_k_frequent_sequences(ddi$description, n = 4, 100)
risk_sev_res = find_top_k_sequences_with_keyword(ddi$description, keyword = 'The risk or severity', n = 6, k = 50)

test_eff_inc = grepl("The therapeutic efficacy of", ddi$description) & grepl("can be increased when used in combination", ddi$description)
test_eff_dec = grepl("The therapeutic efficacy of", ddi$description) & grepl("can be decreased when used in combination", ddi$description)
test_anticinc = grepl("may increase the anticoagulant activities", ddi$description)
test_anticdec = grepl("may decrease the anticoagulant activities", ddi$description)
test_antihinc = grepl("may increase the antihypertensive activities", ddi$description)
test_antihdec = grepl("may decrease the antihypertensive activities", ddi$description)
test_nsd = grepl("(CNS depressant)", ddi$description)
test_scinc = grepl("The serum concentration of", ddi$description) & grepl("can be increased", ddi$description)
test_scdec = grepl("The serum concentration of", ddi$description) & grepl("can be decreased", ddi$description)
test_meinc = grepl("The metabolism of", ddi$description) & grepl("can be increased", ddi$description)
test_medec = grepl("The metabolism of", ddi$description) & grepl("can be decreased", ddi$description)
test_absde = grepl("a decrease in the absorption ", ddi$description)
test_all_risk = grepl("risk or severity", ddi$description)
#
list_tox_test = list()
for(i in 1:length(risk_sev_res$ngrams)){
  list_tox_test[[length(list_tox_test) + 1]] = grepl(risk_sev_res$ngrams[i], ddi$description)
}
names(list_tox_test) = risk_sev_res$ngrams
#
sum(test_eff_inc)
sum(test_eff_dec)
sum(test_anticinc)
sum(test_anticdec)
sum(test_antihinc)
sum(test_antihdec)
sum(test_nsd)
sum(test_scinc)
sum(test_scdec)
sum(test_meinc)
sum(test_medec)
sum(test_absde)
sum(test_all_risk)
#
ddi$effVsNeff = rep('Unk', nrow(ddi))
ddi$effVsNeff[test_eff_inc] = "IncEff"
ddi$effVsNeff[test_eff_dec] = "DecEff"
#
ddi$phamk = rep('Unk', nrow(ddi))
ddi$phamk[test_meinc | test_scdec | test_absde] = "PDecEff"
ddi$phamk[test_medec | test_scinc] = "PIncEff"
#
#### For all cancers
# The risk or severity of adverse effect (not specific)
# The risk or severity of QTc
# The risk or severity of hypertension
# The risk or severity of Tachycardia
# The risk or severity of hyperglycemia
# The risk or severity of hypotension
# The risk or severity of renal
# The risk or severity of infection
# The risk or severity of CNS
# The risk or severity of edema
# The risk or severity of neutropenia
# The risk or severity of myelosuppression
# The risk or severity of cardiotoxicity
# The risk or severity of liver 
# The risk or severity of hyperthermia

set_adv_ids = c(1, 2, 4, 6, 10, 11, 12, 15, 16, 20, 28, 33, 35, 39, 48, 49) 

#### Breast cancer adverse effects
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of peripheral neuropathy 
# The risk or severity of electrolyte imbalance
names(list_tox_test)[c(set_adv_ids, 25, 31, 47)]
sapply(list_tox_test[c(set_adv_ids, 25, 31, 47)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 25, 31, 47)]))
#
ddi$breastADV = rep('Unk', nrow(ddi))
ddi$breastADV[Reduce('|',list_tox_test[c(set_adv_ids, 25, 47)])] = "Adv"

#### Lung cancer adverse effects
# The risk or severity of bleeding
# The risk or severity of Thrombosis
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of respiratory
# The risk or severity of hyponatremia
# The risk or severity of electrolyte imbalance
names(list_tox_test)[c(set_adv_ids, 3, 19, 25, 31, 34, 46)]
sapply(list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)]))
#
ddi$lungADV = rep('Unk', nrow(ddi))
ddi$lungADV[Reduce('|',list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)])] = "Adv"

#### Prostate cancer adverse effects
# The risk or severity of fluid retention
# The risk or severity of urinary retention
# The risk or severity of hypercalcemia
names(list_tox_test)[c(set_adv_ids, 31, 38, 42, 43)]
sapply(list_tox_test[c(set_adv_ids, 31, 38, 42, 43)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 31, 38, 42, 43)]))
#
ddi$prostateADV = rep('Unk', nrow(ddi))
ddi$prostateADV[Reduce('|',list_tox_test[c(set_adv_ids, 31, 38, 42, 43)])] = "Adv"                
                   
#### Ovarian cancer adverse effects  
# The risk or severity of bleeding
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of urinary retention
# The risk or severity of hypercalcemia
# The risk or severity of respiratory
names(list_tox_test)[c(set_adv_ids, 3, 25, 42, 43, 46)]
sapply(list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)]))
#
ddi$ovarianADV = rep('Unk', nrow(ddi))
ddi$ovarianADV[Reduce('|',list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)])] = "Adv"                

#### Kidney/renal cancer adverse effects  
# The risk or severity of bleeding
# The risk or severity of nephrotoxicity
# The risk or severity of electrolyte
# The risk or severity of fluid retention
# The risk or severity of hyperkalemia
# The risk or severity of hyponatremia
# The risk or severity of urinary
# The risk or severity of hypercalcemia
names(list_tox_test)[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)]
sapply(list_tox_test[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)]))
#
ddi$kidneyADV = rep('Unk', nrow(ddi))
ddi$kidneyADV[Reduce('|',list_tox_test[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)])] = "Adv"   
#
#### Skin cancer adverse effects  
# The risk or severity of respiratory: if skin cancer metastasizes to the lungs, some treatments might have respiratory side effects.
names(list_tox_test)[c(set_adv_ids, 46)]
sapply(list_tox_test[c(set_adv_ids, 46)], sum)
sum(Reduce('|',list_tox_test[c(set_adv_ids, 46)]))
#
ddi$skinADV = rep('Unk', nrow(ddi))
ddi$skinADV[Reduce('|',list_tox_test[c(set_adv_ids, 46)])] = "Adv"   
#
#

#######################################################
####  3. Merging DrugCombo-FIMM and DrugBank info. ####
#######################################################

# Merging the dataframes based on matching drug ID combinations
fimm_db_all = merge(FimmDrugComb_drugCombCat2, ddi, 
                    by.x = c("Drug1_DBank_id", "Drug2_DBank_id"), 
                    by.y = c("Drug1_DBank_id", "Drug2_DBank_id"), 
                    all.x = TRUE)
# 
breast_tissue = which(fimm_db_all$tissue == "BREAST")
lung_tissue = which(fimm_db_all$tissue == "LUNG")
kidney_tissue = which(fimm_db_all$tissue == "KIDNEY")
prostate_tissue = which(fimm_db_all$tissue == "PROSTATE")
ovarian_tissue = which(fimm_db_all$tissue == "OVARY")
skin_tissue = which(fimm_db_all$tissue == "SKIN")
#
length(which(complete.cases(fimm_db_all)))
table(fimm_db_all$totSyn, fimm_db_all$effVsNeff)
table(fimm_db_all$totSyn, fimm_db_all$phamk)
print('Association between total synergistic score and ADV info')
print('<breast>')
table(fimm_db_all$totSyn[breast_tissue], fimm_db_all$phamk[breast_tissue])
table(fimm_db_all$totSyn[breast_tissue], fimm_db_all$breastADV[breast_tissue])
print('<lung>')
table(fimm_db_all$totSyn[lung_tissue], fimm_db_all$phamk[lung_tissue])
table(fimm_db_all$totSyn[lung_tissue], fimm_db_all$lungADV[lung_tissue])
print('<kidney>')
table(fimm_db_all$totSyn[kidney_tissue], fimm_db_all$phamk[kidney_tissue])
table(fimm_db_all$totSyn[kidney_tissue], fimm_db_all$kidneyADV[kidney_tissue])
print('<prostate>')
table(fimm_db_all$totSyn[prostate_tissue], fimm_db_all$phamk[prostate_tissue])
table(fimm_db_all$totSyn[prostate_tissue], fimm_db_all$prostateADV[prostate_tissue])
print('<ovarian>')
table(fimm_db_all$totSyn[ovarian_tissue], fimm_db_all$phamk[ovarian_tissue])
table(fimm_db_all$totSyn[ovarian_tissue], fimm_db_all$ovarianADV[ovarian_tissue])
print('<skin>')
table(fimm_db_all$totSyn[skin_tissue], fimm_db_all$phamk[skin_tissue])
table(fimm_db_all$totSyn[skin_tissue], fimm_db_all$skinADV[skin_tissue])
#

# BreastCancer specific
breast_dp = fimm_db_all[breast_tissue,]
# LungCancer specific
lung_dp = fimm_db_all[lung_tissue,]
# ProstateCancer specific
prostate_dp = fimm_db_all[prostate_tissue,]
# KidneyCancer specific
kidney_dp = fimm_db_all[kidney_tissue,]
# OvarianCancer specific
ovarian_dp = fimm_db_all[ovarian_tissue,]
# SkinCancer specific
skin_dp = fimm_db_all[skin_tissue,]

################################## Load protein-protein interactions from BioGrid
#library(biomaRt)
#library(data.table)

#BioGrid_ppi <- read.delim("BIOGRID-ALL-4.4.224.tab3.txt", row.names = 1, header = T, na.strings = "-")
#BioGrid_ppi <- BioGrid_ppi[,-c(19:36)]
#BioGrid_ppi <- BioGrid_ppi[which(BioGrid_ppi$Organism.ID.Interactor.A == 9606 & BioGrid_ppi$Organism.ID.Interactor.B == 9606),]
#
#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#
#entrez_ids <- unique(c(BioGrid_ppi$Entrez.Gene.Interactor.A,
#                       BioGrid_ppi$Entrez.Gene.Interactor.B))
#converted <- getBM(
#  attributes = c('entrezgene_id', 'ensembl_gene_id'),
#  filters = 'entrezgene_id',
#  values = entrez_ids,
#  mart = mart
#)
#BioGrid_ppi = as.data.table(BioGrid_ppi)
#converted = as.data.table(converted)
#etkey(converted, entrezgene_id)
#setkey(BioGrid_ppi, Entrez.Gene.Interactor.A)
# Join on entrez_id1 and add the corresponding Ensembl ID
#BioGrid_ppi <- BioGrid_ppi[converted, ensembl_gene_id1 := ensembl_gene_id]
#setkey(BioGrid_ppi, Entrez.Gene.Interactor.B)
# Join on entrez_id2 and add the corresponding Ensembl ID
#BioGrid_ppi <- BioGrid_ppi[converted, ensembl_gene_id2 := ensembl_gene_id]

# Convert to graph object
#Biogrid_ppi_Net <- graph_from_data_frame(na.exclude(BioGrid_ppi[which(BioGrid_ppi$Experimental.System.Type == "physical"), 
#                                                                c("ensembl_gene_id1", "ensembl_gene_id2")], directed = FALSE))
#Biogrid_ppi_Net <- as.undirected(Biogrid_ppi_Net, mode = "collapse")
#Biogrid_ppi_Net <- simplify(Biogrid_ppi_Net, remove.loops = TRUE, remove.multiple	= TRUE) # remove loops and multi-edges
#print(paste("Network size (vertices, edges):", vcount(Biogrid_ppi_Net), ecount(Biogrid_ppi_Net)))

####################################################
####  4. Build the PPI network with StringDB    ####
####################################################

String_ppi <- read.table("9606.protein.links.detailed.v11.5.txt.gz", header = TRUE)

# Filter all interactions with scores greater than 0 and sourced from database/experimets
String_ppi <- String_ppi[(String_ppi$database > 0 | String_ppi$experimental > 0), ] #

String_proteins <- sort(unique(c(String_ppi$protein1, String_ppi$protein2)))

# Create mapping from Ensembl Peptide ID ID to Ensembl Gene ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensemblPeptideID_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                                          mart = mart, 
                                          filters = "ensembl_peptide_id", 
                                          values = gsub("^9606.", "", String_proteins))
ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id <- paste0("9606.", ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)

String_ppi$Node1_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein1, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]
String_ppi$Node2_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein2, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]

# Filter all interactions with scores greater than 0 and sourced from database/experiments
#String_ppi_Net <- String_ppi[(String_ppi$database > 200 | String_ppi$experimental > 200 | String_ppi$fusion > 200), ] #
String_ppi_Net <- String_ppi[(String_ppi$combined_score > 400), ]
dim(String_ppi_Net)

String_ppi_Net <- na.exclude(String_ppi_Net[, c("protein1", "protein2", "Node1_ensembl_gene_id", "Node2_ensembl_gene_id", "combined_score")])
String_ppi_Net$Node1_type <- "gene"
String_ppi_Net$Node2_type <- "gene"
String_ppi_Net$Edge_type <- "undirected"
String_ppi_Net <- String_ppi_Net[, c("Node1_ensembl_gene_id", "Node1_type", "Node2_ensembl_gene_id", "Node2_type", "Edge_type")]

# Convert to graph object
String_ppi_Net <- graph_from_data_frame(String_ppi_Net[, c("Node1_ensembl_gene_id", "Node2_ensembl_gene_id")], directed = FALSE)
String_ppi_Net <- simplify(String_ppi_Net, remove.loops = TRUE, remove.multiple	= TRUE) # remove loops and multi-edges
# Find the connected components
string_comps <- components(String_ppi_Net)
# Identify the largest component
largest_string_comp <- which.max(string_comps$csize)
# Extract the largest component
String_ppi_Net <- induced.subgraph(String_ppi_Net, which(string_comps$membership == largest_string_comp))
print(paste("Network size (vertices, edges):", vcount(String_ppi_Net), ecount(String_ppi_Net)))

save(String_ppi_Net, file="STRINGgraph.RData")
load("STRINGgraph.RData")

################################## Gene sets / Libraries for the enrichment step
lib_files = list.files("Enrich_new")
breast_libs <- readRDS("Enrich_new/Disease2Gene_BreastCancer_lib.rds")
breast_libs <- breast_libs[c(1:2,10:12,14,17:19,25,26,34,35,40:46)] # we should avoid selecting gene sets associated to cancer subtypes

adr_libs <- readRDS("Enrich_new/drugWithdrawal_Adr2Gene_lib.rds")

################################################
####  5. Load drug target interactions and  #### 
####     overlap with the PPI network       ####
################################################

library(dplyr)
library(igraph)

# Read drug target interactions
drug_target_ixn <- readRDS("DrugBank_Drug_Target_Net.rds")

#load("STRINGgraph.RData")
print(dim(drug_target_ixn))

# Filtering the DataFrame based on the genes in the network
drug_target_ixn <- drug_target_ixn %>%
  filter(Node2_ensembl_gene_id %in% V(String_ppi_Net)$name) # String_ppi_Net / Biogrid_ppi_Net
print(dim(drug_target_ixn))

# Create a named list where the names are drug names and the list elements are aggregated genes
aggregated_genes_list <- drug_target_ixn %>%
  group_by(Node1_drugbank_drug_id) %>%
  summarise(aggregated_genes = paste(Node2_ensembl_gene_id, collapse = ","))

# Merge this information into breast_dp
breast_dat <- breast_dp %>%
  left_join(aggregated_genes_list, by = c("Drug1_DBank_id" = "Node1_drugbank_drug_id")) %>%
  dplyr::rename(aggregated_genes1 = aggregated_genes) %>%
  left_join(aggregated_genes_list, by = c("Drug2_DBank_id" = "Node1_drugbank_drug_id")) %>%
  dplyr::rename(aggregated_genes2 = aggregated_genes) %>%
  filter(!is.na(aggregated_genes1), !is.na(aggregated_genes2)) %>%  # Remove rows where either drug has zero targets
  dplyr::rowwise() %>%
  dplyr::mutate(
    aggregated_genes = list(unique(unlist(strsplit(c(aggregated_genes1, aggregated_genes2), ",")))),
    num_unique_genes = length(unlist(aggregated_genes))
  ) %>%
  dplyr::select(-aggregated_genes1, -aggregated_genes2) #%>%

# Convert the list column to a comma-separated string column
breast_dat$aggregated_genes <- sapply(breast_dat$aggregated_genes, function(x) paste(x, collapse = ","))

# TO BE REPEATED FOR EACH CANCER SET - I suggest writing a general function that would generate the dataset based on the 

###############################################
####  6.    Extend known drug targets      ####
###############################################
library(OmnipathR)

inter_PS <- import_omnipath_interactions(resources=c("PhosphoSite","phosphoELM","PhosphoNetworks","PhosphoSite_MIMP",
                                                     "PhosphoSite_ProtMapper","phosphoELM_MIMP","PhosphoPoint",
                                                     "PhosphoSite_KEA","PhosphoSite_noref"))
inter_SIGNOR <- import_omnipath_interactions(resources=c("SIGNOR","SIGNOR_ProtMapper","SIGNOR_CollecTRI"))
inter_NPA <- import_omnipath_interactions(resources=c("NetPath"))
inter_RIs <- import_transcriptional_interactions()

#TRRUST_DoRothEA, CancerDrugsDB, Reactome_LRdb

# Initialize the biomaRt object for human genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert Ensembl IDs to gene symbols
convert_ensembl_to_symbol <- function(ids_str) {
  ids <- unlist(strsplit(ids_str, ","))
  symbols_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                      filters = "ensembl_gene_id", 
                      values = ids, 
                      mart = ensembl)
  return(paste(symbols_df$hgnc_symbol, collapse = ","))
}
convert_ensembl_to_symbol2 <- function(ids_str) {
  ids <- unlist(strsplit(ids_str, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, keys=ids, 
                                     columns="SYMBOL", keytype="ENSEMBL"))
  return(paste(mapping$SYMBOL, collapse = ","))
}

# Function to convert gene symbols to Ensembl IDs
convert_symbol_to_ensembl <- function(ids_str) {
  ids <- unlist(strsplit(ids_str, ","))
  symbols_df <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                      filters = "hgnc_symbol", 
                      values = ids, 
                      mart = ensembl)
  return(paste(symbols_df$ensembl_gene_id, collapse = ","))
}
convert_symbol_to_ensembl2 <- function(ids_str) {
  ids <- unlist(strsplit(ids_str, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, keys=ids, 
                                     columns="ENSEMBL", keytype="SYMBOL"))
  return(paste(mapping_df$ensembl_gene_id, collapse = ","))
}

# Function to extend known drug targets with any 'interaction' database
extend_targets_inter_database <- function(ids_str, inter_data) {
  known_drug_targets <- unlist(strsplit(ids_str, ","))
  # Extend drug targets using the provided interaction database
  extended_targets <- inter_data[inter_data$source_genesymbol %in% known_drug_targets,]
  extended_targets <- extended_targets[which(extended_targets$curation_effort > 1 & extended_targets$is_directed),]
  all_targets <- unique(c(extended_targets$source_genesymbol, extended_targets$target_genesymbol, known_drug_targets))
  return(list("ext_targets" = paste(all_targets, collapse = ","), 
              "ext_tar_cnt" = length(all_targets)))
}

# Function to extend drug targets by using a kegg pathway
extend_targets_kegg <- function(initial_targets, kegg_data) {
  
  extended_targets <- unlist(strsplit(initial_targets, ","))
  
  # Continue until no new targets are found
  new_targets_found <- TRUE
  while (new_targets_found) {
    new_targets_found <- FALSE
    
    # Find interactions where the source gene is in our list of targets
    new_targets <- kegg_data %>%
      filter(genesymbol_source %in% extended_targets) %>%
      pull(genesymbol_target)
    
    # Add unique new targets to our list
    original_length <- length(extended_targets)
    extended_targets <- unique(c(extended_targets, new_targets))
    
    if (length(extended_targets) > original_length) {
      new_targets_found <- TRUE
    }
  }
  return(extended_targets)
}
extend_targets_kegg_list <- function(initial_targets, list_kegg) {
  # Loop to extend targets for each KEGG data frame
  results <- list()
  for (i in seq_along(list_kegg)) {
    kegg_data <- list_kegg[[i]]
    extended_targets <- extend_targets_kegg(initial_targets, kegg_data)
    results[[names(list_kegg)[i]]] <- extended_targets
  }
  return(list("ext_kegg_targets" = paste(Reduce('union',results), collapse = ","), 
              "ext_kegg_tar_cnt" = length(Reduce('union',results))))
  return()
}

# Convert targets for each drug
breast_dat$target_symbols <- sapply(breast_dat$aggregated_genes, convert_ensembl_to_symbol2)

# Using SIGNOR
signor_result_list <- lapply(breast_dat$target_symbols, function(x) extend_targets_inter_database(x, inter_SIGNOR))
breast_dat$ext_sig_targets <- sapply(signor_result_list, function(x) x$ext_targets)
breast_dat$ext_sig_tar_cnt <- sapply(signor_result_list, function(x) x$ext_tar_cnt)

# Using RNs
reg_result_list <- lapply(breast_dat$target_symbols, function(x) extend_targets_inter_database(x, inter_RIs))
breast_dat$ext_reg_targets <- sapply(reg_result_list, function(x) x$ext_targets)
breast_dat$ext_reg_tar_cnt <- sapply(reg_result_list, function(x) x$ext_tar_cnt)

# Using Phospho
pho_result_list <- lapply(breast_dat$target_symbols, function(x) extend_targets_inter_database(x, inter_PS))
breast_dat$ext_pho_targets <- sapply(pho_result_list, function(x) x$ext_targets)
breast_dat$ext_pho_tar_cnt <- sapply(pho_result_list, function(x) x$ext_tar_cnt)

# Using NetPath
npa_result_list <- lapply(breast_dat$target_symbols, function(x) extend_targets_inter_database(x, inter_NPA))
breast_dat$ext_npa_targets <- sapply(npa_result_list, function(x) x$ext_targets)
breast_dat$ext_npa_tar_cnt <- sapply(npa_result_list, function(x) x$ext_tar_cnt)

# Using KEGG
kegg_ids = read.delim("kegg_cancer.txt")[,1]

downloadKPandCO <- function(drug_targets, ids) {
  #ids = sapply(ids, function(x) gsub("map", "hsa", ids))
  overlap = rep(0,length(ids))
  list_kegg_pt <- list()
  nn <- c()
  for(i in 1:length(ids)) {
    print(ids[i])
    hsap = kegg_pathway_download(ids[i], process = FALSE)
    kegg_nets = kegg_process(hsap$entries, hsap$relations)
    list_kegg_pt[[length(list_kegg_pt)+1]] = kegg_nets
    print(length(intersect(drug_targets, union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target))))
    print(length(union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target)))
    if(length(union(kegg_nets$genesymbol_source, kegg_nets$genesymbol_target)) > 0)
      nn <- c(nn, ids[i])
  }
  names(list_kegg_pt) = nn
  return(list_kegg_pt)
}

# Sample list of Ensembl gene IDs
ensembl_ids_all_drug <- unique(drug_target_ixn$Node2_ensembl_gene_id)

# Map Ensembl IDs to Gene Symbols
gene_symbols <- select(org.Hs.eg.db, 
                       keys=ensembl_ids_all_drug, 
                       columns="SYMBOL", 
                       keytype="ENSEMBL")

kegg_dat = downloadKPandCO(gene_symbols$SYMBOL, kegg_ids)

# Extend the initial targets for each KEGG data frame
kegg_result_list <- lapply(breast_dat$target_symbols, function(x) extend_targets_kegg_list(x, kegg_dat))
breast_dat$ext_kegg_targets <- sapply(kegg_result_list, function(x) x$ext_kegg_targets)
breast_dat$ext_kegg_tar_cnt <- sapply(kegg_result_list, function(x) x$ext_kegg_tar_cnt)

##################################################
####  7. Random walk with restart and FGSEA.  ####
##################################################

library(dnet)
library(fgsea)

statisticsRWR <- function(visiting_prob_matrix) {
  # Initialize empty vectors to store IQR and number of outliers for each column
  iqr_values <- numeric(ncol(visiting_prob_matrix))
  num_outliers <- numeric(ncol(visiting_prob_matrix))
  
  # Loop through each column to compute IQR and number of outliers
  for (i in 1:ncol(visiting_prob_matrix)) {
    column_data <- visiting_prob_matrix[, i]
    
    # Calculate the quartiles
    lower_quartile <- quantile(column_data, 0.25)
    upper_quartile <- quantile(column_data, 0.75)
    
    # Calculate IQR
    iqr_values[i] <- IQR(column_data)
    
    # Calculate outliers: those outside 1.5 * IQR from the quartiles
    lower_bound <- lower_quartile - 1.5 * iqr_values[i]
    upper_bound <- upper_quartile + 1.5 * iqr_values[i]
    outliers <- column_data[column_data < lower_bound | column_data > upper_bound]
    
    # Store the number of outliers
    num_outliers[i] <- length(outliers)
  }
  
  # Print IQR values
  #print(paste("IQR values for each column: ", iqr_values))
  
  # Print number of outliers
  #print(paste("Number of outliers for each column: ", num_outliers))
  
  return(list(IQRs = iqr_values, OUTs = num_outliers))
}
# Find index of maximum distance
func_RWR_threshold <- function(probabilities){
  sorted_probabilities <- sort(probabilities, decreasing = TRUE)

  # Calculate distances to line connecting first and last points
  n <- length(sorted_probabilities)
  line_endpoint1 <- c(1, sorted_probabilities[1])
  line_endpoint2 <- c(n, sorted_probabilities[n])
  line_slope <- (line_endpoint2[2] - line_endpoint1[2]) / (line_endpoint2[1] - line_endpoint1[1])
  line_intercept <- line_endpoint1[2] - line_slope * line_endpoint1[1]
  
  # Calculate perpendicular distances from each point to the line
  distances <- abs((line_slope * (1:n) - sorted_probabilities + line_intercept) / sqrt(line_slope^2 + 1))
  
  # Find index of maximum distance, this is your elbow point
  elbow_point <- which.max(distances)
  
  # To get the actual value at the elbow point
  elbow_value <- sorted_probabilities[elbow_point]
  
  # Count the number of points before the elbow point
  points_before_elbow <- elbow_point - 1
  
  # Print the elbow point and value, and the number of points before it
  #print(paste("Elbow point is at index:", elbow_point, "with value:", elbow_value))
  #print(paste("Number of points before the elbow point:", points_before_elbow))
  
  return(list(ELB = elbow_value, CNT = points_before_elbow))
}

# Convert the target sets to a matrix, each column corresponds to a drug pair
breast_seed_mat1 <- do.call(cbind, lapply(breast_dat$aggregated_genes, function(x) {
  target_set <- unlist(strsplit(x, ","))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% target_set, 1, 0))
}))

# Convert the target sets to a matrix, each column corresponds to a drug pair
breast_seed_mat2 <- do.call(cbind, lapply(breast_dat$ext_kegg_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys=target_set, 
                                     columns="ENSEMBL", 
                                     keytype="SYMBOL"))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

# Convert the target sets to a matrix, each column corresponds to a drug pair - SIGNOR
breast_seed_mat3 <- do.call(cbind, lapply(breast_dat$ext_sig_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys=target_set, 
                                     columns="ENSEMBL", 
                                     keytype="SYMBOL"))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

# Convert the target sets to a matrix, each column corresponds to a drug pair - Regulatory Network (mainly DOROTHEA)
breast_seed_mat4 <- do.call(cbind, lapply(breast_dat$ext_reg_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys=target_set, 
                                     columns="ENSEMBL", 
                                     keytype="SYMBOL"))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

# Convert the target sets to a matrix, each column corresponds to a drug pair - Phospho info
breast_seed_mat5 <- do.call(cbind, lapply(breast_dat$ext_pho_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys=target_set, 
                                     columns="ENSEMBL", 
                                     keytype="SYMBOL"))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

# Convert the target sets to a matrix, each column corresponds to a drug pair - NetPath
breast_seed_mat6 <- do.call(cbind, lapply(breast_dat$ext_npa_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys=target_set, 
                                     columns="ENSEMBL", 
                                     keytype="SYMBOL"))
  sapply(V(String_ppi_Net)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

colnames(breast_seed_mat1) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")
colnames(breast_seed_mat2) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")
colnames(breast_seed_mat3) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")
colnames(breast_seed_mat4) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")
colnames(breast_seed_mat5) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")
colnames(breast_seed_mat6) = paste(breast_dat$Drug1_DBank_id, breast_dat$Drug2_DBank_id, sep="_")

# Run RWR using dRWR
breast_result1 <- dRWR(String_ppi_Net, breast_seed_mat1, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result1) = colnames(breast_seed_mat1)
rownames(breast_result1) = rownames(breast_seed_mat1)

breast_result2 <- dRWR(String_ppi_Net, breast_seed_mat2, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result2) = colnames(breast_seed_mat2)
rownames(breast_result2) = rownames(breast_seed_mat2)

breast_result3 <- dRWR(String_ppi_Net, breast_seed_mat3, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result3) = colnames(breast_seed_mat3)
rownames(breast_result3) = rownames(breast_seed_mat3)

breast_result4 <- dRWR(String_ppi_Net, breast_seed_mat4, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result4) = colnames(breast_seed_mat4)
rownames(breast_result4) = rownames(breast_seed_mat4)

breast_result5 <- dRWR(String_ppi_Net, breast_seed_mat5, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result5) = colnames(breast_seed_mat5)
rownames(breast_result5) = rownames(breast_seed_mat5)

breast_result6 <- dRWR(String_ppi_Net, breast_seed_mat5, 
                       normalise = "row",
                       restart = 0.5, 
                       normalise.affinity.matrix = "none",
                       multicores = 4)

colnames(breast_result6) = colnames(breast_seed_mat6)
rownames(breast_result6) = rownames(breast_seed_mat6)
#
breast_rwr_stat1 <- statisticsRWR(breast_result1)
breast_rwr_stat2 <- statisticsRWR(breast_result2)
breast_rwr_stat3 <- statisticsRWR(breast_result3)
breast_rwr_stat4 <- statisticsRWR(breast_result4)
breast_rwr_stat5 <- statisticsRWR(breast_result5)
breast_rwr_stat5 <- statisticsRWR(breast_result6)
breast_rwr_stat5 <- statisticsRWR(breast_result7)
#
hist(breast_rwr_stat1$IQRs, main = "knonw drug targets - IQRs", col = "blue")
hist(breast_rwr_stat1$OUTs, main = "knonw drug targets - OUTs", col = "lightblue")
#
hist(breast_rwr_stat2$IQRs, main = "TF-based extended - IQRs", col = "red")
hist(breast_rwr_stat2$OUTs, main = "TF-based extended - OUTs", col = "pink")
#
hist(breast_rwr_stat3$IQRs, main = "KEGG-based extended - IQRs", col = "orange")
hist(breast_rwr_stat3$OUTs, main = "KEGG-based extended - OUTs", col = "yellow")
#
hist(breast_rwr_stat4$IQRs, main = "SIGNOR-based extended - IQRs", col = "darkgreen")
hist(breast_rwr_stat4$OUTs, main = "SIGNOR-based extended - OUTs", col = "green")
#
hist(breast_rwr_stat5$IQRs, main = "SL3-based extended - IQRs", col = "darkgrey")
hist(breast_rwr_stat5$OUTs, main = "SL3-based extended - OUTs", col = "grey")
#
hist(sapply(apply(breast_result1, 2, func_RWR_threshold), function(x) x$ELB))
hist(sapply(apply(breast_result2, 2, func_RWR_threshold), function(x) x$ELB))
hist(sapply(apply(breast_result3, 2, func_RWR_threshold), function(x) x$ELB))
#
hist(sapply(apply(breast_result1, 2, func_RWR_threshold), function(x) x$CNT), col = "blue")
hist(sapply(apply(breast_result2, 2, func_RWR_threshold), function(x) x$CNT), col = "red")
hist(sapply(apply(breast_result3, 2, func_RWR_threshold), function(x) x$CNT), col = "orange")

#
# Differences between consecutive probabilities
#sorted_probabilities <- sort(breast_result1[,10], decreasing = TRUE)

#plot(sorted_probabilities, type = "b", pch = 19, xlab = "Index", ylab = "Sorted Probabilities",
#     main = "Elbow Method to Determine Optimal Number")
# Add elbow point to plot
# Calculate distances to line connecting first and last points
#n <- length(sorted_probabilities)
#line_endpoint1 <- c(1, sorted_probabilities[1])
#line_endpoint2 <- c(n, sorted_probabilities[n])
#line_slope <- (line_endpoint2[2] - line_endpoint1[2]) / (line_endpoint2[1] - line_endpoint1[1])
#line_intercept <- line_endpoint1[2] - line_slope * line_endpoint1[1]
#distances <- abs((line_slope * (1:n) - sorted_probabilities + line_intercept) / sqrt(line_slope^2 + 1))
# Find index of maximum distance
#elbow_point <- which.max(distances)
# Add elbow point to plot
#points(elbow_point, sorted_probabilities[elbow_point], col = "red", pch = 4, cex = 2)

#dense_submat <- as.matrix(breast_result3)
#boxplot(dense_submat, las = 2, outline = FALSE)

library(BiocParallel)
runFGSEA <- function(rwr_results, libs, elb_points) {
  #BPPARAM <- MulticoreParam(workers = 4)
  elb_points <- sapply(apply(rwr_results, 2, func_RWR_threshold), function(x) x$ELB)
  enrichment_result_mat <- matrix(0, nrow = length(libs), ncol = ncol(rwr_results))
  rownames(enrichment_result_mat) = names(libs)
  for(i in 1:ncol(rwr_results)) {
    rankedGeneList <- sort(rwr_results[which(rwr_results[,i] > elb_points[i]),i], decreasing = TRUE)
    print(length(rankedGeneList))
    skip_iteration <- FALSE  # Initialize flag variable
    tryCatch({
              enrichment_result <- fgseaMultilevel(pathways = libs,
                                                   stats = rankedGeneList,
                                                   minSize = 5, maxSize = 500, 
                                                   scoreType = 'pos', BPPARAM = SerialParam())
              enrichment_result$NES[which(enrichment_result$pval > 0.01)] <- 0
              enrichment_result_mat[enrichment_result$pathway,i] = enrichment_result$NES
    }, error = function(e) {
      if (grepl("GSEA statistic is not defined when all genes are selected", conditionMessage(e))) {
        print(paste("Skipping iteration", i, "due to GSEA statistic error."))
        skip_iteration <- TRUE  # Skip to next iteration of the loop
      } else {
        # Handle other errors here if needed, or re-throw the error
        stop(e)
      }
    })
    if (skip_iteration) {
      next  # Skip to the next iteration of the loop
    }
  }
  return(enrichment_result_mat)
}

breast_fgsea1 <- runFGSEA(breast_result1, c(breast_libs, adr_libs))
breast_fgsea2 <- runFGSEA(breast_result2, c(breast_libs, adr_libs))
breast_fgsea3 <- runFGSEA(breast_result3, c(breast_libs, adr_libs))
breast_fgsea4 <- runFGSEA(breast_result4, c(breast_libs, adr_libs))
breast_fgsea5 <- runFGSEA(breast_result5, c(breast_libs, adr_libs))
breast_fgsea6 <- runFGSEA(breast_result6, c(breast_libs, adr_libs))

######################################################################
####  8. Visaulizations of the created DrugBank/FIMM data space.  ####
######################################################################

library(FactoMineR)
library(factoextra)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

pca_breast_lib1 = PCA(t(breast_fgsea2[-which(apply(breast_fgsea2,1,var) < 0.01),]), 
                      scale.unit = FALSE, ncp = 5, graph = FALSE)
fviz_pca_ind(pca_breast_lib1, axes = c(1, 2), col.ind = breast_dat$drug_class, label = "none")
fviz_pca_ind(pca_breast_lib1, axes = c(1, 2), col.ind = as.factor(breast_dat$totSyn), label = "none")


plotNESBars <- function(nes_result, category, title) {
  
  df <- as.data.frame(t(nes_result[-which(apply(nes_result,1,var) < 0.01),]))
  
  df_long <- pivot_longer(df, cols = colnames(df), 
                          names_to = "ScoreType", values_to = "Value")
  
  df_long$categories <- rep(category, ncol(df))
  
  summary_stats <- df_long %>%
    group_by(ScoreType, categories) %>%
    summarise(mean_value = mean(Value), se = sd(Value) / sqrt(n()), .groups = 'drop')
  
  summary_stats = summary_stats[which(summary_stats$categories != "NA"),]
  
  print(summary_stats)
  
  # Wrap the long names
  summary_stats$ScoreType <- str_wrap(summary_stats$ScoreType, width = 20)
  
  # Plot using ggplot2
  gp <- ggplot(summary_stats, aes(x = categories, y = mean_value, fill = categories)) +
               geom_bar(stat = "identity", position = "dodge") +
               geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.2, position = position_dodge(0.9)) +
               facet_wrap(~ ScoreType) +
               labs(y = "Mean Value +/- SE", title = title) +
               theme(strip.text = element_text(size = 8))

  return(gp)
}

plotNESBars2 <- function(nes_result, category, title) {
  
  print(length(which(apply(nes_result,1,var) < 0.01)))
  
  df <- as.data.frame(t(nes_result[-which(apply(nes_result,1,var) < 0.01),]))
  
  df_long <- pivot_longer(df, cols = colnames(df), 
                          names_to = "ScoreType", values_to = "Value")
  
  df_long$categories <- rep(category, ncol(df))
  
  summary_stats <- df_long %>%
    group_by(ScoreType, categories) %>%
    summarise(mean_value = mean(Value), se = sd(Value) / sqrt(n()), .groups = 'drop')
  
  # Wrap the long names
  summary_stats$ScoreType <- str_wrap(summary_stats$ScoreType, width = 20)
  
  # Assuming categories are factors, convert to numeric for plotting trend
  summary_stats$categories_numeric <- as.numeric(as.factor(summary_stats$categories))
  
  # Plot using ggplot2
  gp <- ggplot(summary_stats, aes(x = categories, y = mean_value, fill = categories)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.2, position = position_dodge(0.9)) +
    geom_line(aes(x = categories_numeric, group = ScoreType), color = "red") + # Trend line
    facet_wrap(~ ScoreType) +
    labs(y = "Mean Value +/- SE", title = title) +
    theme(strip.text = element_text(size = 8))
  
  return(gp)
}

plotNESBars3 <- function(nes_result, category, title) {
  
  df <- as.data.frame(t(nes_result[-which(apply(nes_result, 1, var) < 0.01),]))
  
  df_long <- pivot_longer(df, cols = colnames(df), 
                          names_to = "ScoreType", values_to = "Value")
  
  df_long$categories <- rep(category, ncol(df))
  
  summary_stats <- df_long %>%
    group_by(ScoreType, categories) %>%
    summarise(mean_value = mean(Value), se = sd(Value) / sqrt(n()), .groups = 'drop')
  
  # Filter out NA categories
  summary_stats = summary_stats[which(summary_stats$categories != "NA"),]
  
  print(summary_stats)
  
  gp <- ggbarplot(summary_stats, x = "categories", y = "mean_value", 
                  fill = "categories", 
                  facet.by = "ScoreType",
                  position = position_dodge(0.9),
                  error.plot = "errorbar",
                  color = "white",
                  palette = "jco",
                  title = title,
                  ylab = "Mean Value +/- SE",
                  add = "mean_se")
  
  return(gp)
}

plotNESBars2(breast_fgsea1[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - direct targets")
plotNESBars2(breast_fgsea2[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - kegg-based ext")
plotNESBars2(breast_fgsea3[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - signor")
plotNESBars2(breast_fgsea4[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - dorothea")
plotNESBars2(breast_fgsea5[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - phospho")
plotNESBars2(breast_fgsea6[names(breast_libs),], as.factor(breast_dat$totSyn), title = "Breast/syn - netpath")

plotNESBars(breast_fgsea1[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - direct targets")
plotNESBars(breast_fgsea2[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - kegg-based ext")
plotNESBars(breast_fgsea3[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - signor")
plotNESBars(breast_fgsea4[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - dorothea")
plotNESBars(breast_fgsea5[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - phospho")
plotNESBars(breast_fgsea6[names(breast_libs),], as.factor(breast_dat$phamk), title = "Breast/pk - netpath")

plotNESBars2(breast_fgsea1[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - direct targets")
plotNESBars2(breast_fgsea2[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - kegg-based ext")
plotNESBars2(breast_fgsea3[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - signor")
plotNESBars2(breast_fgsea4[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - signal link 3")
plotNESBars2(breast_fgsea5[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - dorothea")
plotNESBars2(breast_fgsea6[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - phospho")
plotNESBars2(breast_fgsea7[names(adr_libs),], as.factor(breast_dat$totSyn), title = "Breast/adr - adv - netpath")

plotNESBars(breast_fgsea1[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - direct targets")
plotNESBars(breast_fgsea2[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - kegg-based ext")
plotNESBars(breast_fgsea3[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - signor")
plotNESBars(breast_fgsea4[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - dorothea")
plotNESBars(breast_fgsea5[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - phospho")
plotNESBars(breast_fgsea6[names(adr_libs),], breast_dat$breastADV, title = "Breast/adr2 - netpath")


# Our defintion of good and bad drugs

id_breast_good = which(breast_dat$totSyn >= 3 & breast_dat$breastADV == "Unk")
id_breast_bad = which(breast_dat$totSyn <= -3 & breast_dat$breastADV == "Adv")
breast_dat$fc = rep("NA", nrow(breast_dat))
breast_dat$fc[id_breast_good] = "Good"
breast_dat$fc[id_breast_bad] = "Bad"
table(breast_dat$fc)

plotNESBars(breast_fgsea1[names(breast_libs),], breast_dat$fc, title = "Breast/syn - direct targets")
plotNESBars(breast_fgsea2[names(breast_libs),], breast_dat$fc, title = "Breast/syn - kegg-based ext")
plotNESBars(breast_fgsea3[names(breast_libs),], breast_dat$fc, title = "Breast/syn - signor")
plotNESBars(breast_fgsea5[names(breast_libs),], breast_dat$fc, title = "Breast/syn - dorothea")
plotNESBars(breast_fgsea6[names(breast_libs),], breast_dat$fc, title = "Breast/syn - phospho")
plotNESBars(breast_fgsea7[names(breast_libs),], breast_dat$fc, title = "Breast/syn - netpath")

plotNESBars(breast_fgsea1[names(adr_libs),], breast_dat$fc, title = "Breast/adr - direct targets")
plotNESBars(breast_fgsea2[names(adr_libs),], breast_dat$fc, title = "Breast/adr - kegg-based ext")
plotNESBars(breast_fgsea3[names(adr_libs),], breast_dat$fc, title = "Breast/adr - signor")
plotNESBars(breast_fgsea5[names(adr_libs),], breast_dat$fc, title = "Breast/adr - dorothea")
plotNESBars(breast_fgsea6[names(adr_libs),], breast_dat$fc, title = "Breast/adr - phospho")
plotNESBars(breast_fgsea7[names(adr_libs),], breast_dat$fc, title = "Breast/adr - netpath")
