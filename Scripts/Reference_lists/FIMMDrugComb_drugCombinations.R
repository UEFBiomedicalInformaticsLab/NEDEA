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
write.csv(FimmDrugComb_cellLine, "Databases/FimmDrugComb/cell_lines.csv", quote = TRUE, row.names = FALSE)

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
write.csv(FimmDrugComb_drugs, "Databases/FimmDrugComb/drugs.csv", quote = TRUE, row.names = FALSE)


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

ddi = read.csv("/Users/vittfo/Downloads/drug_drug_interactions.csv") # we need to read the drug-drug infor file from DrugBank
colnames(ddi)[c(1,4)] <- c("Drug1_DBank_id", "Drug2_DBank_id") # pelase rename these column names as you wish

# The following code section was used to study the most frequentely reported set of words (item set) in the strings "description" 
# matching a given key word. 
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
                              
# now we can define our labels
                      
# Indicate the drug pairs associated with a confirmed increase or decrease in efficacy (of which there are very few cases).                       
ddi$effVsNeff = rep('Unk', nrow(ddi))
ddi$effVsNeff[test_eff_inc] = "IncEff"
ddi$effVsNeff[test_eff_dec] = "DecEff"

# Using pharmacokinetics information to identify a potential decrease or increase in efficacy for each drug pair
ddi$phamk = rep('Unk', nrow(ddi))
ddi$phamk[test_meinc | test_scdec | test_absde] = "PDecEff"
ddi$phamk[test_medec | test_scinc] = "PIncEff"

# When analyzing drug-drug interactions (DDIs), understanding the potential implications of these interactions on drug efficacy 
# is crucial for therapeutic decisions. Your strategy uses three key pharmacokinetic parameters to infer potential alterations in drug efficacy:
# - Serum Concentration Changes: Alterations in the serum concentration of a drug can directly affect its therapeutic outcomes. 
#   An increase in serum concentration might enhance its efficacy, while a decrease could reduce its therapeutic effect. 
#   By observing these changes, you can reasonably speculate on the potential impact on drug efficacy.
# - Metabolism Alterations: Changes in the metabolism of a drug can indirectly influence its serum concentration. 
#   An increase in metabolism generally results in decreased serum concentration, leading to potentially reduced efficacy. 
#   Conversely, decreased metabolism can lead to increased serum concentrations, which may enhance efficacy but also the risk of adverse effects.
# - Absorption Decreases: Absorption is the gateway to systemic circulation for orally administered drugs. 
#   A decrease in absorption directly affects the amount of the drug entering systemic circulation, 
#   which can reduce its therapeutic effect.

# Now we select toxicity endpoints which are often observed in cancer for each drug pair 

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
ovararian_tissue = which(fimm_db_all$tissue == "OVARY")
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
table(fimm_db_all$totSyn[ovararian_tissue], fimm_db_all$phamk[ovararian_tissue])
table(fimm_db_all$totSyn[ovararian_tissue], fimm_db_all$ovarianADV[ovararian_tissue])
print('<skin>')
table(fimm_db_all$totSyn[skin_tissue], fimm_db_all$phamk[skin_tissue])
table(fimm_db_all$totSyn[skin_tissue], fimm_db_all$skinADV[skin_tissue])
#

# how we select the good drug combinations

# BreastCancer specific
breast_dp = fimm_db_all[breast_tissue,]
# good drug combination
length(which((breast_dp$totSyn >= 3 & breast_dp$breastADV == "Unk") | 
             (breast_dp$phamk == "PIncEff" & breast_dp$breastADV == "Unk" & breast_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
             (breast_dp$effVsNeff== "IncEff" & breast_dp$breastADV == "Unk" & breast_dp$totSyn >= -2)))
# bad drug combination
length(which((breast_dp$totSyn <= -3 & breast_dp$breastADV == "Adv") | 
             (breast_dp$phamk == "PDecEff" & breast_dp$breastADV == "Unk" & breast_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
             (breast_dp$effVsNeff== "DecEff" & breast_dp$breastADV == "Unk" & breast_dp$totSyn <= -3)))


# LungCancer specific
lung_dp = fimm_db_all[lung_tissue,]
# good drug combination
length(which((lung_dp$totSyn >= 3 & lung_dp$lungADV == "Unk") | 
               (lung_dp$phamk == "PIncEff" & lung_dp$lungADV == "Unk" & lung_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
               (lung_dp$effVsNeff== "IncEff" & lung_dp$lungADV == "Unk" & lung_dp$totSyn >= -2)))
# bad drug combination
length(which((lung_dp$totSyn <= -3 & lung_dp$lungADV == "Adv") | 
               (lung_dp$phamk == "PDecEff" & lung_dp$lungADV == "Unk" & lung_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
               (lung_dp$effVsNeff== "DecEff" & lung_dp$lungADV == "Unk" & lung_dp$totSyn <= -3)))


# ProstateCancer specific
prostate_dp = fimm_db_all[prostate_tissue,]
# good drug combination
length(which((prostate_dp$totSyn >= 3 & prostate_dp$prostateADV == "Unk") | 
               (prostate_dp$phamk == "PIncEff" & prostate_dp$prostateADV == "Unk" & prostate_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
               (prostate_dp$effVsNeff== "IncEff" & prostate_dp$prostateADV == "Unk" & prostate_dp$totSyn >= -2)))
# bad drug combination
length(which((prostate_dp$totSyn <= -3 & prostate_dp$prostateADV == "Adv") | 
               (prostate_dp$phamk == "PDecEff" & prostate_dp$prostateADV == "Unk" & prostate_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
               (prostate_dp$effVsNeff== "DecEff" & prostate_dp$prostateADV == "Unk" & prostate_dp$totSyn <= -3)))


# KidneyCancer specific
kidney_dp = fimm_db_all[kidney_tissue,]
# good drug combination
length(which((kidney_dp$totSyn >= 3 & kidney_dp$kidneyADV == "Unk") | 
               (kidney_dp$phamk == "PIncEff" & kidney_dp$kidneyADV == "Unk" & kidney_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
               (kidney_dp$effVsNeff== "IncEff" & kidney_dp$kidneyADV == "Unk" & kidney_dp$totSyn >= -2)))
# bad drug combination
length(which((kidney_dp$totSyn <= -3 & kidney_dp$kidneyADV == "Adv") | 
               (kidney_dp$phamk == "PDecEff" & kidney_dp$kidneyADV == "Unk" & kidney_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
               (kidney_dp$effVsNeff== "DecEff" & kidney_dp$kidneyADV == "Unk" & kidney_dp$totSyn <= -3)))


# OvarianCancer specific
ovarian_dp = fimm_db_all[ovarian_tissue,]
# good drug combination
length(which((ovarian_dp$totSyn >= 3 & ovarian_dp$ovarianADV == "Unk") | 
               (ovarian_dp$phamk == "PIncEff" & ovarian_dp$ovarianADV == "Unk" & ovarian_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
               (ovarian_dp$effVsNeff== "IncEff" & ovarian_dp$ovarianADV == "Unk" & ovarian_dp$totSyn >= -2)))
# bad drug combination
length(which((ovarian_dp$totSyn <= -3 & ovarian_dp$ovarianADV == "Adv") | 
               (ovarian_dp$phamk == "PDecEff" & ovarian_dp$ovarianADV == "Unk" & ovarian_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
               (ovarian_dp$effVsNeff== "DecEff" & ovarian_dp$ovarianADV == "Unk" & ovarian_dp$totSyn <= -3)))


# SkinCancer specific
skin_dp = fimm_db_all[skin_tissue,]
# good drug combination
length(which((skin_dp$totSyn >= 3 & skin_dp$skinADV == "Unk") | 
               (skin_dp$phamk == "PIncEff" & skin_dp$skinADV == "Unk" & skin_dp$totSyn >= -2) |     # there is indecision in the synergistic score or it is positive
               (skin_dp$effVsNeff== "IncEff" & skin_dp$skinADV == "Unk" & skin_dp$totSyn >= -2)))
# bad drug combination
length(which((skin_dp$totSyn <= -3 & skin_dp$skinADV == "Adv") | 
               (skin_dp$phamk == "PDecEff" & skin_dp$skinADV == "Unk" & skin_dp$totSyn <= -3) |     # there is indecision in the synergistic score or it is positive
               (skin_dp$effVsNeff== "DecEff" & skin_dp$skinADV == "Unk" & skin_dp$totSyn <= -3)))


# then you can add the code to split the dataset based on the tissue columns

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
