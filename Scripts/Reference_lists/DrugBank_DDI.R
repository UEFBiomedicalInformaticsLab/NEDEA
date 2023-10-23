set.seed(5081)




# Script to analyse the drug-drug interactions (DDIs) reported in DrugBank
# Notes:
# (a) Identifies the most frequent type of DDIs
# (b) Literature information on cancer type specific adverse reactions used for drug combination classification




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


# Read the DrugBank drug-drug interactions
# DrugBank_ddi <- read.csv("Databases/DrugBank/drug_drug_interactions.csv")
DrugBank_ddi <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_ddi <- DrugBank_ddi$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")


# Define function to identify the top k frequent word sequence of length n
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


# Define function to identify the top k frequent word sequence of length n containing certain keyword
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




# Obtain the most frequent word sequences
top_freq <- find_top_k_frequent_sequences(DrugBank_ddi$description, n = 4, k = 100)
risk_sev_res <- find_top_k_sequences_with_keyword(DrugBank_ddi$description, 
                                                  keyword = "The risk or severity", n = 6, k = 50)



test_eff_inc <- grepl("The therapeutic efficacy of", DrugBank_ddi$description) & grepl("can be increased when used in combination", DrugBank_ddi$description)
test_eff_dec <- grepl("The therapeutic efficacy of", DrugBank_ddi$description) & grepl("can be decreased when used in combination", DrugBank_ddi$description)
test_anticinc <- grepl("may increase the anticoagulant activities", DrugBank_ddi$description)
test_anticdec <- grepl("may decrease the anticoagulant activities", DrugBank_ddi$description)
test_antihinc <- grepl("may increase the antihypertensive activities", DrugBank_ddi$description)
test_antihdec <- grepl("may decrease the antihypertensive activities", DrugBank_ddi$description)
test_nsd <- grepl("(CNS depressant)", DrugBank_ddi$description)
test_scinc <- grepl("The serum concentration of", DrugBank_ddi$description) & grepl("can be increased", DrugBank_ddi$description)
test_scdec <- grepl("The serum concentration of", DrugBank_ddi$description) & grepl("can be decreased", DrugBank_ddi$description)
test_meinc <- grepl("The metabolism of", DrugBank_ddi$description) & grepl("can be increased", DrugBank_ddi$description)
test_medec <- grepl("The metabolism of", DrugBank_ddi$description) & grepl("can be decreased", DrugBank_ddi$description)
test_absde <- grepl("a decrease in the absorption ", DrugBank_ddi$description)
test_all_risk <- grepl("risk or severity", DrugBank_ddi$description)

list_tox_test = list()
for(i in 1:length(risk_sev_res$ngrams)){
  list_tox_test[[length(list_tox_test) + 1]] <- grepl(risk_sev_res$ngrams[i], DrugBank_ddi$description)
}
names(list_tox_test) <- risk_sev_res$ngrams




# sum(test_eff_inc)
# sum(test_eff_dec)
# sum(test_anticinc)
# sum(test_anticdec)
# sum(test_antihinc)
# sum(test_antihdec)
# sum(test_nsd)
# sum(test_scinc)
# sum(test_scdec)
# sum(test_meinc)
# sum(test_medec)
# sum(test_absde)
# sum(test_all_risk)




# Categorize drug combinations based on increased/decreased therapeutic efficacy
DrugBank_ddi$class_therapeuticEfficacy <- rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$class_therapeuticEfficacy[test_eff_inc] <- "increased"
DrugBank_ddi$class_therapeuticEfficacy[test_eff_dec] <- "decreased"


# Categorize drug combinations based on effect of metabolism
DrugBank_ddi$class_metabolicEffect <- rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$class_metabolicEffect[test_meinc | test_scdec | test_absde] <- "decreased"
DrugBank_ddi$class_metabolicEffect[test_medec | test_scinc] <- "increased"





# Literature based identification of cancer specific adverse reactions

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

set_adv_ids <- c(1, 2, 4, 6, 10, 11, 12, 15, 16, 20, 28, 33, 35, 39, 48, 49) 





#### Breast cancer adverse effects
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of peripheral neuropathy 
# The risk or severity of electrolyte imbalance

names(list_tox_test)[c(set_adv_ids, 25, 31, 47)]
sapply(list_tox_test[c(set_adv_ids, 25, 31, 47)], sum)
sum(Reduce("|",list_tox_test[c(set_adv_ids, 25, 31, 47)]))
DrugBank_ddi$ADR_BreastCancer <- rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_BreastCancer[Reduce("|",list_tox_test[c(set_adv_ids, 25, 47)])] <- "adr_positive"





#### Lung cancer adverse effects
# The risk or severity of bleeding
# The risk or severity of Thrombosis
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of respiratory
# The risk or severity of hyponatremia
# The risk or severity of electrolyte imbalance

names(list_tox_test)[c(set_adv_ids, 3, 19, 25, 31, 34, 46)]
sapply(list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)], sum)
sum(Reduce("|",list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)]))
DrugBank_ddi$ADR_LungCancer = rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_LungCancer[Reduce("|",list_tox_test[c(set_adv_ids, 3, 19, 25, 31, 34, 46)])] <- "adr_positive"






#### Prostate cancer adverse effects
# The risk or severity of fluid retention
# The risk or severity of urinary retention
# The risk or severity of hypercalcemia

names(list_tox_test)[c(set_adv_ids, 31, 38, 42, 43)]
sapply(list_tox_test[c(set_adv_ids, 31, 38, 42, 43)], sum)
sum(Reduce("|",list_tox_test[c(set_adv_ids, 31, 38, 42, 43)]))
DrugBank_ddi$ADR_ProstateCancer = rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_ProstateCancer[Reduce("|",list_tox_test[c(set_adv_ids, 31, 38, 42, 43)])] <- "adr_positive"                





#### Ovarian cancer adverse effects  
# The risk or severity of bleeding
# The risk or severity of Cardiac Arrhythmia
# The risk or severity of urinary retention
# The risk or severity of hypercalcemia
# The risk or severity of respiratory

names(list_tox_test)[c(set_adv_ids, 3, 25, 42, 43, 46)]
sapply(list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)], sum)
sum(Reduce("|",list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)]))
DrugBank_ddi$ADR_OvaryCancer = rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_OvaryCancer[Reduce("|",list_tox_test[c(set_adv_ids, 3, 25, 42, 43, 46)])] <- "adr_positive"                






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
sum(Reduce("|",list_tox_test[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)]))
DrugBank_ddi$ADR_KidneyCancer = rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_KidneyCancer[Reduce("|",list_tox_test[c(set_adv_ids, 3, 5, 8, 31, 34, 38, 42, 43)])] <- "adr_positive"   






#### Skin cancer adverse effects  
# The risk or severity of respiratory: if skin cancer metastasizes to the lungs, some treatments might have respiratory side effects.

names(list_tox_test)[c(set_adv_ids, 46)]
sapply(list_tox_test[c(set_adv_ids, 46)], sum)
sum(Reduce("|",list_tox_test[c(set_adv_ids, 46)]))
DrugBank_ddi$ADR_SkinCancer = rep("unknown", nrow(DrugBank_ddi))
DrugBank_ddi$ADR_SkinCancer[Reduce("|",list_tox_test[c(set_adv_ids, 46)])] <- "adr_positive"   




if(!dir.exists("InputFiles/Reference_list/"))dir.create("InputFiles/Reference_list/", recursive = TRUE)
saveRDS(DrugBank_ddi, "InputFiles/Reference_list/DrugBank_DDI_processed.rds")




print(warnings())