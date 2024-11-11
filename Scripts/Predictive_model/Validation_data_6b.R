set.seed(5081)



# Generate validation data 6b
# Notes:
# (a) Based on ctrdata



# Load libraries
library(unixtools)
library(optparse)
library(fuzzyjoin)
library(stringdist)
library(tidyverse)
library(ctrdata)
library(igraph)
library(org.Hs.eg.db)
library(OmnipathR)
library(readxl)
source("Scripts/Functions/Functions_drug_target.R")


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$disease %in% c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  print_help(opt_parser)
  stop("--disease must be one of the following: BreastCancer, KidneyCancer, LungCancer, OvaryCancer, ProstateCancer, SkinCancer", call. = FALSE)
}


# Define global options for this script
disease <- opt$disease


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease, "\n"))


#####


# Read the ATC codes of the drugs from Drug Bank
DrugBank_drug_ATC <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_drug_ATC$drugs$atc_codes
colnames(DrugBank_drug_ATC) <- gsub("drugbank-id", "DrugBank_drug_ID", colnames(DrugBank_drug_ATC))

# Read the drug combination used in the training 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")


# Add the ATC codes at level_1
# Using many-to-many mapping to map all possible ATC codes to a single drug
train_drugCombs_cat <- train_drugCombs_cat %>%
  left_join(DrugBank_drug_ATC %>%
              dplyr::select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug1_ATC_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  left_join(DrugBank_drug_ATC %>%
              dplyr::select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug2_ATC_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  distinct()

possible_ATC <- as.character(na.exclude(unique(c(train_drugCombs_cat$Drug1_ATC_code_1, train_drugCombs_cat$Drug2_ATC_code_1))))


#####


# Generate the database with the predetermined quury

query <- switch(disease,
                "BreastCancer" = "https://clinicaltrials.gov/search?cond=Breast%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int",
                "KidneyCancer" = "https://clinicaltrials.gov/search?cond=Kidney%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int",
                "LungCancer" = "https://clinicaltrials.gov/search?cond=Lung%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int",
                "OvaryCancer" = "https://clinicaltrials.gov/search?cond=Ovarian%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int",
                "ProstateCancer" = "https://clinicaltrials.gov/search?cond=Prostate%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int",
                "SkinCancer" = "https://clinicaltrials.gov/search?cond=Skin%20Cancer&aggFilters=phase:4%203%202,results:with,status:com%20ter,studyType:int"
)


query <- ctrGetQueryUrl(url = query)


db <- nodbi::src_sqlite(
  dbname = paste0("tmp_dir/6b_", disease, "_clinicalTrials_sql"),
  collection = disease
)

ctrLoadQueryIntoDb(
  queryterm = query,
  euctrresults = FALSE,
  con = db
)


#####


# Extract the intervention information

result_intervention <- dbGetFieldsIntoDf(fields = "protocolSection.armsInterventionsModule.armGroups",
                                         con = db)

result_intervention$protocolSection.armsInterventionsModule.armGroups <- lapply(result_intervention$protocolSection.armsInterventionsModule.armGroups, function(x){
  
  if(is.data.frame(x)){
    
    if(! all(sapply(x$interventionNames, is.character))){
      
      x$interventionNames <- sapply(x$interventionNames, as.character)
      x <- x[which(sapply(x$interventionNames, length) != 0), ] # Removing arms without interventions to remove data precessing error
      
    }
    
    x <- unnest(x, cols = colnames(x)) %>% 
      group_by(label) %>% 
      summarise(across(everything(), ~ toString(unique(.x))))
    
  } else {
    x <- as_tibble(x)
    if(nrow(x) > 1){
      x <- x %>% 
        group_by(label) %>% 
        summarise(across(everything(), ~ toString(unique(.x))))
    }
  }
  x
})

result_intervention <- result_intervention %>%
  unnest(protocolSection.armsInterventionsModule.armGroups)

# Keep only trial arms that contain interventions/treatments
result_intervention <- result_intervention %>% filter(type %in% c("ACTIVE_COMPARATOR", "EXPERIMENTAL"))


# Filter trials/interventions with only two drugs
intervention_types <- c("Drug", "Behavioral", "Biological", "Device", 
                        "Dietary", "Other", "Procedure", "Radiation")

for (type in intervention_types) {
  count_column <- paste0(tolower(type), "_count")
  result_intervention[[count_column]] <- sapply(result_intervention$interventionNames, function(x) {
    str_count(x, paste0(type, ":"))
  }, USE.NAMES = FALSE)
}

result_intervention <- result_intervention %>% 
  filter(drug_count == 2 &
           behavioral_count == 0 &
           biological_count == 0 &
           device_count == 0 &
           dietary_count == 0 &
           other_count == 0 &
           procedure_count == 0 &
           radiation_count == 0) %>%
  dplyr::select(!ends_with("_count"))

# Remove arms with placebo
result_intervention <- result_intervention[grep(pattern = "placebo", 
                                                x = result_intervention$interventionNames, 
                                                ignore.case = TRUE, 
                                                invert = TRUE), ]

# #####
# 
# 
# # Extract the groups information
# result_groups <- dbGetFieldsIntoDf(fields = c("resultsSection.participantFlowModule.groups"),
#                                    con = db)
# 
# result_groups$col_df <- unlist(lapply(result_groups$resultsSection.participantFlowModule.groups, is.data.frame))
# 
# 
# tmp1 <- result_groups %>% 
#   filter(col_df %in% "TRUE") %>% 
#   unnest(resultsSection.participantFlowModule.groups)%>%
#   # select(!col_df) %>% 
#   dplyr::select("_id", "id", "title", "description")
# 
# tmp2 <- result_groups %>% 
#   filter(col_df %in% "FALSE") %>% 
#   mutate("resultsSection.participantFlowModule.groups" = lapply(resultsSection.participantFlowModule.groups, as_tibble) )%>% 
#   unnest(resultsSection.participantFlowModule.groups) %>%
#   # select(!col_df) %>% 
#   dplyr::select("_id", "id", "title", "description")
# 
# result_groups <- rbind(tmp1, tmp2)
# rm(list = c("tmp1", "tmp2"))
# 
# 
# #####


# Extract the adverse event information
result_adverseEvents <- dbGetFieldsIntoDf(fields = c("resultsSection.adverseEventsModule.eventGroups"),
                                          con = db)

result_adverseEvents <- result_adverseEvents %>% filter(`_id` %in% result_intervention$`_id`) # Keep only trial of prefiltered interventions

result_adverseEvents$col_df <- unlist(lapply(result_adverseEvents$resultsSection.adverseEventsModule.eventGroups, is.data.frame))

tmp1 <- result_adverseEvents %>% 
  filter(col_df %in% "TRUE") %>% 
  unnest(resultsSection.adverseEventsModule.eventGroups)%>%
  # select(!col_df) %>% 
  dplyr::select(c("_id", "id", "title", "description", 
                  "deathsNumAffected", "deathsNumAtRisk", 
                  "seriousNumAffected", "seriousNumAtRisk"))

tmp2 <- result_adverseEvents %>% 
  filter(col_df %in% "FALSE") %>% 
  mutate("resultsSection.adverseEventsModule.eventGroups" = lapply(resultsSection.adverseEventsModule.eventGroups, as_tibble) ) %>% 
  unnest(resultsSection.adverseEventsModule.eventGroups) %>%
  # select(!col_df) %>% 
  dplyr::select(c("_id", "id", "title", "description", 
                  "deathsNumAffected", "deathsNumAtRisk", 
                  "seriousNumAffected", "seriousNumAtRisk"))

result_adverseEvents <- rbind(tmp1, tmp2)
rm(list = c("tmp1", "tmp2"))


#####


# Merge the data

# valid_drugCombs_cat <- result_adverseEvents %>% 
#   dplyr::select(!c("id")) %>% 
#   left_join(result_intervention, by = c("_id" = "_id", "title" = "label", "description" = "description")) %>%
#   filter(!is.na(interventionNames))

# valid_drugCombs_cat <- fuzzy_left_join(x = result_adverseEvents %>% dplyr::select(-id), 
#                                        y = result_intervention, 
#                                        by = c("_id" = "_id", "title" = "label", "description" = "description"), 
#                                        match_fun = list("_id" = function(t1,  t2){ t1 == t2 }, 
#                                                         "title" = function(t1, t2){ stringdist(t1, t2, method = "jw") < 0.2 }, 
#                                                         "description" = function(t1, t2){ stringdist(t1, t2, method = "jw") < 0.2 })
# ) %>% 
#   mutate( title_dist = stringdist(title, label, method = "jw"),
#           desc_dist = stringdist(`description.x`, `description.y`, method = "jw")) %>%
#   group_by(`_id.x`, `title`, `description.x`) %>%
#   slice_min(order_by = title_dist + desc_dist, n = 1) %>%
#   ungroup() %>%
#   dplyr::select(-title_dist, -desc_dist, -`description.y`, -`_id.y`) %>%
#   dplyr::rename("description" = "description.x", "_id" = "_id.x") %>%
#   filter(!is.na(interventionNames))



valid_drugCombs_cat <- fuzzy_left_join(x = result_adverseEvents %>% dplyr::select(-id), 
                                       y = result_intervention, 
                                       by = c("_id" = "_id", "title" = "label"), 
                                       match_fun = list("_id" = function(t1,  t2){ t1 == t2 }, 
                                                        "title" = function(t1, t2){ stringdist(t1, t2, method = "jw") < 0.2 })
) %>% 
  mutate( title_dist = stringdist(title, label, method = "jw") ) %>%
  group_by(`_id.x`, `title`) %>%
  slice_min(order_by = title_dist, n = 1) %>%
  ungroup() %>%
  dplyr::select(-title_dist, -`description.y`, -`_id.y`) %>%
  dplyr::rename("description" = "description.x", "_id" = "_id.x") %>%
  filter(!is.na(interventionNames))


# valid_drugCombs_cat <- valid_drugCombs_cat[grep(pattern = "Behavioral|Biological|Device|Dietary|Other|Procedure|Radiation", 
#                                                 x = valid_drugCombs_cat$interventionNames, 
#                                                 ignore.case = TRUE, 
#                                                 invert = FALSE), ]

# valid_drugCombs_cat$keep <- sapply(valid_drugCombs_cat$interventionNames, function(x){  str_count(x,  "Drug:")  }, 
#                                    USE.NAMES = FALSE)
# 
# valid_drugCombs_cat <- valid_drugCombs_cat[valid_drugCombs_cat$keep == 2, ]
# 
# 
# valid_drugCombs_cat <- valid_drugCombs_cat[grep(pattern = "placebo", 
#                                                 x = valid_drugCombs_cat$interventionNames, 
#                                                 ignore.case = TRUE, 
#                                                 invert = TRUE), ]

valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  separate(interventionNames, into = c("Drug1_name", "Drug2_name"), sep = ", Drug: ") 

valid_drugCombs_cat$Drug1_name <- gsub("Drug: ", "", valid_drugCombs_cat$Drug1_name)
# valid_drugCombs_cat$Drug2_name <- gsub("Drug: ", "", valid_drugCombs_cat$Drug2_name)

valid_drugCombs_cat$Drug1_name <- toupper(valid_drugCombs_cat$Drug1_name)
valid_drugCombs_cat$Drug2_name <- toupper(valid_drugCombs_cat$Drug2_name)


#####


valid_drugCombs_cat[is.na(valid_drugCombs_cat$deathsNumAffected), "deathsNumAffected"] <- 0
valid_drugCombs_cat[is.na(valid_drugCombs_cat$deathsNumAtRisk), "deathsNumAtRisk"] <- 0
valid_drugCombs_cat[is.na(valid_drugCombs_cat$seriousNumAffected), "seriousNumAffected"] <- 0
valid_drugCombs_cat[is.na(valid_drugCombs_cat$seriousNumAtRisk), "seriousNumAtRisk"] <- 0

valid_drugCombs_cat <- valid_drugCombs_cat[valid_drugCombs_cat$seriousNumAtRisk > 10 | valid_drugCombs_cat$seriousNumAtRisk > 10, ]

valid_drugCombs_cat$death_ratio <- (valid_drugCombs_cat$deathsNumAffected/valid_drugCombs_cat$deathsNumAtRisk) * 100
valid_drugCombs_cat$seriousAE_ratio <- (valid_drugCombs_cat$seriousNumAffected/valid_drugCombs_cat$seriousNumAtRisk) * 100


valid_drugCombs_cat[valid_drugCombs_cat$death_ratio %in% "NaN", "death_ratio"] <- 0
valid_drugCombs_cat[valid_drugCombs_cat$seriousAE_ratio %in% "NaN", "seriousAE_ratio"] <- 0


valid_drugCombs_cat$class_EffAdv <- NA
valid_drugCombs_cat$class_EffAdv <- ifelse( (valid_drugCombs_cat$death_ratio > 25 | valid_drugCombs_cat$seriousAE_ratio > 25)  , "Adv", "Eff")


#####


# Use the drug names to map to the Drugbank drug IDs

DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")

DrugBank_drug_info <- DrugBank_data$drugs$general_information
DrugBank_drug_info <- DrugBank_drug_info[, c("primary_key", "name", "type")]
colnames(DrugBank_drug_info)[1] <- "DrugBank_id"

DrugBank_drug_info <- DrugBank_drug_info %>% 
  left_join(DrugBank_data$drugs$external_identifiers %>% 
              filter(!resource %in% "Wikipedia") %>% 
              dplyr::select(c("parent_key", "identifier")), 
            by = c("DrugBank_id" = "parent_key")) %>% 
  left_join(DrugBank_data$drugs$synonyms %>% 
              dplyr::select(c("drugbank-id", "synonym")), 
            by = c("DrugBank_id" = "drugbank-id"), 
            relationship = "many-to-many") %>%
  left_join(DrugBank_data$products %>% 
              dplyr::select(c("parent_key", "name")) %>% 
              distinct(), 
            by = c("DrugBank_id" = "parent_key"), 
            relationship = "many-to-many")

DrugBank_drug_info <- DrugBank_drug_info %>% 
  pivot_longer(cols = c("name.x", "name.y", "identifier", "synonym"), 
               names_to = "id_type", 
               values_to = "names") %>%
  dplyr::select(!c("id_type")) %>%
  distinct()

DrugBank_drug_info$names <- toupper(DrugBank_drug_info$names)

valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  dplyr::rename("ct_type" = "type", "NCT_ID"  = "_id") %>% 
  left_join(DrugBank_drug_info, 
            by = c("Drug1_name" = "names"), 
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug1_DrugBank_id" = "DrugBank_id", "Drug1_type" = "type") %>% 
  left_join(DrugBank_drug_info , 
            by = c("Drug2_name" = "names"),
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug2_DrugBank_id" = "DrugBank_id", "Drug2_type" = "type") 


valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  filter(!(is.na(Drug1_DrugBank_id) | is.na(Drug2_DrugBank_id))) %>% 
  filter((Drug1_type == "small molecule") & (Drug2_type == "small molecule"))

valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "NCT_ID", "class_EffAdv")]


#####


# Add the ATC codes at level_1
# Using many-to-many mapping to map all possible ATC codes to a single drug
valid_drugCombs_cat <- valid_drugCombs_cat %>%
  left_join(DrugBank_drug_ATC %>%
              dplyr::select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug1_ATC_", .)),
            by = c("Drug1_DrugBank_id" = "Drug1_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  left_join(DrugBank_drug_ATC %>%
              dplyr::select(code_1, DrugBank_drug_ID) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste0("Drug2_ATC_", .)),
            by = c("Drug2_DrugBank_id" = "Drug2_ATC_DrugBank_drug_ID"),
            relationship = "many-to-many") %>%
  distinct()


# Remove those with missing ATC
valid_drugCombs_cat <- valid_drugCombs_cat[!( is.na(valid_drugCombs_cat$Drug1_ATC_code_1) |
                                                is.na(valid_drugCombs_cat$Drug2_ATC_code_1) ), ]


# Filter to keep only those within the training framework
valid_drugCombs_cat <- valid_drugCombs_cat %>% filter((Drug1_ATC_code_1 %in% possible_ATC) & (Drug2_ATC_code_1 %in% possible_ATC))

valid_drugCombs_cat <- valid_drugCombs_cat %>% dplyr::select(!c(Drug1_ATC_code_1, Drug2_ATC_code_1)) %>% distinct()

if(nrow(valid_drugCombs_cat) == 0){ stop("No drug combinations found") }


#####


# Read the drug combination used in the training 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")

# Check for overlapping drug combinations
remove_rows <- c()
for(i in 1:nrow(valid_drugCombs_cat)){
  drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id", drop = TRUE]
  drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id", drop = TRUE]
  
  tmp1 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug1 & train_drugCombs_cat$Drug2_DrugBank_id == drug2, ]
  tmp2 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug2 & train_drugCombs_cat$Drug2_DrugBank_id == drug1, ]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    remove_rows <- c(remove_rows, i)
  }
}

if(length(remove_rows) > 0){
  valid_drugCombs_cat <- valid_drugCombs_cat[-remove_rows,]
}

if(nrow(valid_drugCombs_cat) == 0){ stop("No drug combinations found") }


#####


# Get list all unique combinations

all_drugs <- sort(unique(c(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id)))

all_drug_combs <- matrix( NA, 
                          nrow = length(all_drugs), 
                          ncol = length(all_drugs), 
                          dimnames = list(all_drugs, all_drugs) )
all_drug_combs[upper.tri(all_drug_combs, diag = FALSE)] <- NA


if(nrow(valid_drugCombs_cat) > 0){
  for(i in 1:nrow(valid_drugCombs_cat)){
    
    drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id", drop = TRUE]
    drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id", drop = TRUE]
    comb_info <- paste( unlist(valid_drugCombs_cat[i, c("NCT_ID", "class_EffAdv"), drop = TRUE], use.names = FALSE) , collapse = "___")
    comb_info <- paste0("[", comb_info, "]")
    
    if(!is.na(all_drug_combs[drug1, drug2])){
      all_drug_combs[drug1, drug2] <- paste(na.exclude(c(all_drug_combs[drug1, drug2], comb_info)), collapse = ";")
    }else{
      all_drug_combs[drug2, drug1] <- paste(na.exclude(c(all_drug_combs[drug2, drug1], comb_info)), collapse = ";")
    }
  }
}

all_drug_combs <- as.data.frame(all_drug_combs) %>% 
  rownames_to_column("Drug1_DrugBank_id") %>% 
  pivot_longer(-Drug1_DrugBank_id, 
               names_to = "Drug2_DrugBank_id", 
               values_to = "drugCombs") %>% 
  filter(!is.na(drugCombs))


all_drug_combs <- all_drug_combs %>% 
  separate_rows(drugCombs, sep = "\\];\\[") %>% 
  separate(drugCombs, into = c("NCT_ID", "class_EffAdv"), sep = "___") 

all_drug_combs$NCT_ID <- str_replace(all_drug_combs$NCT_ID, "^\\[", "")
all_drug_combs$class_EffAdv <- str_replace(all_drug_combs$class_EffAdv, "\\]$", "")

all_drug_combs$comb_name <- paste(all_drug_combs$Drug1_DrugBank_id, all_drug_combs$Drug2_DrugBank_id, sep = "_")

valid_drugCombs_cat <- all_drug_combs %>% group_by(comb_name) %>%  summarise(across(everything(), ~ toString(unique(.x))))

valid_drugCombs_cat[grepl("Adv, Eff|Eff, Adv", valid_drugCombs_cat$class_EffAdv), ]$class_EffAdv <- "Adv" # If any combination gets both categories, set it as adverse


#####


# Add the generic names of the drug combinations
valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  left_join(DrugBank_data$drugs$general_information %>% 
              dplyr::select(c("primary_key", "name")), 
            by = c("Drug1_DrugBank_id" = "primary_key"), 
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug1_name" = "name") %>% 
  left_join(DrugBank_data$drugs$general_information %>% 
              dplyr::select(c("primary_key", "name")) , 
            by = c("Drug2_DrugBank_id" = "primary_key"),
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug2_name" = "name") 


#####


# Read the network
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
input_network_nodes <- suppressMessages(select(org.Hs.eg.db, 
                                               keys = V(input_network)$name, 
                                               columns = "SYMBOL", 
                                               keytype = "ENSEMBL"))
input_network_nodes <- input_network_nodes$SYMBOL
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))


# Read drug target interactions
drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_associations.rds")


# Keep drug targets that are included in the input network
drug_target_ixn <- drug_target_ixn %>% 
  filter(ensembl_gene_id %in% V(input_network)$name) 


# Aggregate the drug targets as a single string
drugTarget_list <- drug_target_ixn %>%
  group_by(drugbank_drug_id) %>%
  summarise(drugTarget_ensembl_id = paste(ensembl_gene_id, collapse = ","))


# Merge the drug targets information with the drug combinations data
cat("\nExtracting targets of the drug combinations\n")

valid_drugCombs_cat <- valid_drugCombs_cat %>%
  left_join(drugTarget_list, by = c("Drug1_DrugBank_id" = "drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_1 = drugTarget_ensembl_id) %>%
  left_join(drugTarget_list, by = c("Drug2_DrugBank_id" = "drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_2 = drugTarget_ensembl_id) %>%
  filter(!is.na(drugTarget_ensembl_id_1), !is.na(drugTarget_ensembl_id_2)) %>%   
  dplyr::rowwise() %>%
  dplyr::mutate(
    drugTarget_ensembl_id =  list(unique(unlist(strsplit(na.exclude(c(drugTarget_ensembl_id_1, drugTarget_ensembl_id_2)), ",")))),
    drugTarget_count = length(unlist(drugTarget_ensembl_id))
  ) %>%
  dplyr::select(-drugTarget_ensembl_id_1, -drugTarget_ensembl_id_2) 


valid_drugCombs_cat <- valid_drugCombs_cat[valid_drugCombs_cat$drugTarget_count > 1, ]

# Simplify the target column as comma separated string
valid_drugCombs_cat$drugTarget_ensembl_id <- sapply(valid_drugCombs_cat$drugTarget_ensembl_id, function(x) paste(x, collapse = ","))


# Convert Ensembl IDs of the targets to gene symbols
valid_drugCombs_cat$drugTarget_geneSymbol <- sapply(valid_drugCombs_cat$drugTarget_ensembl_id, convert_ensembl_to_symbol, USE.NAMES = FALSE)


#####


# Extract known gene/protein interactions from various databases 
ixns_PS <- import_omnipath_interactions(resources = c("PhosphoSite","phosphoELM","PhosphoNetworks","PhosphoSite_MIMP",
                                                      "PhosphoSite_ProtMapper","phosphoELM_MIMP","PhosphoPoint",
                                                      "PhosphoSite_KEA","PhosphoSite_noref"))
ixns_SIGNOR <- import_omnipath_interactions(resources = c("SIGNOR","SIGNOR_ProtMapper","SIGNOR_CollecTRI"))
ixns_NPA <- import_omnipath_interactions(resources = c("NetPath"))
ixns_RIs <- import_transcriptional_interactions()


# Filter interactions to keep only nodes in the input network
ixns_PS <- ixns_PS[ixns_PS$source_genesymbol %in% input_network_nodes & ixns_PS$target_genesymbol %in% input_network_nodes, ]
ixns_SIGNOR <- ixns_SIGNOR[ixns_SIGNOR$source_genesymbol %in% input_network_nodes & ixns_SIGNOR$target_genesymbol %in% input_network_nodes, ]
ixns_NPA <- ixns_NPA[ixns_NPA$source_genesymbol %in% input_network_nodes & ixns_NPA$target_genesymbol %in% input_network_nodes, ]
ixns_RIs <- ixns_RIs[ixns_RIs$source_genesymbol %in% input_network_nodes & ixns_RIs$target_genesymbol %in% input_network_nodes, ]


# Extend the drug targets based on protein phosphorylation
cat("\nExtending drug targets based on protein phosphorylation information\n")
ps_result_list <- lapply(valid_drugCombs_cat$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_PS))
valid_drugCombs_cat$ext_PS_targets <- sapply(ps_result_list, function(x) x$ext_targets)
valid_drugCombs_cat$ext_PS_tar_cnt <- sapply(ps_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from SIGNOR database
cat("\nExtending drug targets based on interactions from SIGNOR database\n")
signor_result_list <- lapply(valid_drugCombs_cat$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_SIGNOR))
valid_drugCombs_cat$ext_SIGNOR_targets <- sapply(signor_result_list, function(x) x$ext_targets)
valid_drugCombs_cat$ext_SIGNOR_tar_cnt <- sapply(signor_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from NetPath database
cat("\nExtending drug targets based on interactions from NetPath database\n")
npa_result_list <- lapply(valid_drugCombs_cat$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_NPA))
valid_drugCombs_cat$ext_NPA_targets <- sapply(npa_result_list, function(x) x$ext_targets)
valid_drugCombs_cat$ext_NPA_tar_cnt <- sapply(npa_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on regulatory network (TF-target interactions)
cat("\nExtending drug targets based on regulatory network (TF-target interactions)\n")
reg_result_list <- lapply(valid_drugCombs_cat$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_RIs))
valid_drugCombs_cat$ext_RI_targets <- sapply(reg_result_list, function(x) x$ext_targets)
valid_drugCombs_cat$ext_RI_tar_cnt <- sapply(reg_result_list, function(x) x$ext_tar_cnt)


# Read the list of KEGG pathways linked to Cancer hallmarks
if(!dir.exists("Databases/CHG/")){dir.create("Databases/CHG/", recursive = TRUE)}
if(!file.exists("Databases/CHG/Table_1.xls")){
  download.file(url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7013921/bin/Table_1.xls",
                destfile = "Databases/CHG/Table_1.xls", method = "wget")
}

CHG_pathways <- read_excel("Databases/CHG/Table_1.xls", skip = 2)
CHG_pathways <- CHG_pathways$Pathway
CHG_pathways <- sort(unique(unlist(str_split(CHG_pathways, "\n"))))
CHG_pathways <- CHG_pathways[CHG_pathways != ""]


# Download protein interactions from KEGG and filter to keep proteins in the network
ixns_KEGG <- download_KEGG_ixns(pathway_ids = CHG_pathways)
ixns_KEGG <- lapply(ixns_KEGG, function(x){x[x$genesymbol_source %in% input_network_nodes & x$genesymbol_target %in% input_network_nodes,]})


# Extend the drug targets using interactions from KEGG 
cat("\nExtending drug targets using KEGG\n")
kegg_result_list <- lapply(valid_drugCombs_cat$drugTarget_geneSymbol, function(x) extend_targets_kegg_list(x, ixns_KEGG))
valid_drugCombs_cat$ext_KEGG_targets <- sapply(kegg_result_list, function(x) x$ext_kegg_targets)
valid_drugCombs_cat$ext_KEGG_tar_cnt <- sapply(kegg_result_list, function(x) x$ext_kegg_tar_cnt)


#####


if(!dir.exists("InputFiles/Validation_data_6b/")){dir.create("InputFiles/Validation_data_6b/", recursive = TRUE)}
saveRDS(valid_drugCombs_cat, file = paste0("InputFiles/Validation_data_6b/drugCombs_validation6b_", disease, ".rds"))

cat(paste0("\nNumber of drug combinations: ", nrow(valid_drugCombs_cat), "\n"))


#####


print(warnings())