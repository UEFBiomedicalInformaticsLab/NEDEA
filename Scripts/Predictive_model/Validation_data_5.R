
# Load libraries
library(tidyverse)





disease <- "KidneyCancer"






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

# Remove those with missing ATC 
train_drugCombs_cat <- train_drugCombs_cat[!(is.na(train_drugCombs_cat$Drug1_ATC_code_1) | is.na(train_drugCombs_cat$Drug2_ATC_code_1)), ]


# Get list of ATC pairs 
all_atc <- sort(unique(c(train_drugCombs_cat$Drug1_ATC_code_1, train_drugCombs_cat$Drug2_ATC_code_1)))

ATC_count_mat <- matrix(0, 
                        nrow = length(all_atc), 
                        ncol = length(all_atc), 
                        dimnames = list(all_atc, all_atc)
)
ATC_count_mat[upper.tri(ATC_count_mat, diag = FALSE)] <- NA



if(nrow(train_drugCombs_cat) > 1){
  for(i in 1:nrow(train_drugCombs_cat)){
    
    atc_1 <- train_drugCombs_cat[i, "Drug1_ATC_code_1"]
    atc_2 <- train_drugCombs_cat[i, "Drug2_ATC_code_1"]
    
    
    if(!is.na(ATC_count_mat[atc_1, atc_2])){
      ATC_count_mat[atc_1, atc_2] <- ATC_count_mat[atc_1, atc_2] + 1
    }else{
      ATC_count_mat[atc_2, atc_1] <- ATC_count_mat[atc_2, atc_1] + 1
      
    }
  }
}

ATC_count_mat <- as.data.frame(ATC_count_mat) %>% 
  rownames_to_column("ATC1") %>% 
  pivot_longer(-ATC1, 
               names_to= "ATC2", 
               values_to = "Count")

ATC_count_mat <- ATC_count_mat %>% filter(!is.na(Count)) %>% filter(Count > 0)

possible_ATC_pairs <- unique(paste(ATC_count_mat$ATC1, ATC_count_mat$ATC2, sep = "_"))


#####


# Download and read the data
if(!dir.exists("Databases/Shi_2024/"))dir.create("Databases/Shi_2024/", recursive = TRUE)
if(!file.exists("Databases/Shi_2024/efficacy_df.csv")){
  stop("Manually download the files 'Clinical trials efficacy results (csv)' and 'Clinical trials safety results (csv)' from 'https://doi.org/10.6084/m9.figshare.c.6860254.v1' ")
}

CT_efficacy <- read.csv("Databases/Shi_2024/efficacy_df.csv", 
                        header = TRUE, 
                        row.names = 1)
row.names(CT_efficacy) <- NULL


CT_safety <- read.csv("Databases/Shi_2024/safety_df.csv", 
                      header = TRUE, 
                      row.names = 1)
row.names(CT_safety) <- NULL


#####


# Process the efficacy data 
switch(disease, 
       "BreastCancer" = { pattern = "Breast Neoplasms" },
       "KidneyCancer" = { pattern = "Carcinoma, Renal Cell|Kidney Neoplasms" },
       "LungCancer" = { pattern = "Carcinoma, Non-Small-Cell Lung|Lung Neoplasms|Small Cell Lung Carcinoma" },
       "OvaryCancer" = { pattern = "Ovarian Neoplasms|Carcinoma, Ovarian" },
       "ProstateCancer" = { pattern = "Prostatic Neoplasms" },
       "SkinCancer" = { pattern = "Skin Neoplasms|Melanoma" },
       )

# CT_efficacy <- CT_efficacy[grep("Breast Neoplasms", CT_efficacy$condition_mesh, ignore.case = TRUE), ]
# CT_efficacy <- CT_efficacy[grep("Carcinoma, Renal Cell|Kidney Neoplasms", CT_efficacy$condition_mesh, ignore.case = TRUE), ]
# CT_efficacy <- CT_efficacy[grep("Carcinoma, Non-Small-Cell Lung|Lung Neoplasms|Small Cell Lung Carcinoma", CT_efficacy$condition_mesh, ignore.case = TRUE), ]
# CT_efficacy <- CT_efficacy[grep("Ovarian Neoplasms|Carcinoma, Ovarian", CT_efficacy$condition_mesh, ignore.case = TRUE), ]
# CT_efficacy <- CT_efficacy[grep("Prostatic Neoplasms", CT_efficacy$condition_mesh, ignore.case = TRUE), ]
# CT_efficacy <- CT_efficacy[grep("Skin Neoplasms|Melanoma", CT_efficacy$condition_mesh, ignore.case = TRUE), ]

CT_efficacy <- CT_efficacy[grep(pattern, CT_efficacy$condition_mesh, ignore.case = TRUE), ]

CT_efficacy <- CT_efficacy[grep("\\['Drug', 'Drug'\\]", CT_efficacy$intervention_type, ignore.case = TRUE), ] # Keep only trial arms involving two drugs

# CT_efficacy <- CT_efficacy[grep("\\['placebo'\\]", CT_efficacy$comparator, ignore.case = TRUE), ] # Keep only comparisons with placebo
# CT_efficacy <- CT_efficacy[grep("primary", CT_efficacy$outcome_type, ignore.case = TRUE), ] # Check only primary outcomes
# CT_efficacy <- CT_efficacy[CT_efficacy$label %in% "positive", ] # Check only positive outcomes i.e., this arm performed better than the comparator group
CT_efficacy <- CT_efficacy  %>% 
  select(c("NCT_ID", "trial_phase", "condition_mesh", 
           "intervention", "comparator", "intervention_group")) %>%
  distinct()


####


# Process the safety data
CT_safety <- CT_safety[CT_safety$NCT_ID %in% CT_efficacy$NCT_ID, ]
CT_safety <- CT_safety[grep("placebo", CT_safety$arm_title, ignore.case = TRUE, invert = TRUE), ] # Keep only arms which is not the comparator (i.e., not the placebo arm)
CT_safety <- CT_safety[CT_safety$category %in% c("Total"), ] # Will check only the total number of patients who died or had serious AE
CT_safety$affected_risk_perc <- (CT_safety$affected/CT_safety$at_risk) * 100 # Calculate the % of patients died or had serious AE of all the patients under observation
CT_safety$is_adverse <- ifelse( (CT_safety$affected_risk_perc < 20), FALSE, TRUE) # Set criteria to define adverse combinations 

# If the arm is adverse based on either mortality or serious AE, it is marked as adverse
CT_safety <- CT_safety %>%
  group_by(NCT_ID, arm_title, matched_group) %>% 
  summarise("is_adverse" = any(is_adverse, na.rm = TRUE), 
            .groups = "keep")


#####


# Merge the efficacy and safety parts and ideity drug combinations for validation
valid_drugCombs_cat <- CT_efficacy %>%
  left_join(CT_safety, 
            by = c("NCT_ID" = "NCT_ID", "intervention_group" = "matched_group"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(is_adverse))

valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  separate(col = "intervention", 
           into = c("Drug1_name", "Drug2_name"), 
           sep = "\\', \\'")

valid_drugCombs_cat$Drug1_name <- gsub("^\\[\\'", "", valid_drugCombs_cat$Drug1_name)
valid_drugCombs_cat$Drug2_name <- gsub("\\'\\]$", "", valid_drugCombs_cat$Drug2_name)
valid_drugCombs_cat$class_EffAdv <- ifelse(valid_drugCombs_cat$is_adverse, "Adv", "Eff")


#####


# Manual modifications on drug names
valid_drugCombs_cat$Drug1_name <- gsub("Atezolizumab \\(MPDL3280A\\), an engineered anti-PDL1 antibody", "Atezolizumab", valid_drugCombs_cat$Drug1_name)
valid_drugCombs_cat$Drug2_name <- gsub("Atezolizumab \\(MPDL3280A\\), an engineered anti-PDL1 antibody", "Atezolizumab", valid_drugCombs_cat$Drug2_name)

valid_drugCombs_cat$Drug1_name <- gsub("AZD8931", "Sapitinib", valid_drugCombs_cat$Drug1_name)
valid_drugCombs_cat$Drug2_name <- gsub("AZD8931", "Sapitinib", valid_drugCombs_cat$Drug2_name)

valid_drugCombs_cat$Drug1_name <- toupper(valid_drugCombs_cat$Drug1_name)
valid_drugCombs_cat$Drug2_name <- toupper(valid_drugCombs_cat$Drug2_name)


#####


# Read the drug info from Drug Bank
DrugBank_drug_info <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_info <- DrugBank_drug_info$drugs$general_information
DrugBank_drug_info <- DrugBank_drug_info[, c("primary_key", "name", "type")]
colnames(DrugBank_drug_info)[1] <- "DrugBank_id"
DrugBank_drug_info$name <- toupper(DrugBank_drug_info$name)

valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  left_join(DrugBank_drug_info, 
            by = c("Drug1_name" = "name")) %>% 
  rename("Drug1_DrugBank_id" = "DrugBank_id", "Drug1_type" = "type") %>% 
  left_join(DrugBank_drug_info, 
            by = c("Drug2_name" = "name")) %>% 
  rename("Drug2_DrugBank_id" = "DrugBank_id", "Drug2_type" = "type")


valid_drugCombs_cat <- valid_drugCombs_cat %>% filter((Drug1_type == "small molecule") | (Drug2_type == "small molecule"))
valid_drugCombs_cat$condition_mesh <- gsub("\\[|\\]|\\'", "", valid_drugCombs_cat$condition_mesh)

valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "Drug1_name", "Drug2_name", 
                                               "NCT_ID", "condition_mesh", "trial_phase", "class_EffAdv")]


#####


# Get list all unique combinations

all_drugs <- sort(unique(c(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id)))

all_drug_combs <- matrix(NA, 
                         nrow = length(all_drugs), 
                         ncol = length(all_drugs), 
                         dimnames = list(all_drugs, all_drugs)
)
all_drug_combs[upper.tri(all_drug_combs, diag = FALSE)] <- NA


if(nrow(valid_drugCombs_cat) > 0){
  for(i in 1:nrow(valid_drugCombs_cat)){
    
    drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id", drop = TRUE]
    drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id", drop = TRUE]
    comb_info <- paste( unlist(valid_drugCombs_cat[i, , drop = TRUE], use.names = FALSE) , collapse = "___")
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


all_drug_combs <- all_drug_combs %>% separate_rows(drugCombs, sep = "\\];\\[") %>% 
  separate(drugCombs, into = c("Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                               "Drug1_name", "Drug2_name", 
                               "NCT_ID", "condition_mesh", "trial_phase", "class_EffAdv"), 
           sep = "___") 

all_drug_combs$Drug1_DrugBank_id <- str_replace(all_drug_combs$Drug1_DrugBank_id, "^\\[", "")
all_drug_combs$class_EffAdv <- str_replace(all_drug_combs$class_EffAdv, "\\]$", "")

all_drug_combs$comb_name <- paste(all_drug_combs$Drug1_DrugBank_id, all_drug_combs$Drug2_DrugBank_id, sep = "_")


valid_drugCombs_cat <- all_drug_combs %>% group_by(comb_name) %>%  summarise(across(everything(), ~ toString(unique(.x))))

valid_drugCombs_cat$class_EffAdv <- ifelse(grepl("Adv", valid_drugCombs_cat$class_EffAdv), "Adv", "Eff")


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
valid_drugCombs_cat <- valid_drugCombs_cat[paste(valid_drugCombs_cat$Drug1_ATC_code_1, valid_drugCombs_cat$Drug2_ATC_code_1, sep = "_") %in% possible_ATC_pairs | 
                                             paste(valid_drugCombs_cat$Drug2_ATC_code_1, valid_drugCombs_cat$Drug1_ATC_code_1, sep = "_") %in% possible_ATC_pairs, ]

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




print(valid_drugCombs_cat)












print(warnings())