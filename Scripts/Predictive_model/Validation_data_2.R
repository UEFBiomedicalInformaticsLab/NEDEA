set.seed(5081)



# Generate validation data 2 (with drugs with ATC codes in the training)
# Notes:
# (a) Includes only adverse drug combinations



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(igraph)
library(org.Hs.eg.db)
library(OmnipathR)
library(readxl)
source("Scripts/Functions/Functions_drug_target.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


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
               names_to = "ATC2", 
               values_to = "Count")

ATC_count_mat <- ATC_count_mat %>% filter(!is.na(Count)) %>% filter(Count > 0)

possible_ATC_pairs <- unique(paste(ATC_count_mat$ATC1, ATC_count_mat$ATC2, sep = "_"))


#####


# Download list of licensed anti-cancer drugs
if(!dir.exists("Databases/Cancer_Drug_Database/")){dir.create("Databases/Cancer_Drug_Database/", recursive = TRUE)}
if(!file.exists("Databases/Cancer_Drug_Database/cancerdrugsdb.txt")){
  download.file(url = "https://sciencedata.anticancerfund.org/pages//cancerdrugsdb.txt",
                destfile = "Databases/Cancer_Drug_Database/cancerdrugsdb.txt", method = "wget")
}

cancer_drug_list <- read.table("Databases/Cancer_Drug_Database/cancerdrugsdb.txt", sep = "\t", 
                               fill = TRUE, header = TRUE, strip.white = TRUE, check.names = FALSE)
cancer_drug_list <- cancer_drug_list[, c("Product", "DrugBank ID", "Indications")]
cancer_drug_list <- cancer_drug_list %>% mutate(DrugBank_drug_id = gsub("<.*>(DB\\d{5})</.*", "\\1", `DrugBank ID`))


# Extract the cancer specific drug list
grep_pattern <- switch (disease,
                        "BreastCancer" = "breast cancer|breast carcinoma",
                        "KidneyCancer" = "renal cell cancer|renal cell carcinoma",
                        "LungCancer" = "lung cancer|lung carcinoma",
                        "OvaryCancer" = "ovarian cancer|ovarian carcinoma|ovarian epithelial cancer",
                        "ProstateCancer" = "prostate Cancer|carcinoma of the prostate",
                        "SkinCancer" = "carcinoma of the skin|melanoma"
)

cancer_drug_list <- cancer_drug_list[grep(pattern = grep_pattern, x = cancer_drug_list$Indications, ignore.case = TRUE), ]


#####


# Read the drug type information
DrugBank_drug_type <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_type <- DrugBank_drug_type$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


# Read the DDI data
DrugBank_ddi <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_ddi <- DrugBank_ddi$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% DrugBank_drug_type$primary_key & DrugBank_ddi$Drug2_DrugBank_id %in% DrugBank_drug_type$primary_key, ] # retain only small molecular drugs


# Extract DDI with toxicity
DrugBank_ddi <- DrugBank_ddi[grep("risk or severity of .+toxicity can be increased|risk or severity of liver damage can be increased|risk or severity of adverse effects can be increased |increase the .+toxic activities", DrugBank_ddi$description, ignore.case = TRUE), ]


DrugBank_ddi$DDI_type <- gsub(pattern = ".*(risk or severity of .+ can be increased).*|.*(increase the .+toxic activities).*", 
                              replacement = "\\1\\2", 
                              x = DrugBank_ddi$description, 
                              ignore.case = TRUE)

sort(table(DrugBank_ddi$DDI_type ), decreasing = TRUE)


# Keep only combinations involving licensed anti-cancer drugs
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% cancer_drug_list$DrugBank_drug_id | DrugBank_ddi$Drug2_DrugBank_id %in% cancer_drug_list$DrugBank_drug_id, ]


#####


# Add the ATC codes at level_1
# Using many-to-many mapping to map all possible ATC codes to a single drug
DrugBank_ddi <- DrugBank_ddi %>%
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
DrugBank_ddi <- DrugBank_ddi[!(is.na(DrugBank_ddi$Drug1_ATC_code_1) | is.na(DrugBank_ddi$Drug2_ATC_code_1)), ]

# Filter to keep only those within the training framework
DrugBank_ddi <- DrugBank_ddi[paste(DrugBank_ddi$Drug1_ATC_code_1, DrugBank_ddi$Drug2_ATC_code_1, sep = "_") %in% possible_ATC_pairs | 
                               paste(DrugBank_ddi$Drug2_ATC_code_1, DrugBank_ddi$Drug1_ATC_code_1, sep = "_") %in% possible_ATC_pairs, ]

DrugBank_ddi <- DrugBank_ddi %>% dplyr::select(!c(Drug1_ATC_code_1, Drug2_ATC_code_1)) %>% distinct()


#####


# Extract unique list of DDIs

DrugBank_ddi$keep <- NA
for(i in 1:nrow(DrugBank_ddi)){
  if(is.na(DrugBank_ddi[i,"keep"])){
    DrugBank_ddi[i,"keep"] <- TRUE
    
    drug1 <- DrugBank_ddi[i, "Drug1_DrugBank_id", drop = TRUE]
    drug2 <- DrugBank_ddi[i, "Drug2_DrugBank_id", drop = TRUE]
    
    DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% drug2 & DrugBank_ddi$Drug2_DrugBank_id %in% drug1, "keep"] <- FALSE
  }
}

valid_drugCombs_cat <- DrugBank_ddi[DrugBank_ddi$keep == "TRUE", ]


# Assign labels
valid_drugCombs_cat$class_EffAdv <- "Adv"
valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_EffAdv", "description", "DDI_type")]
colnames(valid_drugCombs_cat)[4] <- "DDI_description"


######


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


valid_drugCombs_cat$comb_name <- paste(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id, sep = "_")


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
    drugTarget_ensembl_id =  list(unique(unlist(strsplit(c(drugTarget_ensembl_id_1, drugTarget_ensembl_id_2), ",")))),
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


# Add annotations to the combinations
valid_drugCombs_cat$Drug1_indications <- cancer_drug_list$Indications[match(valid_drugCombs_cat$Drug1_DrugBank_id, cancer_drug_list$DrugBank_drug_id)]
valid_drugCombs_cat$Drug2_indications <- cancer_drug_list$Indications[match(valid_drugCombs_cat$Drug2_DrugBank_id, cancer_drug_list$DrugBank_drug_id)]


#####


# Remove DDI types with 10 or less than 10 drug combinations
keep_DDI_type <- names(which(table(valid_drugCombs_cat$DDI_type, useNA = "ifany") > 10))
valid_drugCombs_cat <- valid_drugCombs_cat[valid_drugCombs_cat$DDI_type %in% keep_DDI_type, ]


#####

if(!dir.exists("InputFiles/Validation_data_2/")){dir.create("InputFiles/Validation_data_2/", recursive = TRUE)}
saveRDS(valid_drugCombs_cat, file = paste0("InputFiles/Validation_data_2/drugCombs_validation2_", disease, ".rds"))

cat(paste0("\nNumber of drug combinations: ", nrow(valid_drugCombs_cat), "\n"))


#####


print(warnings())