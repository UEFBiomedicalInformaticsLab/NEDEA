set.seed(5081)



# Generate validation data 4
# Notes:
# (a) 



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


# Read the drug combinations 
valid_drugCombs_cat <- read.csv("Databases/Manual_curation/Curated_drug_combination_dataset.csv", header = TRUE)
valid_drugCombs_cat <- valid_drugCombs_cat[valid_drugCombs_cat$Disease == disease, ]


#####


# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
DrugBank_ddi <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", paste0("ADR_", disease))]
DrugBank_ddi$comb_name <- paste(DrugBank_ddi$Drug1_DrugBank_id, DrugBank_ddi$Drug2_DrugBank_id, sep = "_")


# Assign labels
# For two drug combinations, if known ADR present, classify as adverse else effective
# For three drug combinations, if known ADR present for at least one pair, classify as adverse else effective

valid_drugCombs_cat$class_EffAdv <- NA
for(i in 1:nrow(valid_drugCombs_cat)){
  
  drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id"]
  drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id"]
  drug3 <- valid_drugCombs_cat[i, "Drug3_DrugBank_id"]
  
  if(is.na(drug3)){
    
    tmp1 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug1 & DrugBank_ddi$Drug2_DrugBank_id == drug2, paste0("ADR_", disease), drop = TRUE]
    if(length(tmp1) == 0){tmp1 <- "unknown"}
    valid_drugCombs_cat[i, "class_EffAdv"] <- ifelse(tmp1 == "adr_positive", "Adv", "Eff")
    rm(tmp1)
    
  }else{
    
    tmp1 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug1 & DrugBank_ddi$Drug2_DrugBank_id == drug2, paste0("ADR_", disease), drop = TRUE]
    if(length(tmp1) == 0){tmp1 <- "unknown"}
    
    tmp2 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug2 & DrugBank_ddi$Drug2_DrugBank_id == drug3, paste0("ADR_", disease), drop = TRUE]
    if(length(tmp2) == 0){tmp2 <- "unknown"}
    
    tmp3 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug3 & DrugBank_ddi$Drug2_DrugBank_id == drug1, paste0("ADR_", disease), drop = TRUE]
    if(length(tmp3) == 0){tmp3 <- "unknown"}
    
    valid_drugCombs_cat[i, "class_EffAdv"] <- ifelse(any(tmp1 == "adr_positive", tmp2 == "adr_positive", tmp3 == "adr_positive"), "Adv", "Eff")
    rm(list = c("tmp1", "tmp2", "tmp3"))
    
  }
}


######


# Read the drug combination used in the training 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")

# Check for overlapping drug combinations
remove_rows <- c()
for(i in 1:nrow(valid_drugCombs_cat)){
  drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id"]
  drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id"]
  
  tmp1 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug1 & train_drugCombs_cat$Drug2_DrugBank_id == drug2, ]
  tmp2 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug2 & train_drugCombs_cat$Drug2_DrugBank_id == drug1, ]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    remove_rows <- c(remove_rows, i)
  }
}

if(length(remove_rows) > 0){
  valid_drugCombs_cat <- valid_drugCombs_cat[-remove_rows,]
}


valid_drugCombs_cat$comb_name <- paste(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id, valid_drugCombs_cat$Drug3_DrugBank_id, sep = "_")
valid_drugCombs_cat$comb_name <- gsub("_NA$", "", valid_drugCombs_cat$comb_name)


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
  left_join(drugTarget_list, by = c("Drug3_DrugBank_id" = "drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_3 = drugTarget_ensembl_id) %>%
  filter(!is.na(drugTarget_ensembl_id_1), !is.na(drugTarget_ensembl_id_2)) %>%   
  dplyr::rowwise() %>%
  dplyr::mutate(
    drugTarget_ensembl_id =  list(unique(unlist(strsplit(na.exclude(c(drugTarget_ensembl_id_1, drugTarget_ensembl_id_2, drugTarget_ensembl_id_3)), ",")))),
    drugTarget_count = length(unlist(drugTarget_ensembl_id))
  ) %>%
  dplyr::select(-drugTarget_ensembl_id_1, -drugTarget_ensembl_id_2, -drugTarget_ensembl_id_2, -drugTarget_ensembl_id_3) 

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


if(!dir.exists("InputFiles/Validation_data_4/")){dir.create("InputFiles/Validation_data_4/", recursive = TRUE)}
saveRDS(valid_drugCombs_cat, file = paste0("InputFiles/Validation_data_4/drugCombs_validation4_", disease, ".rds"))

cat(paste0("\nNumber of drug combinations: ", nrow(valid_drugCombs_cat), "\n"))



print(warnings())