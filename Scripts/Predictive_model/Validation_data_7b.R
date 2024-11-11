set.seed(5081)



# Generate validation data 7b
# Notes:
# (a) Based on the paper: A pan-cancer screen identifies drug combination benefit in cancer cell lines at the individual and population level



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
DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_ATC <- DrugBank_data$drugs$atc_codes
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


# Read the list of cell lines used in training dataset
training_cellLines <- read.csv("OutputFiles/Tables/FIMM_used_cell_lines.csv", header = TRUE)

query_tissue <- switch( disease,
                        "BreastCancer" = "breast",
                        "KidneyCancer" = "kidney",
                        "LungCancer" = "lung",
                        "OvaryCancer" = "ovary",
                        "ProstateCancer" = "prostate",
                        "SkinCancer" = "skin" )

training_cellLines <- training_cellLines %>% filter(tissue_name %in% query_tissue)
training_cellLines <- unique(training_cellLines$cell_line_name)
rm(query_tissue)


#####


# Download the data
if(!dir.exists("Databases/Vis_2024/")){dir.create("Databases/Vis_2024", recursive = TRUE)}

if(!file.exists("Databases/Vis_2024/Table_S3.xlsx")){
  download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2666379124004087-mmc4.xlsx",
                destfile = "Databases/Vis_2024/Table_S3.xlsx", method = "wget")
}


# Read the data
ECB_data <- openxlsx::read.xlsx(xlsxFile = "Databases/Vis_2024/Table_S3.xlsx", sheet = "Table3")

query_tissue <- switch( disease,
                        "BreastCancer" = "Breast",
                        "KidneyCancer" = "Kidney",
                        "LungCancer" = "Lung",
                        "OvaryCancer" = "Ovary",
                        # "ProstateCancer" = "",
                        "SkinCancer" = "Skin" )

ECB_data <- ECB_data %>% filter(tissue %in% query_tissue)


# Check the number of endpoints per combination + cell line group
tmp1 <- ECB_data %>% group_by(LIBRARY_ID, ANCHOR_ID, CELL_LINE_NAME) %>% summarise(ep_count = length(unique(Endpoint)))
print( paste0( "Number of endpoints per combination per cell line:: Min = ", min(tmp1$ep_count), " Max = ", max(tmp1$ep_count), " Avg = ", mean(tmp1$ep_count) ) )

# Check the number of cell lines per combination 
tmp1 <- ECB_data %>% group_by(LIBRARY_ID, ANCHOR_ID) %>% summarise(cellline_count = length(unique(CELL_LINE_NAME)))
print( paste0( "Number of cell lines per combination :: Min = ", min(tmp1$cellline_count), " Max = ", max(tmp1$cellline_count), " Avg = ", round(mean(tmp1$cellline_count), 0) ) )


# Keep only cell lines used in training
print("Number of cell lines in training:")
print( table(unique(ECB_data$CELL_LINE_NAME) %in% training_cellLines))

ECB_data <- ECB_data %>% filter(CELL_LINE_NAME %in% training_cellLines)


# Check the number of endpoints per combination + cell line group
tmp1 <- ECB_data %>% group_by(LIBRARY_ID, ANCHOR_ID, CELL_LINE_NAME) %>% summarise(ep_count = length(unique(Endpoint)))
print( paste0( "Number of endpoints per combination per cell line:: Min = ", min(tmp1$ep_count), " Max = ", max(tmp1$ep_count), " Avg = ", mean(tmp1$ep_count) ) )

# Check the number of cell lines per combination 
tmp1 <- ECB_data %>% group_by(LIBRARY_ID, ANCHOR_ID) %>% summarise(cellline_count = length(unique(CELL_LINE_NAME)))
print( paste0( "Number of cell lines per combination :: Min = ", min(tmp1$cellline_count), " Max = ", max(tmp1$cellline_count), " Avg = ", round(mean(tmp1$cellline_count), 0) ) )


# Merge across cell line with each combination to get the average ECB
# Higher value indicates higher changes of being effective across all cell lines

ECB_data <- ECB_data %>% 
  dplyr::select(c(LIBRARY_NAME, ANCHOR_NAME, CELL_LINE_NAME, ECB)) %>% 
  group_by(LIBRARY_NAME, ANCHOR_NAME) %>%
  summarise(mean_ECB = round(mean(ECB), 2))

ECB_data$ANCHOR_NAME <- toupper(ECB_data$ANCHOR_NAME)
ECB_data$LIBRARY_NAME <- toupper(ECB_data$LIBRARY_NAME)


#####


# Use the drug names to map to the Drugbank drug IDs

# DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")

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

valid_drugCombs_cat <- ECB_data %>% 
  left_join(DrugBank_drug_info, 
            by = c("LIBRARY_NAME" = "names"), 
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug1_DrugBank_id" = "DrugBank_id", "Drug1_type" = "type", "Drug1_name" = "LIBRARY_NAME") %>% 
  left_join(DrugBank_drug_info , 
            by = c("ANCHOR_NAME" = "names"), 
            relationship = "many-to-many") %>% 
  dplyr::rename("Drug2_DrugBank_id" = "DrugBank_id", "Drug2_type" = "type", "Drug2_name" = "ANCHOR_NAME") 


valid_drugCombs_cat <- valid_drugCombs_cat %>% 
  filter(!(is.na(Drug1_DrugBank_id) | is.na(Drug2_DrugBank_id))) %>% 
  filter((Drug1_type == "small molecule") & (Drug2_type == "small molecule"))


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
    comb_info <- paste( unlist(valid_drugCombs_cat[i, c("mean_ECB"), drop = TRUE], use.names = FALSE) , collapse = "___")
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
  separate(drugCombs, into = c("mean_ECB"), sep = "___") 

all_drug_combs$mean_ECB <- as.numeric(str_replace_all(all_drug_combs$mean_ECB, pattern = "^\\[|\\]$", ""))


# Some combinations have different mean ECB since based on their position as anchor 
# or library, they were tested on different number of cell lines
# Use mean here
all_drug_combs <- all_drug_combs %>% 
  group_by(Drug1_DrugBank_id, Drug2_DrugBank_id) %>% 
  summarise(mean_ECB = mean((mean_ECB)))

all_drug_combs$comb_name <- paste(all_drug_combs$Drug1_DrugBank_id, all_drug_combs$Drug2_DrugBank_id, sep = "_")

valid_drugCombs_cat <- all_drug_combs
rm(all_drug_combs)


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


# Read the DDI data
DrugBank_ddi <- DrugBank_data$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")


# Keep only DDI involving the drugs in the validation dataset
DrugBank_ddi <- DrugBank_ddi %>% filter(Drug1_DrugBank_id %in% unique(c(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id)) & 
                                          Drug2_DrugBank_id %in% unique(c(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id)))



# Extract DDI with toxicity
DrugBank_ddi <- DrugBank_ddi[grep("risk or severity of .+toxicity can be increased|risk or severity of liver damage can be increased|risk or severity of adverse effects can be increased |increase the .+toxic activities", DrugBank_ddi$description, ignore.case = TRUE), ]


DrugBank_ddi$DDI_type <- gsub(pattern = ".*(risk or severity of .+ can be increased).*|.*(increase the .+toxic activities).*", 
                              replacement = "\\1\\2", 
                              x = DrugBank_ddi$description, 
                              ignore.case = TRUE)

sort(table(DrugBank_ddi$DDI_type ), decreasing = TRUE)

# Merge the data
valid_drugCombs_cat <- valid_drugCombs_cat %>%
  left_join(DrugBank_ddi %>% dplyr::select(c(Drug1_DrugBank_id, Drug2_DrugBank_id, DDI_type)),
            by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"))


# Assign class
valid_drugCombs_cat$class_EffAdv <- NA
valid_drugCombs_cat[(valid_drugCombs_cat$mean_ECB > 0) & is.na(valid_drugCombs_cat$DDI_type), "class_EffAdv"] <- "Eff"
valid_drugCombs_cat[(valid_drugCombs_cat$mean_ECB > 0) & !is.na(valid_drugCombs_cat$DDI_type), "class_EffAdv"] <- "Adv"
valid_drugCombs_cat <- valid_drugCombs_cat[!is.na(valid_drugCombs_cat$class_EffAdv), ]


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


if(!dir.exists("InputFiles/Validation_data_7b/")){dir.create("InputFiles/Validation_data_7b/", recursive = TRUE)}
saveRDS(valid_drugCombs_cat, file = paste0("InputFiles/Validation_data_7b/drugCombs_validation7b_", disease, ".rds"))

cat(paste0("\nNumber of drug combinations: ", nrow(valid_drugCombs_cat), "\n"))


#####


print(warnings())