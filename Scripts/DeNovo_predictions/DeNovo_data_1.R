set.seed(5081)



# Script to generate drug combinations for de novo predictions



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


# Read the drug type information
DrugBank_drug_type <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_type <- DrugBank_drug_type$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


# Download list of licensed anti-cancer drugs
if(!dir.exists("Databases/Cancer_Drug_Database/")){dir.create("Databases/Cancer_Drug_Database/", recursive = TRUE)}
if(!file.exists("Databases/Cancer_Drug_Database/cancerdrugsdb.txt")){
  download.file(url = "https://sciencedata.anticancerfund.org/pages//cancerdrugsdb.txt",
                destfile = "Databases/Cancer_Drug_Database/cancerdrugsdb.txt", method = "wget")
}

cancer_drug_list <- read.table("Databases/Cancer_Drug_Database/cancerdrugsdb.txt", sep = "\t", 
                               fill = TRUE, header = TRUE, strip.white = TRUE, check.names = FALSE)

cancer_drug_list <- cancer_drug_list[, !colnames(cancer_drug_list) %in% c("")]
cancer_drug_list <- cancer_drug_list %>% mutate(DrugBank_drug_id = gsub("<.*>(DB\\d{5})</.*", "\\1", `DrugBank ID`))
cancer_drug_list <- cancer_drug_list[, c("DrugBank_drug_id", "EMA", "FDA", "WHO", "Year", "Generic", "Product", "Indications")] # 312
cancer_drug_list <- cancer_drug_list[cancer_drug_list$EMA == "Y" & cancer_drug_list$FDA == "Y" , ] # 197
cancer_drug_list <- cancer_drug_list[cancer_drug_list$Generic == "Y", ] # 48
cancer_drug_list <- cancer_drug_list[cancer_drug_list$DrugBank_drug_id %in% DrugBank_drug_type$primary_key, ] # retain only small molecular drugs


# Generate all possible combinations
denovo_drugCombs <- as.data.frame(t(combn(cancer_drug_list$DrugBank_drug_id, m = 2, simplify = TRUE)))
colnames(denovo_drugCombs) <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")
denovo_drugCombs <- denovo_drugCombs[denovo_drugCombs$Drug1_DrugBank_id != denovo_drugCombs$Drug2_DrugBank_id, ]


# Add annotations to the drug combinations
denovo_drugCombs <- denovo_drugCombs %>% 
  left_join(cancer_drug_list %>% rename_with(.cols = everything(), 
                                             .fn = ~ paste0("Drug1_", .)), 
            by = c("Drug1_DrugBank_id" = "Drug1_DrugBank_drug_id")) %>%
  left_join(cancer_drug_list %>% rename_with(.cols = everything(), 
                                             .fn = ~ paste0("Drug2_", .)), 
            by = c("Drug2_DrugBank_id" = "Drug2_DrugBank_drug_id"))


####


# Read the drug combination used in the training 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")

# Check for overlapping drug combinations
remove_rows <- c()
for(i in 1:nrow(denovo_drugCombs)){
  drug1 <- denovo_drugCombs[i, "Drug1_DrugBank_id", drop = TRUE]
  drug2 <- denovo_drugCombs[i, "Drug2_DrugBank_id", drop = TRUE]
  
  tmp1 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug1 & train_drugCombs_cat$Drug2_DrugBank_id == drug2, ]
  tmp2 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug2 & train_drugCombs_cat$Drug2_DrugBank_id == drug1, ]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    remove_rows <- c(remove_rows, i)
  }
}

if(length(remove_rows) > 0){
  denovo_drugCombs <- denovo_drugCombs[-remove_rows,]
}


#####


# Download the C-DCDB combinations
if(!dir.exists("Databases/CDCDB/")){dir.create("Databases/CDCDB/", recursive = TRUE)}
if(!file.exists("Databases/CDCDB/09.04.2024.zip")){
  cat("\n\nNOTE: Manually download the file '09.04.2024.zip' from the 'Downloads' page in https://icc.ise.bgu.ac.il/medical_ai/CDCDB/ and place it within 'Databases/CDCDB/'. \n\n")
}
if(!file.exists("Databases/CDCDB/web_preview.csv")){
  unzip(zipfile = "Databases/CDCDB/09.04.2024.zip", exdir = "Databases/CDCDB/")
}


# Read the drug combinations
CDCDB_drugCombs <- read.csv("Databases/CDCDB/web_preview.csv")
CDCDB_drugCombs <- CDCDB_drugCombs[CDCDB_drugCombs$source == "clinicaltrials.gov", c("drugs", "drugbank_identifiers", "source_id")] # Keep only those from clinical trials
CDCDB_drugCombs <- CDCDB_drugCombs %>%
  filter(str_count(drugbank_identifiers, ";") %in% c(1)) %>%
  filter(!str_detect(drugbank_identifiers, "NA|PLACEBO")) 

CDCDB_drugCombs <- separate(CDCDB_drugCombs,
                            col = "drugbank_identifiers", 
                            into = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"), 
                            sep = ";", fill = "right")


# Add info about condition for the trials
CDCDB_conditions <- read.csv("Databases/CDCDB/conditions_df.csv")
CDCDB_drugCombs$condition <- CDCDB_conditions$condition_downcase[match(CDCDB_drugCombs$source_id, CDCDB_conditions$nct_id)]


# Extract the cancer specific drug combinations that have been in clinical trials
grep_pattern <- switch (disease,
                        "BreastCancer" = "breast cancer|breast carcinoma",
                        "KidneyCancer" = "renal cell cancer|renal cell carcinoma",
                        "LungCancer" = "lung cancer|lung carcinoma",
                        "OvaryCancer" = "ovarian cancer|ovarian carcinoma|ovarian epithelial cancer|carcinoma, ovarian epithelial",
                        "ProstateCancer" = "prostate Cancer|prostate carcinoma|prostate adenocarcinoma",
                        "SkinCancer" = "skin cancer|melanoma"
)

CDCDB_drugCombs <- CDCDB_drugCombs[grep(pattern = grep_pattern, x = CDCDB_drugCombs$condition, ignore.case = TRUE), ]


#####


# Remove the new drug pairs which has already been used in clinical trials
remove_rows <- c()
for(i in 1:nrow(denovo_drugCombs)){
  drug1 <- denovo_drugCombs[i, "Drug1_DrugBank_id", drop = TRUE]
  drug2 <- denovo_drugCombs[i, "Drug2_DrugBank_id", drop = TRUE]
  
  tmp1 <- CDCDB_drugCombs[CDCDB_drugCombs$Drug1_DrugBank_id == drug1 & CDCDB_drugCombs$Drug2_DrugBank_id == drug2, ]
  tmp2 <- CDCDB_drugCombs[CDCDB_drugCombs$Drug1_DrugBank_id == drug2 & CDCDB_drugCombs$Drug2_DrugBank_id == drug1, ]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    remove_rows <- c(remove_rows, i)
  }
}

if(length(remove_rows) > 0){
  denovo_drugCombs <- denovo_drugCombs[-remove_rows,]
}
rm(remove_rows)


#####


# Read the DDI data
DrugBank_ddi <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_ddi <- DrugBank_ddi$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% DrugBank_drug_type$primary_key & DrugBank_ddi$Drug2_DrugBank_id %in% DrugBank_drug_type$primary_key, ] # retain only small molecular drugs


# Extract DDI with toxicity
DrugBank_ddi <- DrugBank_ddi[grep("risk or severity of .+toxicity can be increased|risk or severity of liver damage can be increased|risk or severity of adverse effects can be increased |increase the .+toxic activities", DrugBank_ddi$description, ignore.case = TRUE), ]


# Add annotation of known DDI
denovo_drugCombs$DDI_description <- NA
for(i in 1:nrow(denovo_drugCombs)){
  drug1 <- denovo_drugCombs[i, "Drug1_DrugBank_id", drop = TRUE]
  drug2 <- denovo_drugCombs[i, "Drug2_DrugBank_id", drop = TRUE]
  
  tmp1 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug1 & DrugBank_ddi$Drug2_DrugBank_id == drug2, "description"]
  tmp2 <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id == drug2 & DrugBank_ddi$Drug2_DrugBank_id == drug1, "description"]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    denovo_drugCombs[i, "DDI_description"] <- paste(unique(c(tmp1, tmp2)), collapse = "; ")
  }
}


# Remove drug combinations with known adverse DDI
denovo_drugCombs <- denovo_drugCombs[is.na(denovo_drugCombs$DDI_description), ]


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
denovo_drugCombs <- denovo_drugCombs %>%
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

denovo_drugCombs <- denovo_drugCombs[denovo_drugCombs$drugTarget_count > 1, ]

# Simplify the target column as comma separated string
denovo_drugCombs$drugTarget_ensembl_id <- sapply(denovo_drugCombs$drugTarget_ensembl_id, function(x) paste(x, collapse = ","))


# Convert Ensembl IDs of the targets to gene symbols
denovo_drugCombs$drugTarget_geneSymbol <- sapply(denovo_drugCombs$drugTarget_ensembl_id, convert_ensembl_to_symbol, USE.NAMES = FALSE)


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
ps_result_list <- lapply(denovo_drugCombs$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_PS))
denovo_drugCombs$ext_PS_targets <- sapply(ps_result_list, function(x) x$ext_targets)
denovo_drugCombs$ext_PS_tar_cnt <- sapply(ps_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from SIGNOR database
cat("\nExtending drug targets based on interactions from SIGNOR database\n")
signor_result_list <- lapply(denovo_drugCombs$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_SIGNOR))
denovo_drugCombs$ext_SIGNOR_targets <- sapply(signor_result_list, function(x) x$ext_targets)
denovo_drugCombs$ext_SIGNOR_tar_cnt <- sapply(signor_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from NetPath database
cat("\nExtending drug targets based on interactions from NetPath database\n")
npa_result_list <- lapply(denovo_drugCombs$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_NPA))
denovo_drugCombs$ext_NPA_targets <- sapply(npa_result_list, function(x) x$ext_targets)
denovo_drugCombs$ext_NPA_tar_cnt <- sapply(npa_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on regulatory network (TF-target interactions)
cat("\nExtending drug targets based on regulatory network (TF-target interactions)\n")
reg_result_list <- lapply(denovo_drugCombs$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_RIs))
denovo_drugCombs$ext_RI_targets <- sapply(reg_result_list, function(x) x$ext_targets)
denovo_drugCombs$ext_RI_tar_cnt <- sapply(reg_result_list, function(x) x$ext_tar_cnt)


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
kegg_result_list <- lapply(denovo_drugCombs$drugTarget_geneSymbol, function(x) extend_targets_kegg_list(x, ixns_KEGG))
denovo_drugCombs$ext_KEGG_targets <- sapply(kegg_result_list, function(x) x$ext_kegg_targets)
denovo_drugCombs$ext_KEGG_tar_cnt <- sapply(kegg_result_list, function(x) x$ext_kegg_tar_cnt)


#####


if(!dir.exists("InputFiles/DeNovo_data_1/")){dir.create("InputFiles/DeNovo_data_1/", recursive = TRUE)}
saveRDS(denovo_drugCombs, file = paste0("InputFiles/DeNovo_data_1/drugCombs_denovo1_", disease, ".rds"))

cat(paste0("\nNumber of drug combinations: ", nrow(denovo_drugCombs), "\n"))



print(warnings())