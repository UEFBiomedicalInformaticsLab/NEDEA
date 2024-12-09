set.seed(5081)


# Script to extract known drug targets from DrugBank and extend them using various databases


# Load libraries
library(optparse)
library(igraph)
library(tidyverse)
library(org.Hs.eg.db)
library(OmnipathR)
library(readxl)
source("Scripts/Functions/Functions_drug_target.R")


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

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))


#####


# Read the FIMM drug combinations
FimmDrugComb_drugCombCat <- readRDS("InputFiles/Reference_list/FimmDrugComb_drugCombinations.rds")


# Extract the disease specific drug combinations
switch(disease,
       "BreastCancer" = {tissue_select <- "breast"},
       "KidneyCancer" = {tissue_select <- "kidney"},
       "LungCancer" = {tissue_select <- "lung"},
       "OvaryCancer" = {tissue_select <- "ovary"},
       "ProstateCancer" = {tissue_select <- "prostate"},
       "SkinCancer" = {tissue_select <- "skin"})

drugCombs <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$tissue_name == tissue_select, ]
cat(paste0("\n\nNumber of drug combinations: ", nrow(drugCombs), "\n\n"))


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
cat("\n\nReading drug-target interactions")
cat(paste0("\n\tNumber of drug-target interactions: ", nrow(drug_target_ixn)))
cat(paste0("\n\tNumber of drugs: ", length(unique(drug_target_ixn$drugbank_drug_id))))
cat(paste0("\n\tNumber of targets: ", length(unique(drug_target_ixn$ensembl_gene_id))))


# Keep drug targets that are included in the input network
drug_target_ixn <- drug_target_ixn %>% 
  filter(ensembl_gene_id %in% V(input_network)$name) 
cat("\nFiltering drug-target interactions to keep targets in the network")
cat(paste0("\n\tNumber of drug-target interactions: ", nrow(drug_target_ixn)))
cat(paste0("\n\tNumber of drugs: ", length(unique(drug_target_ixn$drugbank_drug_id))))
cat(paste0("\n\tNumber of targets: ", length(unique(drug_target_ixn$ensembl_gene_id))))


# Aggregate the drug targets as a single string
drugTarget_list <- drug_target_ixn %>%
  group_by(drugbank_drug_id) %>%
  summarise(drugTarget_ensembl_id = paste(ensembl_gene_id, collapse = ","))


# Merge the drug targets information with the drug combinations data
cat("\nExtracting targets of the drug combinations\n")
drugCombs_targets <- drugCombs %>%
  left_join(drugTarget_list, by = c("Drug1_DrugBank_id" = "drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_1 = drugTarget_ensembl_id) %>%
  left_join(drugTarget_list, by = c("Drug2_DrugBank_id" = "drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_2 = drugTarget_ensembl_id) %>%
  filter(!is.na(drugTarget_ensembl_id_1), !is.na(drugTarget_ensembl_id_2)) %>%  # Remove rows where either drug has zero targets
  dplyr::rowwise() %>%
  dplyr::mutate(
    drugTarget_ensembl_id = list(unique(unlist(strsplit(c(drugTarget_ensembl_id_1, drugTarget_ensembl_id_2), ",")))),
    drugTarget_count = length(unlist(drugTarget_ensembl_id))
  ) %>%
  dplyr::select(-drugTarget_ensembl_id_1, -drugTarget_ensembl_id_2) 


# Simplify the target column as comma separated string
drugCombs_targets$drugTarget_ensembl_id <- sapply(drugCombs_targets$drugTarget_ensembl_id, function(x) paste(x, collapse = ","))


# Convert Ensembl IDs of the targets to gene symbols
drugCombs_targets$drugTarget_geneSymbol <- sapply(drugCombs_targets$drugTarget_ensembl_id, convert_ensembl_to_symbol, USE.NAMES = FALSE)
 

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
ps_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_PS))
drugCombs_targets$ext_PS_targets <- sapply(ps_result_list, function(x) x$ext_targets)
drugCombs_targets$ext_PS_tar_cnt <- sapply(ps_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from SIGNOR database
cat("\nExtending drug targets based on interactions from SIGNOR database\n")
signor_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_SIGNOR))
drugCombs_targets$ext_SIGNOR_targets <- sapply(signor_result_list, function(x) x$ext_targets)
drugCombs_targets$ext_SIGNOR_tar_cnt <- sapply(signor_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on interactions from NetPath database
cat("\nExtending drug targets based on interactions from NetPath database\n")
npa_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_NPA))
drugCombs_targets$ext_NPA_targets <- sapply(npa_result_list, function(x) x$ext_targets)
drugCombs_targets$ext_NPA_tar_cnt <- sapply(npa_result_list, function(x) x$ext_tar_cnt)


# Extend the drug targets based on regulatory network (TF-target interactions)
cat("\nExtending drug targets based on regulatory network (TF-target interactions)\n")
reg_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_ixns_database(x, ixns_RIs))
drugCombs_targets$ext_RI_targets <- sapply(reg_result_list, function(x) x$ext_targets)
drugCombs_targets$ext_RI_tar_cnt <- sapply(reg_result_list, function(x) x$ext_tar_cnt)


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
kegg_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_kegg_list(x, ixns_KEGG))
drugCombs_targets$ext_KEGG_targets <- sapply(kegg_result_list, function(x) x$ext_kegg_targets)
drugCombs_targets$ext_KEGG_tar_cnt <- sapply(kegg_result_list, function(x) x$ext_kegg_tar_cnt)


#####


# Save the basic information for the extracted drug combinations
drugCombs_1 <- drugCombs_targets[, c("Drug1_DrugBank_id", "drug_row", "Drug2_DrugBank_id", "drug_col", 
                                     "tissue_name", 
                                     "avgSynSM", "avgSynZIP", "avgSynLoewe", "avgSynHSA", "avgSynBliss", 
                                     "meanSynS", "sdSynS", "meanSynZIP", "sdSynZIP", "meanSynLoewe", 
                                     "sdSynLoewe", "meanSynHSA", "sdSynHSA", "meanSynBliss", "sdSynBliss", 
                                     "cS", "cZIP", "cLoewe", "cHSA", "cBliss", 
                                     "Syn_level", "class_synergyScore")]

if(!dir.exists("InputFiles/Drug_combination_data/")){dir.create("InputFiles/Drug_combination_data/", recursive = TRUE)}
saveRDS(drugCombs_1, file = paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))



# Save the targets for the extracted drug combinations
drugCombs_2 <- drugCombs_targets[, c("Drug1_DrugBank_id", "drug_row", "Drug2_DrugBank_id", "drug_col", 
                                     "drugTarget_ensembl_id", "drugTarget_count", "drugTarget_geneSymbol", 
                                     "ext_PS_targets", "ext_PS_tar_cnt", 
                                     "ext_SIGNOR_targets", "ext_SIGNOR_tar_cnt", 
                                     "ext_NPA_targets", "ext_NPA_tar_cnt", 
                                     "ext_RI_targets", "ext_RI_tar_cnt", 
                                     "ext_KEGG_targets", "ext_KEGG_tar_cnt")]

if(!dir.exists("InputFiles/Drug_combination_targets/")){dir.create("InputFiles/Drug_combination_targets/", recursive = TRUE)}
saveRDS(drugCombs_2, file = paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))


#####


print(warnings())