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




# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}


# Define global options for this script 
disease <- opt$disease

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))


# Read the drug combinations
drugCombs <- readRDS("InputFiles/Drug_combinations/Drug_combinations.rds")
drugCombs <- drugCombs[[disease]]
cat(paste0("\n\nNumber of drug combinations: ", nrow(drugCombs), "\n\n"))


# Read the network
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))


# Read drug target interactions
drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")
cat("\n\nReading drug-target interactions")
cat(paste0("\n\tNumber of drug-target interactions: ", nrow(drug_target_ixn)))
cat(paste0("\n\tNumber of drugs: ", length(unique(drug_target_ixn$Node1_drugbank_drug_id))))
cat(paste0("\n\tNumber of targets: ", length(unique(drug_target_ixn$Node2_ensembl_gene_id))))


# Keep drug targets that are included in the input network
drug_target_ixn <- drug_target_ixn %>% 
  filter(Node2_ensembl_gene_id %in% V(input_network)$name) 
cat("\nFiltering drug-target interactions to keep targets in the network")
cat(paste0("\n\tNumber of drug-target interactions: ", nrow(drug_target_ixn)))
cat(paste0("\n\tNumber of drugs: ", length(unique(drug_target_ixn$Node1_drugbank_drug_id))))
cat(paste0("\n\tNumber of targets: ", length(unique(drug_target_ixn$Node2_ensembl_gene_id))))


# Aggregate the drug targets as a single string
drugTarget_list <- drug_target_ixn %>%
  group_by(Node1_drugbank_drug_id) %>%
  summarise(drugTarget_ensembl_id = paste(Node2_ensembl_gene_id, collapse = ","))


# Merge the drug targets information with the drug combinations data
cat("\nExtracting targets of the drug combinations\n")
drugCombs_targets <- drugCombs %>%
  left_join(drugTarget_list, by = c("Drug1_DrugBank_id" = "Node1_drugbank_drug_id")) %>%
  dplyr::rename(drugTarget_ensembl_id_1 = drugTarget_ensembl_id) %>%
  left_join(drugTarget_list, by = c("Drug2_DrugBank_id" = "Node1_drugbank_drug_id")) %>%
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


# Extract known gene/protein interactions from various databases 
ixns_PS <- import_omnipath_interactions(resources = c("PhosphoSite","phosphoELM","PhosphoNetworks","PhosphoSite_MIMP",
                                                      "PhosphoSite_ProtMapper","phosphoELM_MIMP","PhosphoPoint",
                                                      "PhosphoSite_KEA","PhosphoSite_noref"))
ixns_SIGNOR <- import_omnipath_interactions(resources = c("SIGNOR","SIGNOR_ProtMapper","SIGNOR_CollecTRI"))
ixns_NPA <- import_omnipath_interactions(resources = c("NetPath"))
ixns_RIs <- import_transcriptional_interactions()


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


# # Extract the gene symbols for all the known targets in drug-target interactions
# all_drug_target_gene_symbols <- select(org.Hs.eg.db, 
#                                        keys = unique(drug_target_ixn$Node2_ensembl_gene_id), 
#                                        columns = "SYMBOL", 
#                                        keytype = "ENSEMBL")

# Download protein interactions from KEGG
ixns_KEGG <- download_KEGG_ixns(pathway_ids = CHG_pathways)


# Extend the drug targets using interactions from KEGG 
cat("\nExtending drug targets using KEGG\n")
kegg_result_list <- lapply(drugCombs_targets$drugTarget_geneSymbol, function(x) extend_targets_kegg_list(x, ixns_KEGG))
drugCombs_targets$ext_KEGG_targets <- sapply(kegg_result_list, function(x) x$ext_kegg_targets)
drugCombs_targets$ext_KEGG_tar_cnt <- sapply(kegg_result_list, function(x) x$ext_kegg_tar_cnt)



if(!dir.exists("InputFiles/Drug_targets/")){dir.create("InputFiles/Drug_targets/", recursive = TRUE)}
saveRDS(drugCombs_targets, file = paste0("InputFiles/Drug_targets/Drug_targets_extended_", disease, ".rds"))



print(warnings())