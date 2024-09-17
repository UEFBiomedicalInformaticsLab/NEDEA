set.seed(5081)


# Build the input files for Drug2Ways


# Load libraries
library(tidyverse)


#####

disease <- "SkinCancer"
drug_target_type <- "known"

# Prepare the gene-gene links

input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
input_network <- igraph::as_data_frame(input_network, "edges")
colnames(input_network) <- c("source", "target")
input_network$polarity <- 1


#####


# Prepare phenotype-gene link

enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))

enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))

enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

# Keep only features that were used in the predictive system
model <- readRDS(file = paste0("OutputFiles/Predictive_model/model_NES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
predictive_system_features <- model$feature_threshold$feature
enrichment_lib <- enrichment_lib[names(enrichment_lib) %in% predictive_system_features] 

phenotype_gene_link <- bind_rows(lapply(enrichment_lib, as.data.frame), 
                                 .id = "source") 
colnames(phenotype_gene_link)[2] <- "target"

phenotype_gene_link$polarity <- 1


#####


# Prepare the drug-target associations

drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "comb_name", "class_EffAdv")]

drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
drugCombs_targets <- drugCombs_targets[drugCombs_targets$comb_name %in% drugCombs_cat$comb_name, ]

drugCombs_targets <- drugCombs_targets %>% 
  select(c(comb_name, drugTarget_ensembl_id)) %>% 
  separate_rows(drugTarget_ensembl_id, sep = ",") %>%
  rename("source" = "comb_name", "target" = "drugTarget_ensembl_id")

drugCombs_targets$polarity <- 1


#####

if(!dir.exists("InputFiles/Compare_external/drug2ways")){
  dir.create("InputFiles/Compare_external/drug2ways", recursive = TRUE)
}

# Export the network
network <- rbind(input_network, phenotype_gene_link, drugCombs_targets)
write.table(input_network, file = "InputFiles/Compare_external/drug2ways/drug2ways_network.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Export the drug list
drugCombs <- data.frame("comb_name" = unique(drugCombs_targets$source))
write.table(drugCombs, file = "InputFiles/Compare_external/drug2ways/drug2ways_drugs.tsv", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


# Export phenotype list
phenotypes <- data.frame("phenotype" = unique(phenotype_gene_link$source))
write.table(phenotypes, file = "InputFiles/Compare_external/drug2ways/drug2ways_phenotypes.tsv", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


#####


print(warnings())