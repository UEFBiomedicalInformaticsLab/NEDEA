library(openxlsx)
library(igraph)
library(org.Hs.eg.db)
library(dnet)
library(foreach)
library(doParallel)
library(tidyverse)
source("Scripts/Functions/Functions_drug_target.R")
source("Scripts/Functions/Functions_parallelprocesses.R")
source("Scripts/Functions/Functions_FGSEA.R")


disease <- "LungCancer"
drug_target_type <- "known"



drugCombs <- read.xlsx("Barabasi_test/41467_2019_9186_MOESM5_ESM.xlsx", sheet = 4)
drugCombs <- drugCombs[1:73,1:2]
Bbsi_drug_list <- unique(c(drugCombs$`Interaction_A(DrugBank_ID)`, drugCombs$`Interaction_B(DrugBank_ID)`))

colnames(drugCombs) <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")


# for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
#   print(disease)
#   drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
#   my_drug_list <- unique(c(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id))
# 
#   print(length(intersect(my_drug_list, Bbsi_drug_list)))
#   print(drugCombs_targets[drugCombs_targets$drugTarget_count < 2,])
# }





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




drug_target_col <- "drugTarget_geneSymbol"

# Extract the targets for the drug combinations
seed_matrix_targets <- drugCombs_targets[,drug_target_col, drop = FALSE]
seed_matrix_targets <- apply(seed_matrix_targets, 1, function(x){paste(x, collapse = ",")})




# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))



# Convert the target sets to a matrix, each column corresponds to a drug pair
rwr_seed_matrix <- do.call(cbind, lapply(seed_matrix_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys = target_set, 
                                     columns = "ENSEMBL", 
                                     keytype = "SYMBOL"))
  sapply(V(rwr_input_network)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

colnames(rwr_seed_matrix) <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")






# Run RWR using dRWR
nproc <- 10
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 

rwr_result <- dRWR(g = rwr_input_network, 
                   setSeeds = rwr_seed_matrix,
                   normalise = "row",
                   restart = 0.5, 
                   normalise.affinity.matrix = "none",
                   multicores = nproc)

colnames(rwr_result) <- colnames(rwr_seed_matrix)
rownames(rwr_result) <- rownames(rwr_seed_matrix)


stopCluster(cl)
unregister_dopar()





fgsea_result_final <- list()
# FGSEA on combined efficacy-safety library
cat("\n--- Executing FGSEA on combined efficacy-safety library")
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries_extended/Disease2Gene_", disease, "_extendedLib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries_extended/curatedAdr2Gene_extendedLib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

# fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
#                                       enrichment_library = enrichment_lib,
#                                       disease = disease,
#                                       drug_target_type = drug_target_type,
#                                       nproc = nproc)




rwr_data = rwr_result
enrichment_library = enrichment_lib





# Create an empty matrix to store the enrichment result
enrichment_result_mat <- matrix(0, 
                                nrow = length(enrichment_library), 
                                ncol = ncol(rwr_data), 
                                dimnames = list(names(enrichment_library), 
                                                colnames(rwr_data)))





drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, 
                                     drugCombs_targets$Drug2_DrugBank_id, 
                                     sep = "_")



# For each drug combination run FGSEA and extract the significant NES
for(drugComb in colnames(rwr_data)){
  
  # Extract the targets of the drug combination
  target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drugComb, drug_target_col, drop = TRUE]
  if(drug_target_type == "all"){
    target_set <- apply(target_set, 1, function(x){paste(x, collapse = ",")})
  }
  target_set <- unique(unlist(strsplit(target_set, ",")))
  suppressMessages(target_set <- select(org.Hs.eg.db, 
                                        keys = target_set, 
                                        columns = "ENSEMBL", 
                                        keytype = "SYMBOL"))
  target_set <- target_set$ENSEMBL
  
  
  # Select the genes for FGSEA
  rwr_data_select <- rwr_data[, drugComb]
  rwr_threshold <- func_RWR_threshold(rwr_data_select[!names(rwr_data_select) %in% target_set])
  rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold$ELB], decreasing = TRUE)
  
  
  # Run FGSEA
  skip_iteration <- FALSE  # Initialize flag variable
  
  tryCatch({
    
    enrichment_result <- fgseaMultilevel(pathways = enrichment_library,
                                         stats = rankedGeneList,
                                         minSize = 5, 
                                         maxSize = 1500, 
                                         scoreType = "pos", 
                                         nproc = nproc,
                                         BPPARAM = MulticoreParam(progressbar = FALSE))
    
    enrichment_result$NES[which(enrichment_result$padj > 0.05)] <- 0
    
    enrichment_result_mat[enrichment_result$pathway, drugComb] <- enrichment_result$NES
    
  }, 
  error = function(e) {
    if (grepl("GSEA statistic is not defined when all genes are selected", conditionMessage(e))) {
      cat(paste("\n\t- Skipping iteration for", drugComb, "due to GSEA statistic error."))
      skip_iteration <- TRUE  # Skip to next iteration of the loop
    } else {
      # Handle other errors here if needed, or re-throw the error
      stop(e)
    }
  })
  if (skip_iteration) {
    next  # Skip to the next iteration of the loop
  }
}

fgsea_result <- enrichment_result_mat





fgsea_result_final[["combinedEfficacySafety"]] <- fgsea_result


tmp1 <- colSums(fgsea_result)
length(tmp1)
length(tmp1[tmp1>0])







plot_data <- fgsea_result


# Find the features to plot

feature_vars <- apply(plot_data, 1, var)
feature_vars <- sort(feature_vars, decreasing = TRUE)
efficacy_feature_select <- names(feature_vars[grep("^\\[DISEASE\\]", names(feature_vars))][1])
safety_feature_select <- names(feature_vars[grep("^\\[ADR\\]", names(feature_vars))][1])

plot_data <- plot_data[row.names(plot_data) %in% c(efficacy_feature_select, safety_feature_select), ,  drop = FALSE]

x_axis_label = str_wrap(efficacy_feature_select, 30)
y_axis_label = str_wrap(safety_feature_select, 30)












plot_data <- as.data.frame(t(plot_data))
colnames(plot_data) <- c("F1", "F2")

ggplot() +
  geom_point(data = plot_data, 
             mapping = aes(x = F1, y = F2),  
             size = 0.5,
             shape = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  labs(title = paste0("Disease: ", disease,
                      "\nTarget: ", drug_target_type),
       x = x_axis_label,
       y = y_axis_label,
       color = "Class")







