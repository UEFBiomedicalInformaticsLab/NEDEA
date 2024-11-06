set.seed(5081)



# Generate the features for validation data 6a and check applicability domain



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(igraph)
library(org.Hs.eg.db)
library(dnet)
library(foreach)
library(doParallel)
library(BiocParallel)
library(fgsea)
source("Scripts/Functions/Functions_RWR.R")
source("Scripts/Functions/Functions_parallelprocesses.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known",
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL,
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer.", call.=FALSE)
  }
}


# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type
nproc <- opt$nproc


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type))
cat(paste0("\nnproc: ", nproc, "\n\n"))


#####


# Read the network on which to execute RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(rwr_input_network), ", edges = ", ecount(rwr_input_network), "\n\n"))


# Read the drug combinations
valid_drugCombs_cat <- readRDS(file = paste0("InputFiles/Validation_data_6a/drugCombs_validation6a_", disease, ".rds"))


if(nrow(valid_drugCombs_cat) == 0){
  stop(paste0("No drug combinations for ", disease, ". No features will be generated."), call.=FALSE)
}



# Select the column containing the drug target based on th user input
switch(drug_target_type,
       "known" = { drug_target_col <- c("drugTarget_geneSymbol") },
       "PS" = { drug_target_col <- c("ext_PS_targets") },
       "SIGNOR" = { drug_target_col <- c("ext_SIGNOR_targets") },
       "NPA" = { drug_target_col <- c("ext_NPA_targets") },
       "RI" = { drug_target_col <- c("ext_RI_targets") },
       "KEGG" = { drug_target_col <- c("ext_KEGG_targets") },
       "all" = { drug_target_col <- c("drugTarget_geneSymbol", "ext_PS_targets", 
                                      "ext_SIGNOR_targets", "ext_NPA_targets", 
                                      "ext_RI_targets", "ext_KEGG_targets") })

valid_drugCombs_cat <- valid_drugCombs_cat[, c("comb_name", 
                                               "Drug1_DrugBank_id", "Drug2_DrugBank_id", 
                                               "Drug1_name", "Drug2_name", "NCT_ID", 
                                               "class_EffAdv", drug_target_col)]


#####


# Extract the targets for the drug combinations
seed_matrix_targets <- valid_drugCombs_cat[,drug_target_col, drop = FALSE]
seed_matrix_targets <- apply(seed_matrix_targets, 1, function(x){paste(x, collapse = ",")})


# Convert the target sets to a matrix, each column corresponds to a drug pair
rwr_seed_matrix <- do.call(cbind, lapply(seed_matrix_targets, function(x) {
  target_set <- unlist(strsplit(x, ","))
  
  suppressMessages(mapping <- select(org.Hs.eg.db, 
                                     keys = target_set, 
                                     columns = "ENSEMBL", 
                                     keytype = "SYMBOL"))
  sapply(V(rwr_input_network)$name, function(y) ifelse(y %in% mapping$ENSEMBL, 1, 0))
}))

colnames(rwr_seed_matrix) <- valid_drugCombs_cat$comb_name


# Run RWR using dRWR
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


#####


# FGSEA on combined efficacy-safety library
cat("\n--- Executing FGSEA on combined efficacy-safety library")
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

# Create an empty matrix to store the enrichment result
enrichment_result_mat <- matrix(0, 
                                nrow = length(enrichment_lib), 
                                ncol = ncol(rwr_result), 
                                dimnames = list(names(enrichment_lib), 
                                                colnames(rwr_result)))

# For each drug combination run FGSEA and extract the significant NES
for(drugComb in colnames(rwr_result)){
  
  # Extract the targets of the drug combination
  target_set <- valid_drugCombs_cat[valid_drugCombs_cat$comb_name %in% drugComb, drug_target_col, drop = TRUE]
  target_set <- unique(unlist(strsplit(target_set, ",")))
  suppressMessages(target_set <- select(org.Hs.eg.db, 
                                        keys = target_set, 
                                        columns = "ENSEMBL", 
                                        keytype = "SYMBOL"))
  target_set <- target_set$ENSEMBL
  
  # Select the genes for FGSEA
  rwr_data_select <- rwr_result[, drugComb]
  rwr_threshold <- func_RWR_threshold(rwr_data_select[!names(rwr_data_select) %in% target_set])
  rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold$ELB], decreasing = TRUE)
  
  # Run FGSEA
  skip_iteration <- FALSE  # Initialize flag variable
  
  tryCatch({
    
    enrichment_result <- fgseaMultilevel(pathways = enrichment_lib,
                                         stats = rankedGeneList,
                                         minSize = 5, 
                                         maxSize = 500, 
                                         scoreType = "pos", 
                                         nproc = nproc,
                                         BPPARAM = MulticoreParam(progressbar = FALSE))
    
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


#####


# Save the results

if(!dir.exists("OutputFiles/Validation_data_6a/Features/")){dir.create("OutputFiles/Validation_data_6a/Features/", recursive = TRUE)}
saveRDS(enrichment_result_mat, file = paste0("OutputFiles/Validation_data_6a/Features/fgseaNES_combinedEfficacySafety_", disease, "_", drug_target_type, ".rds"))


#####


# Read the training set drug combination 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]


# Read the FGSEA results
train_fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
train_fgsea_result <- train_fgsea_result[["combinedEfficacySafety"]]

valid_fgsea_result <- enrichment_result_mat


# Read the important features and select the top to plot
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"))

safety_feature_select <- selected_features[grep("^\\[ADR\\]", selected_features$feature), ]
if(nrow(safety_feature_select) > 0){
  safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)][1]
}else{ safety_feature_select <- c() }

efficacy_feature_select <- selected_features[grep("^\\[DISEASE\\]", selected_features$feature), ]
if(nrow(efficacy_feature_select) > 0){
  efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)][1]
}else{ efficacy_feature_select <- c() }


# Get the FGSEA results for the train and test data
train_fgsea_result <- train_fgsea_result[row.names(train_fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                         colnames(train_fgsea_result) %in% train_drugCombs_cat$comb_name, drop = FALSE]

valid_fgsea_result <- valid_fgsea_result[row.names(valid_fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                         colnames(valid_fgsea_result) %in% valid_drugCombs_cat$comb_name, drop = FALSE]


train_fgsea_result <- as.data.frame(t(train_fgsea_result))
train_fgsea_result$category <- train_drugCombs_cat$class_EffAdv[match(row.names(train_fgsea_result), train_drugCombs_cat$comb_name)]
train_fgsea_result$data_group <- "Training"

valid_fgsea_result <- as.data.frame(t(valid_fgsea_result))
valid_fgsea_result$category <- valid_drugCombs_cat$class_EffAdv[match(row.names(valid_fgsea_result), valid_drugCombs_cat$comb_name)]
valid_fgsea_result$data_group <- "Validation"


# Plot

plot_data <- rbind(train_fgsea_result, valid_fgsea_result)

colnames(plot_data)[colnames(plot_data) %in% efficacy_feature_select] <- "F1"
colnames(plot_data)[colnames(plot_data) %in% safety_feature_select] <- "F2"

x_axis_label = str_wrap(efficacy_feature_select, 40)
y_axis_label = str_wrap(safety_feature_select, 40)



if(!dir.exists("OutputFiles/Plots/validation_data_applicibility_domain/")){
  dir.create("OutputFiles/Plots/validation_data_applicibility_domain/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/validation_data_applicibility_domain/plot_validation6a_AD_combinedEfficacySafety_", disease, "_", drug_target_type, ".tiff"),
     width = 7, height = 6,
     units = "cm", compression = "lzw", res = 1200)


ggplot() +
  geom_point(data = plot_data, 
             mapping = aes(x = F1, y = F2, color = category, shape = data_group),  
             size = 0.5, 
             stroke = 0.1) +
  scale_shape_manual(values = c("Training" = 1, "Validation" = 3)) + 
  scale_color_manual(values = c("Eff" = "#77DD77", "Adv" = "#FF6961")) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "right",
        legend.box.spacing = unit(0.1, "cm"),
        legend.key = element_rect(fill = NA, linewidth = 0), 
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 2.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        legend.text = element_text(size = 2, margin = margin(t = 0, r = 0, b = 0, l = 1)),
        legend.margin = margin(1,0.5,0.5,0.5),
        legend.spacing.y = unit(0.15, "cm")  ) +
  guides(color = guide_legend(override.aes = list(size = 1)),
         shape = guide_legend(override.aes = list(size = 1))) +
  labs(title = disease,
       x = x_axis_label,
       y = y_axis_label,
       color = "Category", 
       shape = "Group")

dev.off()



print(warnings())