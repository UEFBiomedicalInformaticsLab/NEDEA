source("Scripts/Functions/Functions_RWR.R")




disease <- "BreastCancer"
drug_target_type <- "KEGG"
drug_comb_name <- "DB07101_DB05482"

show_efficacy <- TRUE
efficacy_lib_select <-"Breast adenocarcinoma (C0858252)[DisGeNET_curated]"

show_safety <- FALSE
safety_lib_select <- "Adverse reaction (08.06.01.018) [ADReCS]"

show_RWR_genes <- FALSE


# Read the network on which to execute RWR
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))




# Read the libraries and extract gene list
if(show_efficacy){
  enrichment_lib_efficacy <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  enrichment_lib_efficacy <- lapply(enrichment_lib_efficacy, function(x){x[x %in% V(input_network)$name]})
  enrichment_lib_efficacy <- enrichment_lib_efficacy[[efficacy_lib_select]]
}


if(show_safety){
  enrichment_lib_safety <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
  enrichment_lib_safety <- lapply(enrichment_lib_safety, function(x){x[x %in% V(input_network)$name]})
  enrichment_lib_safety <- enrichment_lib_safety[[safety_lib_select]]
}







# Extract the drug targets 
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, 
                                     drugCombs_targets$Drug2_DrugBank_id, 
                                     sep = "_")
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

target_set <- drugCombs_targets[drugCombs_targets$comb_name %in% drug_comb_name, drug_target_col, drop = TRUE]
target_set <- unlist(strsplit(target_set, ","))
suppressMessages(target_set <- select(org.Hs.eg.db, 
                                      keys = target_set, 
                                      columns = "ENSEMBL", 
                                      keytype = "SYMBOL"))
target_set <- target_set$ENSEMBL

target_set <- target_set[target_set %in% V(input_network)$name]












# Extract the genes selected for FGSEA

if(show_RWR_genes){
  rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
  if(!file.exists(rwr_result_file_path)){
    stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
  }
  rwr_data <- readRDS(file = rwr_result_file_path)
  rwr_data_select <- rwr_data[, drug_comb_name]
  rwr_threshold <- func_RWR_threshold(rwr_data_select[!names(rwr_data_select) %in% target_set])
  rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold$ELB], decreasing = TRUE)
}



# Prepare the node labels
node_labels <- as.data.frame(matrix(NA, 
                                    nrow = vcount(input_network), 
                                    ncol = 5, 
                                    dimnames = list(c(V(input_network)$name), 
                                                    c("Efficacy", "Safety", "Targets", "RWR_genes", "RWR_probs")
                                    )
))

View(node_labels)


node_labels[target_set, ]$Targets <- "target"

if(show_efficacy){node_labels[enrichment_lib_efficacy, ]$Efficacy <- "efficacy"}
if(show_safety){node_labels[enrichment_lib_safety, ]$Safety <- "safety"}
if(show_RWR_genes){
  node_labels[names(rankedGeneList), ]$RWR_genes <- "rwrGenes"
  rwr_data_select <- as.data.frame(rwr_data_select)
  node_labels$RWR_probs <- rwr_data_select$rwr_data_select[match(row.names(node_labels), row.names(rwr_data_select))]
}

node_labels <- node_labels %>% 
  unite(col = label, Efficacy, Safety, Targets, RWR_genes, 
        sep = "_", remove = FALSE, na.rm = TRUE)



if(nrow(node_labels[grep("_", node_labels$label), ]) > 0 ){
  node_labels[grep("_", node_labels$label), ]$label <- "overlap"
}
node_labels[node_labels$label == "", ]$label <- "others"


# Add color to nodes 


# node_labels$color <- switch(node_labels$label, 
#                             "efficacy" = "green", 
#                             "safety" = "red",
#                             "overlap" = "orange",
#                             "rwrGenes" = "blue")



node_labels$color <- unlist(sapply(node_labels$label, 
                                   FUN = function(x){switch(x, 
                                                            "efficacy" = "green", 
                                                            "safety" = "red",
                                                            "overlap" = "orange",
                                                            "rwrGenes" = "cyan", 
                                                            "target" = "blue",
                                                            "others" = "lightblue"
                                   )}, 
                                   USE.NAMES = FALSE))

# Map the node attributes to the graph
for(attr in colnames(node_labels)){
  input_network <- set_vertex_attr(graph = input_network,
                                   name = attr,
                                   index = V(input_network)[V(input_network)$name %in% row.names(node_labels)],
                                   value = node_labels[,attr])
}




select_nodes <- node_labels[node_labels$label != "others", ]
subnet <- induced_subgraph(graph = input_network, 
                           vids = V(input_network)[V(input_network)$name %in% row.names(select_nodes)])





visIgraph(subnet, 
          physics = FALSE, 
          smooth = FALSE, 
          randomSeed = 5081) %>%
  visIgraphLayout(layout = "layout_nicely") %>% 
  visOptions(width = 1920, 
             height = 1200,
             highlightNearest = list("hover" = TRUE),
             clickToUse = TRUE)







