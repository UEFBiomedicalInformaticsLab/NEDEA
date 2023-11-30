

library(tidyverse)
library(ComplexHeatmap)
library(igraph)
library(org.Hs.eg.db)
source("Scripts/Functions/Functions_RWR.R")



disease <- "BreastCancer"
drug_target_type <- "known"
drug_comb_name <- "DB01204_DB00619"
lib_name <-"[Efficacy_LungCancer] Adenocarcinoma of lung (disorder) (C0152013)[DisGeNET_curated]"


# Read the network on which to execute RWR
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))



# library(RCy3)
# createNetworkFromIgraph(input_network, title = "input_network", collection = "Network")

# # Generate the adjacency matrix
# network_adj_mat <- as_adjacency_matrix(graph = input_network, 
#                                        type = "both", 
#                                        sparse = FALSE)
# 
# 
# 






# Read the libraries and extract gene list
enrichment_lib_efficacy <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
enrichment_lib_efficacy <- lapply(enrichment_lib_efficacy, function(x){x[x %in% V(input_network)$name]})
enrichment_lib_efficacy <- unique(unlist(enrichment_lib_efficacy, use.names = FALSE))

enrichment_lib_safety <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
enrichment_lib_safety <- lapply(enrichment_lib_safety, function(x){x[x %in% V(input_network)$name]})
enrichment_lib_safety <- unique(unlist(enrichment_lib_safety, use.names = FALSE))



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



# Extract the genes selected for FGSEA
rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
if(!file.exists(rwr_result_file_path)){
  stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
}
rwr_result <- readRDS(file = rwr_result_file_path)
rwr_threshold <- sapply(apply(rwr_result, 2, func_RWR_threshold), function(x){x$ELB})
rwr_data_select <- rwr_result[, drug_comb_name]
rwr_threshold <- func_RWR_threshold(probabilities = rwr_data_select)
rwr_threshold <- rwr_threshold$ELB
rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold], decreasing = TRUE)









# Prepare the node labels
node_labels <- as.data.frame(matrix(NA, 
                                    nrow = vcount(input_network), 
                                    ncol = 5, 
                                    dimnames = list(c(V(input_network)$name), 
                                                    c("Efficacy", "Safety", "Targets", "rankedGenes", "RWR_probs")
                                    )
))

View(node_labels)

node_labels[enrichment_lib_efficacy, ]$Efficacy <- "efficacy"
node_labels[enrichment_lib_safety, ]$Safety <- "safety"
node_labels[target_set, ]$Targets <- "target"
node_labels[names(rankedGeneList), ]$rankedGenes <- "rankedGene"

rwr_data_select <- as.data.frame(rwr_data_select)
node_labels$RWR_probs <- rwr_data_select$rwr_data_select[match(row.names(node_labels), row.names(rwr_data_select))]

node_labels <- node_labels %>% unite(col = label, Efficacy, Safety, Targets, rankedGenes, sep = "_", remove = FALSE, na.rm = TRUE)





# Map the node attributes to the graph
for(attr in colnames(node_labels)){
  input_network <- set_vertex_attr(graph = input_network,
                                   name = attr,
                                   index = V(input_network)[V(input_network)$name %in% row.names(node_labels)],
                                   value = node_labels[,attr])
}








select_nodes <- node_labels[node_labels$label != "", ]
plot_net <- induced_subgraph(graph = input_network, vids = V(input_network)[V(input_network)$name %in% row.names(select_nodes)])





visIgraph(plot_net, 
          physics = FALSE, 
          smooth = FALSE, 
          randomSeed = 5081) %>%
  visIgraphLayout(layout = "layout_nicely") %>% 
  visOptions(width = 1920, 
             height = 1200,
             highlightNearest = list("hover" = TRUE),
             clickToUse = TRUE)





rglplot(plot_net)










library(RCy3)
createNetworkFromIgraph(plot_net, title = "Plot Network", collection = "Network")
loadTableData(data = node_labels)
# setNodeSizeBypass(node.names = V(plot_net)$name, new.sizes = 35)
# layoutNetwork(layout.name='circular')
style.name = "dataStyle"
createVisualStyle(style.name)
setVisualStyle(style.name)

setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(1, style.name)
setNodeLabelMapping('display name', style.name)




# ggraph(input_network, 
#        layout = "igraph", 
#        algorithm = "fr") +
#   geom_edge_link() + 
#   geom_node_point() + 
#   theme_void() +
#   theme(legend.position = "none") 


dist_mat <- dist(network_adj_mat)
tsne_res <- Rtsne(dist_mat, is_distance = TRUE)


dist_mat_1 <- distances(graph = input_network, 
                        v = V(input_network), 
                        to = V(input_network), 
                        mode = "all", 
                        algorithm = "dijkstra")




tsne_res_1 <- Rtsne(dist_mat_1, is_distance = TRUE)



plot_data <- as.data.frame(tsne_res_1$Y)
row.names(plot_data) <- row.names(dist_mat_1)
colnames(plot_data) <- c("tsne_dim1", "tsne_dim2")


plot_data <- merge(plot_data , node_labels, by = 0)
plot_data <- column_to_rownames(plot_data, "Row.names")
plot_data[plot_data$label == "", ]$label <- "others"
plot_data[grep("_", plot_data$label), ]$label <- "overlap"


plot_data$label <- as.factor(plot_data$label)

point_colour <- c("efficacy" = "green",
                  "safety" = "red",
                  "target" = "blue", 
                  "rankedGene" = "cyan", 
                  "others" = "black",
                  "overlap" = "orange")


plot <- ggplot(data = plot_data, 
               aes(x = tsne_dim1, y = tsne_dim2, color = label, order = RWR_probs)) + 
  geom_point(shape = 3) + 
  scale_color_manual(values = point_colour)





library(plotly)
ggplotly(plot)




hclust_res <- hclust(d = as.dist(dist_mat_1))















