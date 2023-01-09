# ppi <- expand.grid(paste0("Protein", seq(1:10)), paste0("Protein", seq(1:10))) # For one component network
# colnames(ppi) <- c("Node1", "Node2")
# ppi <- ppi[ppi$Node1 != ppi$Node2, ]

# ppi <- ppi[sample(1:nrow(ppi), 30), ] # Select only 30% edges

# ppi <- graph_from_data_frame(ppi, directed = FALSE)
# 
# ppi <- delete_edges(ppi, E(ppi))

# Network with two components
ppi1 <- expand.grid(paste0("Protein", 1:10), paste0("Protein", 1:10))
ppi2 <- expand.grid(paste0("Protein", 11:20), paste0("Protein", 11:20))
ppi <- rbind(ppi1, ppi2)
colnames(ppi) <- c("Node1", "Node2")
ppi <- ppi[ppi$Node1 != ppi$Node2, ]
ppi <- ppi[sample(1:nrow(ppi), 30), ] # Select only 30% edges
ppi <- graph_from_data_frame(ppi, directed = FALSE)
edge_density(ppi)
count_components(ppi)


dti <- expand.grid(paste0("Drug", seq(1:4)), paste0("Protein", seq(1:10)))
colnames(dti) <- c("Node1_drugbank_drug_id", "Node2_ensembl_gene_id")
dti <- dti[sample(1:nrow(dti), 5), ]
dti


drug1 <- "Drug1"
drug2 <- "Drug4"

rwr_input_seed <- func_RWR_seed_from_DTI(rwr_input_network = ppi, 
                                         drugs = c(drug1, drug2), 
                                         drug_target_ixn = dti)


union_seed <- ifelse(rowSums(rwr_input_seed) > 0, 1, 0)
intersect_seed <- ifelse(rowSums(rwr_input_seed) == 2, 1, 0)
rwr_input_seed <- cbind(rwr_input_seed, union_seed, intersect_seed)



# Execute dNet RWR
rwr_result <- suppressMessages(dRWR(g = ppi, 
                                    setSeeds = rwr_input_seed, 
                                    normalise = "none", 
                                    restart = 0.5, 
                                    normalise.affinity.matrix = "none"))


# Prepare result for export
rwr_result <- as.data.frame(as.matrix(rwr_result))
colnames(rwr_result) <- colnames(rwr_input_seed)
rwr_result$node_name <- V(ppi)$name
rwr_result$is_union_seed <- as.logical(as.data.frame(rwr_input_seed)$union_seed[match(rwr_result$node_name, row.names(rwr_input_seed))])

rwr_result$addEffect <- rowSums(rwr_result[,c(drug1, drug2)])
rwr_result <- rwr_result[order(rwr_result$union_seed, decreasing = TRUE), ]
rwr_result <- rwr_result[, c("node_name", "is_union_seed",drug1, drug2, "addEffect", "union_seed", "intersect_seed")]
row.names(rwr_result) <- NULL
rwr_result